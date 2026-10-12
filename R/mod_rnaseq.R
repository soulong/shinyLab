# shinyLab - RNA-seq Module
# RNA-seq differential expression pipeline (based on shinyCat)
# Workflow: upload counts + metadata -> QC & DESeq2 -> Wald comparisons -> ORA/GSEA enrichment

library(shiny)
library(shinydashboard)
library(tidyverse)
library(readxl)
library(rhandsontable)
library(DT)
library(shinyjs)
library(DESeq2)
library(corrplot)
library(ggforce)
library(ggrepel)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(furrr)


# =============================================================================
# Constants
# =============================================================================

rnaseq_file_accept <- c(".csv", ".xlsx")

# metadata columns reserved for bookkeeping; they are never offered as design
# factors (matched by name, so column order in the uploaded file is free)
rnaseq_metadata_reserved_cols <- c("sample_name", "batch", "include")

# worker pool for parallel enrichment (multisession), set up once at startup
rnaseq_n_workers <- function() {
  cores <- parallel::detectCores()
  max(1, floor(ifelse(is.na(cores), 2, cores) / 2))
}
future::plan(future::multisession, workers = rnaseq_n_workers())


# =============================================================================
# Pure Computation Functions
# =============================================================================

#' Build a DESeq2 object from raw counts and colData
#'
#' @param counts data.frame with gene_id/gene_name columns followed by samples
#' @param coldata data.frame (row names = sample_name) with batch and condition
#' @param filter_symbol logical, remove rows with empty gene_name
#' @param min_counts numeric, keep genes with rowMeans > min_counts
#' @return a DESeq2 DESeqDataSet object
rnaseq_build_dds <- function(counts, coldata, filter_symbol = TRUE, min_counts = 10) {
  if (filter_symbol) {
    counts <- dplyr::filter(counts, !is.na(gene_name))
  } else if ("gene_id" %in% colnames(counts)) {
    counts <- dplyr::mutate(counts, gene_name = ifelse(is.na(gene_name), gene_id, gene_name))
  }

  counts <- counts[, c("gene_name", rownames(coldata))] %>%
    dplyr::distinct(gene_name, .keep_all = TRUE) %>%
    tibble::column_to_rownames(var = "gene_name")

  counts <- counts[rowMeans(counts) > min_counts, ]

  if (length(unique(coldata$batch)) > 1) {
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ batch + condition)
  } else {
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
  }

  DESeq(dds, test = "Wald", quiet = TRUE)
}


#' Wald test (with ashr shrinkage) for one comparison
#'
#' @param dds DESeqDataSet from rnaseq_build_dds
#' @param contrast character, contrast condition
#' @param control character, control condition
#' @return tibble with gene_name, baseMean, stat, pvalue, padj, shrunk log2FoldChange/lfcSE
rnaseq_wald_result <- function(dds, contrast, control) {
  compare_condition <- c("condition", contrast, control)
  result_unshrink <- DESeq2::results(dds, contrast = compare_condition)
  result_shrink <- lfcShrink(dds, contrast = compare_condition, res = result_unshrink,
                             type = "ashr", quiet = TRUE)

  res <- as_tibble(result_unshrink, rownames = "gene_name") %>%
    dplyr::select(gene_name, baseMean, stat, pvalue, padj) %>%
    left_join(as_tibble(result_shrink, rownames = "gene_name") %>%
                dplyr::select(gene_name, log2FoldChange, lfcSE),
              by = "gene_name") %>%
    relocate(log2FoldChange, .before = pvalue) %>%
    arrange(pvalue) %>%
    mutate(baseMean = round(baseMean, 0),
           log2FoldChange = round(log2FoldChange, 5),
           lfcSE = round(lfcSE, 5),
           stat = round(stat, 5))
  res
}


#' Build ranked gene lists (by Wald stat) for GSEA
#'
#' @param results_list named list of wald result tibbles (per comparison)
#' @return named list of named numeric vectors (gene_name -> stat), sorted desc
rnaseq_gsea_ranked <- function(results_list) {
  purrr::map(results_list, function(x) {
    x %>%
      dplyr::distinct(gene_name, .keep_all = TRUE) %>%
      dplyr::arrange(dplyr::desc(stat)) %>%
      dplyr::select(gene_name, stat) %>%
      tibble::deframe()
  })
}


# memo cache: full MSigDB table per species, avoids re-downloading on database switches
rnaseq_msigdb_cache <- new.env()


#' Build TERM2GENE table from MSigDB
#'
#' Supports both the classic msigdbr schema (gs_cat/gs_subcat/gene_symbol) and
#' the new schema (gs_collection/gs_subcollection/db_gene_symbol).
#'
#' @param species "hs" or "mm"
#' @param database one of "HALLMARK", "GOBP", "KEGG"
#' @return tibble with gs_name and gene_symbol columns
rnaseq_msigdb_geneset <- function(species, database) {
  sp <- if_else(species == "hs", "human", "mouse")
  if (is.null(rnaseq_msigdb_cache[[sp]])) {
    rnaseq_msigdb_cache[[sp]] <- msigdbr::msigdbr(species = sp)
  }
  msigdb <- rnaseq_msigdb_cache[[sp]]

  if ("gs_cat" %in% colnames(msigdb)) {
    # classic schema (gs_cat/gs_subcat/gene_symbol)
    genesets <- if (database == "HALLMARK") {
      dplyr::filter(msigdb, gs_cat == "H")
    } else if (database == "KEGG") {
      dplyr::filter(msigdb, gs_subcat == "CP:KEGG")
    } else {
      dplyr::filter(msigdb, gs_subcat == "GO:BP")
    }
    symbol_col <- "gene_symbol"
  } else {
    # new schema (gs_collection/gs_subcollection/db_gene_symbol)
    genesets <- if (database == "HALLMARK") {
      dplyr::filter(msigdb, gs_collection == "H")
    } else if (database == "KEGG") {
      dplyr::filter(msigdb, gs_subcollection %in% c("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS"))
    } else {
      dplyr::filter(msigdb, gs_subcollection == "GO:BP")
    }
    symbol_col <- "db_gene_symbol"
  }

  dplyr::select(genesets, gs_name, gene_symbol = dplyr::all_of(symbol_col)) %>%
    dplyr::distinct()
}


# empty ORA result skeleton (keeps downstream bind_rows/plots stable)
rnaseq_ora_empty <- function() {
  tibble(ID = character(), Description = character(), GeneRatio = character(),
         BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
         qvalue = numeric(), geneID = character(), Count = integer(),
         direction = character())
}


#' Run enricher (ORA) for one gene list, returning a direction-tagged result tibble
#'
#' @param universe background genes actually tested (all genes in the Wald
#'   results). Without it, enricher would fall back to every gene in the
#'   gene-set table, which biases the enrichment p-values.
rnaseq_ora_direction <- function(genes, universe, genesets, pval_cutoff, min_setSize, max_setSize, direction) {
  if (length(genes) == 0) return(rnaseq_ora_empty())
  res <- suppressWarnings(tryCatch({
    x <- enricher(genes, universe = universe, pvalueCutoff = pval_cutoff, qvalueCutoff = 0.7,
                  minGSSize = min_setSize, maxGSSize = max_setSize,
                  TERM2GENE = genesets)
    x@result
  }, error = function(e) NULL))
  if (is.null(res) || nrow(res) == 0) return(rnaseq_ora_empty())
  res$direction <- direction
  as_tibble(res)
}


#' Run ORA (enricher) for all comparisons, up- and down-regulated genes separately
#'
#' @param results_sig named list of significant wald result tibbles
#' @param universe character vector of all tested genes (ORA background)
#' @param genesets TERM2GENE tibble from rnaseq_msigdb_geneset
#' @return named list of result tibbles with a direction column
rnaseq_run_ora <- function(results_sig, universe, genesets, pval_cutoff, min_setSize, max_setSize) {
  out <- furrr::future_map(results_sig, function(res) {
    up <- dplyr::filter(res, log2FoldChange > 0) %>% pull(gene_name)
    down <- dplyr::filter(res, log2FoldChange < 0) %>% pull(gene_name)
    bind_rows(
      rnaseq_ora_direction(up, universe, genesets, pval_cutoff, min_setSize, max_setSize, "Up"),
      rnaseq_ora_direction(down, universe, genesets, pval_cutoff, min_setSize, max_setSize, "Down")
    )
  }, .options = furrr::furrr_options(seed = TRUE))
  names(out) <- names(results_sig)
  out
}


#' Run GSEA for all comparisons in parallel (multisession, half of available cores)
#'
#' @param ranked_list named list of ranked gene vectors from rnaseq_gsea_ranked
#' @param genesets TERM2GENE tibble from rnaseq_msigdb_geneset
#' @param pval_cutoff numeric, pvalue cutoff for keeping significant terms
#' @param min_setSize numeric, minimum gene set size
#' @param max_setSize numeric, maximum gene set size
#' @return named list of clusterProfiler gseaResult objects (NULL entries dropped)
rnaseq_run_gsea <- function(ranked_list, genesets, pval_cutoff = 0.05,
                            min_setSize = 10, max_setSize = 500) {
  out <- furrr::future_map(ranked_list, function(genelist) {
    genelist <- genelist[is.finite(genelist)]
    if (length(genelist) < 10) return(NULL)
    suppressWarnings(tryCatch(
      clusterProfiler::GSEA(genelist, minGSSize = min_setSize, maxGSSize = max_setSize,
                            pvalueCutoff = pval_cutoff, TERM2GENE = genesets,
                            seed = TRUE, verbose = FALSE),
      error = function(e) NULL
    ))
  }, .options = furrr::furrr_options(seed = TRUE))
  names(out) <- names(ranked_list)
  out[!purrr::map_lgl(out, is.null)]
}


#' Volcano plot for a wald result tibble
rnaseq_plot_volcano <- function(res,
                                plot_title = "Volcano of DEGs",
                                size = 1.5, alpha = 0.6,
                                baseMean_thred = 10,
                                fc_thred = 2,
                                padj_thred = 0.05,
                                xlims = NULL, ylims = NULL,
                                legend = FALSE,
                                ns_resampling = 500,
                                color = c("#377eb8", "gray50", "#e41a1c")) {
  df <- filter(res, baseMean > baseMean_thred, !is.na(padj)) %>%
    mutate(`-log10(padj)` = -log10(padj))
  # take sig genes
  df.sig <- filter(df, padj < padj_thred, abs(log2FoldChange) > log2(fc_thred)) %>%
    mutate(direction = if_else(log2FoldChange > 0, "up", "down"))
  # sample the non-significant background cloud: the complement of df.sig, so
  # padj-significant but sub-threshold-LFC genes still render (as "ns")
  # instead of disappearing from the plot
  df.unsig <- filter(df, !(padj < padj_thred & abs(log2FoldChange) > log2(fc_thred))) %>%
    slice_sample(n = ns_resampling) %>%
    mutate(direction = "ns")
  # merge sig and un-sig
  df.merge <- bind_rows(df.sig, df.unsig)
  if (!is.null(xlims))
    df.merge <- mutate(df.merge, log2FoldChange = if_else(log2FoldChange > xlims[2], xlims[2],
                                                          if_else(log2FoldChange < xlims[1], xlims[1], log2FoldChange)))
  if (!is.null(ylims))
    df.merge <- mutate(df.merge, `-log10(padj)` = if_else(`-log10(padj)` > ylims[2], ylims[2],
                                                          if_else(`-log10(padj)` < ylims[1], ylims[1], `-log10(padj)`)))
  # plot
  plot <- ggplot(df.merge, aes(log2FoldChange, `-log10(padj)`, fill = direction)) +
    geom_point(size = size, alpha = alpha, shape = 21, color = "transparent", show.legend = legend) +
    scale_fill_manual(values = c(down = color[1], ns = color[2], up = color[3])) +
    geom_vline(xintercept = c(-log2(fc_thred), log2(fc_thred)), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(padj_thred), linetype = "dashed", color = "grey50") +
    labs(x = "log2(FoldChange)", y = "-log10(padj)",
         title = plot_title,
         subtitle = str_glue("padj_threshold: {padj_thred}, FoldChange_threshold: {fc_thred}, baseMean_threshold: {baseMean_thred}")) +
    theme_bw(14)
  if (!is.null(xlims)) plot <- plot + xlim(xlims)
  if (!is.null(ylims)) plot <- plot + ylim(ylims)

  return(plot)
}


#' ORA dot plot (top n terms per direction, faceted, lollipop style)
rnaseq_plot_ora_dot <- function(ora_df, n_top = 10) {
  df <- ora_df %>%
    dplyr::filter(!is.na(p.adjust)) %>%
    group_by(direction) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice(seq_len(min(n_top, n()))) %>%
    ungroup() %>%
    mutate(`-log10(p.adjust)` = -log10(p.adjust)) %>%
    dplyr::arrange(direction, p.adjust) %>%
    mutate(Description = factor(Description, levels = unique(Description)))

  ggplot(df, aes(`-log10(p.adjust)`, Description)) +
    geom_segment(aes(yend = Description), xend = 0, linewidth = 0.7, color = "grey50", alpha = 0.7) +
    geom_point(aes(size = Count, fill = `-log10(p.adjust)`), shape = 21) +
    facet_wrap(vars(direction), scales = "free_y", ncol = 1) +
    scale_size_continuous(range = c(2, 6)) +
    scale_fill_gradient(low = "#fee8c8", high = "#ef6548") +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 60)) +
    theme_bw() +
    labs(x = "-log10(p.adjust)", y = "", size = "DEG count")
}


#' GSEA dot plot (top n terms per direction, merged and ordered by NES)
rnaseq_plot_gsea_dot <- function(gsea_df, n_top = 10) {
  df <- gsea_df %>%
    dplyr::filter(!is.na(p.adjust)) %>%
    mutate(direction = if_else(NES >= 0, "Up-regulated", "Down-regulated"),
           Core_Count = map_int(core_enrichment, ~ length(str_split(.x, fixed("/"), simplify = TRUE))),
           `-log10(p.adjust)` = -log10(p.adjust)) %>%
    group_by(direction) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice(seq_len(min(n_top, n()))) %>%
    ungroup() %>%
    dplyr::arrange(NES) %>%
    mutate(Description = factor(Description, levels = unique(Description)))

  ggplot(df, aes(NES, Description)) +
    geom_segment(aes(yend = Description), xend = 0, linewidth = 0.7, color = "grey50", alpha = 0.7) +
    geom_point(aes(size = Core_Count, fill = NES), shape = 21) +
    scale_size_continuous(range = c(2, 6)) +
    scale_fill_gradient2(low = "#377eb8", mid = "white", high = "#e41a1c", midpoint = 0) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    theme_classic() +
    labs(y = "", size = "Core genes")
}


#' ggplot2 heatmap (replaces pheatmap)
#'
#' @param mat numeric matrix, rows = genes, columns = samples
#' @param scale "row" or "column", scale data before plotting
#' @param show_rownames logical, show y-axis gene labels
#' @param cluster_rows logical, order rows by hclust on scaled values
#' @param row_label_size numeric, size of y-axis labels
#' @param sample_colors named character vector (sample -> color) or NULL
#' @return ggplot object
rnaseq_plot_heatmap <- function(mat, scale = "row", show_rownames = TRUE,
                                cluster_rows = FALSE, row_label_size = 7,
                                sample_colors = NULL) {
  if (scale == "row") {
    mat <- t(scale(t(mat)))
  } else if (scale == "column") {
    mat <- scale(mat)
  }

  if (cluster_rows) {
    mat <- mat[hclust(dist(mat))$order, , drop = FALSE]
  }

  df <- as.data.frame(mat) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "value") %>%
    dplyr::mutate(gene = factor(gene, levels = unique(gene)))

  p <- ggplot(df, aes(sample, gene, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, na.value = "grey90") +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = row_label_size),
          panel.grid = element_blank(),
          legend.position = "right")

  if (!show_rownames) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }

  if (!is.null(sample_colors)) {
    # element_text() cannot color labels individually; replace the axis-b grob
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    gtab <- ggplot2::ggplotGrob(p)
    lbl_cols <- unname(sample_colors[unique(df$sample)])
    samples <- unique(df$sample)
    n_lab <- length(samples)
    x_npc <- seq(0.5 / n_lab, 1 - 0.5 / n_lab, length.out = n_lab)
    axis_idx <- which(gtab$layout$name == "axis-b")
    axis_lbl <- grid::textGrob(label = samples, x = x_npc, y = grid::unit(0.25, "npc"),
                               rot = 45, just = "right",
                               gp = grid::gpar(col = lbl_cols, fontsize = 10))
    gtab$grobs[[axis_idx]] <- grid::gTree(children = grid::gList(axis_lbl))
    p <- gtab
  }

  p
}


#' Standard DT table with buttons/scroller (downloads all rows)
rnaseq_dt <- function(data) {
  DT::renderDataTable(data, server = FALSE,
                      filter = "top", selection = "single", rownames = FALSE,
                      style = "bootstrap", class = "cell-border stripe",
                      extensions = c("Buttons", "Scroller"),
                      options = list(dom = "Brtip",
                                     columnDefs = list(list(className = "dt-center", targets = "_all")),
                                     buttons = c(I("colvis"), "copy", "excel"),
                                     deferRender = TRUE, scrollY = "150px", scroller = TRUE,
                                     scrollX = TRUE, searching = TRUE))
}


# =============================================================================
# Sub-module 1: Data Upload (counts + metadata + colData)
# =============================================================================

rnaseq_upload_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(status = "primary", solidHeader = TRUE, width = 12,
        fileInput(ns("counts_file"),
                  tagList("Raw counts (.csv / .xlsx)",
                          hint_tip("Gene x sample count matrix. Must include a 'gene_name' column; a 'gene_id' column is used as a fallback name for genes whose symbol is empty.")),
                  accept = rnaseq_file_accept),
        uiOutput(ns("counts_sheet_ui")),
        fileInput(ns("metadata_file"),
                  tagList("Sample metadata (.csv / .xlsx)",
                          hint_tip("Required columns: sample_name, batch, include. Any other column can be ticked as a design factor (e.g. treatment); the condition tested is the chosen factors pasted together. Rows with include = yes are analysed.")),
                  accept = rnaseq_file_accept),
        uiOutput(ns("metadata_sheet_ui")),
        shinyjs::disabled(actionButton(ns("upload"), "Upload",
                                       title = "Load both files; the metadata table becomes editable below.")),
        helpText("Upload both files, edit the metadata if needed, then continue with step 2.")
    ),
    uiOutput(ns("upload_tables_ui"))
  )
}


rnaseq_upload_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      # ----- upload counts file -----
      counts_file <- reactive({
        req(input$counts_file)
        input$counts_file
      })
      counts_filetype <- reactive(if_else(tools::file_ext(counts_file()$name) == "xlsx", "xlsx", "csv"))

      output$counts_sheet_ui <- renderUI({
        req(counts_file())
        if (counts_filetype() == "xlsx") {
          counts_sheets <- readxl::excel_sheets(counts_file()$datapath)
          radioButtons(ns("counts_data_sheet_selected"), label = NULL,
                       choiceNames = counts_sheets, choiceValues = counts_sheets,
                       selected = dplyr::last(counts_sheets), inline = TRUE)
        } else NULL
      })

      counts <- eventReactive(input$upload, {
        if (is.null(counts_file()$datapath)) return(NULL)
        showNotification("Read counts data ...", duration = 2, type = "message")
        if (counts_filetype() == "xlsx") {
          readxl::read_xlsx(counts_file()$datapath, sheet = input$counts_data_sheet_selected)
        } else readr::read_csv(counts_file()$datapath)
      })

      # ----- upload metadata file -----
      metadata_file <- reactive({
        req(input$metadata_file)
        input$metadata_file
      })
      metadata_filetype <- reactive(if_else(tools::file_ext(metadata_file()$name) == "xlsx", "xlsx", "csv"))

      output$metadata_sheet_ui <- renderUI({
        req(metadata_file())
        if (metadata_filetype() == "xlsx") {
          metadata_sheets <- readxl::excel_sheets(metadata_file()$datapath)
          radioButtons(ns("metadata_sheet_selected"), label = NULL,
                       choiceNames = metadata_sheets, choiceValues = metadata_sheets,
                       selected = dplyr::first(metadata_sheets), inline = TRUE)
        } else NULL
      })

      metadata <- eventReactive(input$upload, {
        if (is.null(metadata_file()$datapath)) return(NULL)
        showNotification("Read metadata ...", duration = 1, type = "message")
        if (metadata_filetype() == "xlsx") {
          readxl::read_xlsx(metadata_file()$datapath, sheet = input$metadata_sheet_selected,
                            col_types = "text")
        } else readr::read_csv(metadata_file()$datapath, col_types = "text")
      })

      # enable upload button when both files are selected
      observe({
        req(counts_file(), metadata_file())
        shinyjs::enable("upload")
      })

      # ----- edit metadata -----
      output$metadata_table <- renderRHandsontable({
        req(metadata())
        tbl <- rhandsontable(metadata(), useTypes = TRUE) %>%
          hot_table(highlightCol = TRUE, highlightRow = TRUE, stretchH = "all") %>%
          hot_cols(columnSorting = TRUE)
        for (nm in intersect(c("sample_name", "batch"), colnames(metadata()))) {
          tbl <- hot_col(tbl, nm, readOnly = TRUE)
        }
        tbl
      })

      metadata_final <- reactive({
        if (!is.null(input$metadata_table)) {
          hot_to_r(input$metadata_table)
        } else metadata()
      })

      # ----- build colData -----
      output$choose_factors_ui <- renderUI({
        req(metadata())
        cn <- colnames(metadata())
        required <- c("sample_name", "batch", "include")
        if (!all(required %in% cn)) {
          showNotification(paste0(
            "Metadata must contain the columns: ", paste0(required, collapse = ", "),
            ". Found columns: ", paste(cn, collapse = ", ")),
            type = "error", duration = 15)
          return(helpText("Metadata is missing required column(s): ",
                          paste(setdiff(required, cn), collapse = ", ")))
        }
        factor_cols <- setdiff(cn, rnaseq_metadata_reserved_cols)
        if (length(factor_cols) == 0) {
          showNotification(
            "Metadata has no column left for a design factor - add a condition column (e.g. treatment) besides sample_name, batch and include.",
            type = "error", duration = 15)
          return(helpText("No design factor available: add a condition column (e.g. treatment) besides sample_name, batch and include."))
        }
        checkboxGroupInput(ns("choose_factors"), label = NULL, inline = TRUE,
                           choiceNames = factor_cols,
                           choiceValues = factor_cols,
                           selected = factor_cols)
      })

      coldata <- reactive({
        req(input$choose_factors)
        data <- metadata_final()
        if (anyDuplicated(data$sample_name) > 0) {
          showNotification("Metadata contains duplicated sample_name values; deduplicate before running DESeq2.",
                           type = "warning", duration = 10)
        }
        data$condition <- apply(data[, input$choose_factors, drop = FALSE], 1, paste0, collapse = ".")
        data <- dplyr::filter(data, tolower(trimws(include)) == "yes") %>%
          .[, c("sample_name", "batch", "condition"), drop = FALSE] %>%
          tibble::column_to_rownames("sample_name")
        data
      })

      output$coldata_table <- renderRHandsontable({
        req(coldata())
        rhandsontable(coldata(), readOnly = TRUE, stretchH = "all", rowHeaderWidth = 200) %>%
          hot_cols(columnSorting = TRUE)
      })

      # ----- table boxes (hidden until data is uploaded) -----
      output$upload_tables_ui <- renderUI({
        req(metadata())
        fluidRow(
          box(status = "primary", solidHeader = TRUE, width = 6,
              rHandsontableOutput(ns("metadata_table"))
          ),
          box(status = "primary", solidHeader = TRUE, width = 6,
              fluidRow(
                column(9, div(align = "center", uiOutput(ns("choose_factors_ui")))),
                column(3)
              ),
              fluidRow(
                column(12, rHandsontableOutput(ns("coldata_table")))
              )
          )
        )
      })

      return(list(counts = counts, coldata = coldata))
    }
  )
}


# =============================================================================
# Sub-module 2: QC & DESeq2 (PCA, correlation, variable genes, dds)
# =============================================================================

rnaseq_qc_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(width = 12, status = "primary", solidHeader = TRUE,
          title = tagList("QC & DESeq2",
                          hint_tip("Run DESeq2 once here; the PCA, correlation and heatmap use variance-stabilised values, and the later tabs reuse the fitted model.")),
          fluidRow(
            column(2,
                   radioButtons(ns("filter_symbol"), tagList("Empty symbols",
                                hint_tip("Remove genes whose gene_name is empty, or backfill them from gene_id.")),
                                selected = "true", inline = TRUE, choices = c(true = "true", false = "false"))
            ),
            column(2,
                   numericInput(ns("min_counts"), tagList("Counts threshold",
                                 hint_tip("Genes with mean count below this value are removed before the test.")),
                                 value = 10, min = 1, step = 1)
            ),
            column(2,
                   radioButtons(ns("pca_show_label"), "Labels on PCA",
                                selected = "true", inline = TRUE, choices = c(true = "true", false = "false"))
            ),
            column(2,
                   radioButtons(ns("calculate_source"), tagList("QC values",
                                hint_tip("Values behind the QC plots: vst = variance-stabilising transform (recommended); counts_norm = DESeq2 size-factor-normalised counts.")),
                                selected = "vst", inline = TRUE, choices = c(vst = "vst", counts_norm = "counts_norm"))
            ),
            column(4, align = "center", br(),
                   actionButton(ns("apply"), "Run DESeq2", class = "btn-default",
                                title = "Fit the DESeq2 model (batch is included as a covariate when it has more than one level).")
            )
          )
      )
    ),
    uiOutput(ns("qc_results_ui"))
  )
}


rnaseq_qc_server <- function(id, counts, coldata) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      # vst() errors on small gene counts; fall back to the classical VST
      vst_transform <- function(dds_obj) {
        tryCatch(
          DESeq2::vst(dds_obj, blind = FALSE),
          error = function(e) DESeq2::varianceStabilizingTransformation(dds_obj, blind = FALSE)
        )
      }

      # ----- run DESeq2 -----
      dds <- eventReactive(input$apply, {
        req(counts(), coldata())
        if (!all(rownames(coldata()) %in% colnames(counts()))) {
          showNotification("Not all selected metadata sample_name are contained in counts data, check it!",
                           duration = 5, type = "error")
          return(NULL)
        }
        if (!"gene_name" %in% colnames(counts())) {
          showNotification("Counts data must contain a 'gene_name' column!",
                           duration = 5, type = "error")
          return(NULL)
        }

        withProgress(message = "Running DESeq2 ...", value = 0.1, {
          rnaseq_build_dds(counts(), coldata(),
                           filter_symbol = as.logical(input$filter_symbol),
                           min_counts = input$min_counts)
        })
      })

      mat <- reactive({
        req(dds())
        if (input$calculate_source == "vst") {
          assay(vst_transform(dds()))
        } else {
          DESeq2::counts(dds(), normalized = TRUE)
        }
      })

      # ----- PCA -----
      pca_p <- reactive({
        req(dds())
        pcaData <- plotPCA(vst_transform(dds()), intgroup = c("condition"), returnData = TRUE)
        percentVar <- round(100 * attr(pcaData, "percentVar"))
        p <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
          geom_point(size = 2) +
          xlab(paste0("PC1: ", percentVar[1], "% variance")) +
          ylab(paste0("PC2: ", percentVar[2], "% variance")) +
          # dashed ellipses share the point colour scale
          stat_ellipse(aes(color = condition), level = 0.9, linetype = 2) +
          coord_fixed() +
          theme_bw(16) +
          theme(legend.position = "top")
        if (as.logical(input$pca_show_label)) p <- p + geom_text_repel(aes(label = name), show.legend = FALSE)
        p
      })

      output$pca_plot <- renderPlot({
        req(pca_p())
        pca_p()
      }, res = 96)

      output$dl_pca <- downloadHandler(
        filename = function() paste0(format(Sys.Date(), "%Y-%m-%d"), "_pca.pdf"),
        content = function(file) {
          req(pca_p())
          export_plot_pdf(pca_p(), file, width = 7, height = 6)
        })

      # ----- sample correlation -----
      cor_plot_fn <- reactive({
        req(mat())
        function() corrplot.mixed(cor(mat()),
                                  lower = "number", upper = "ellipse",
                                  order = "AOE",
                                  tl.pos = "lt", tl.col = "black",
                                  number.cex = 0.7)
      })

      output$cor_plot <- renderPlot(cor_plot_fn()())

      output$dl_cor <- downloadHandler(
        filename = function() paste0(format(Sys.Date(), "%Y-%m-%d"), "_correlation.pdf"),
        content = function(file) {
          req(mat())
          export_plot_pdf(cor_plot_fn(), file, width = 7, height = 7)
        })

      # ----- top variable genes heatmap -----
      top_var_p <- reactive({
        req(mat())
        n_top <- min(2000, nrow(mat()))
        top_variable_genes <- order(matrixStats::rowVars(mat()), decreasing = TRUE)[seq_len(n_top)]
        cd <- as.data.frame(colData(dds()))
        sample_condition <- cd[colnames(mat()), "condition"]
        cond_colors <- setNames(scales::hue_pal()(length(unique(sample_condition))), unique(sample_condition))
        rnaseq_plot_heatmap(mat()[top_variable_genes, ],
                            cluster_rows = TRUE,
                            show_rownames = FALSE,
                            sample_colors = setNames(cond_colors[sample_condition], colnames(mat())))
      })

      output$top_var_plot <- renderPlot({
        req(top_var_p())
        p <- top_var_p()
        if (inherits(p, "gtable")) grid::grid.draw(p) else p
      })

      output$dl_topvar <- downloadHandler(
        filename = function() paste0(format(Sys.Date(), "%Y-%m-%d"), "_variable_genes_heatmap.pdf"),
        content = function(file) {
          req(top_var_p())
          export_plot_pdf(top_var_p(), file, width = 8, height = 10)
        })

      # ----- QC result boxes (hidden until DESeq2 has run) -----
      output$qc_results_ui <- renderUI({
        req(dds())
        fluidRow(
          box(width = 6, status = "primary", solidHeader = TRUE,
              div(style = "position:relative",
                  tags$div(style = "position:absolute;top:5px;right:5px;z-index:10;",
                           downloadButton(ns("dl_pca"), label = NULL, icon = icon("download"), class = "btn-xs", title = "Download PCA plot (PDF)")),
                  plotOutput(ns("pca_plot"), height = "400px")),
              br(),
              div(style = "position:relative",
                  tags$div(style = "position:absolute;top:5px;right:5px;z-index:10;",
                           downloadButton(ns("dl_cor"), label = NULL, icon = icon("download"), class = "btn-xs", title = "Download sample-correlation plot (PDF)")),
                  plotOutput(ns("cor_plot")))
          ),
          box(width = 6, status = "primary", solidHeader = TRUE,
              div(style = "position:relative",
                  tags$div(style = "position:absolute;top:5px;right:5px;z-index:10;",
                           downloadButton(ns("dl_topvar"), label = NULL, icon = icon("download"), class = "btn-xs", title = "Download top-variable-genes heatmap (PDF)")),
                  plotOutput(ns("top_var_plot"), height = "815px"))
          )
        )
      })

      return(list(dds = dds, mat = mat))
    }
  )
}


# =============================================================================
# Sub-module 3: Differential Expression (Wald comparisons, tables, volcano)
# =============================================================================

rnaseq_dge_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(width = 12, status = "primary", solidHeader = TRUE,
          title = tagList("Differential Expression",
                          hint_tip("Add one or more condition comparisons, run the Wald tests, then browse significant genes and the volcano plot. The thresholds below also drive the volcano cut lines.")),
          fluidRow(
            column(4, align = "left",
                   numericInput(ns("baseMean_threshold"), tagList("Min baseMean",
                                hint_tip("Keep genes with DESeq2 baseMean above this value.")),
                                100, 1, Inf, width = "100%")
            ),
            column(4, align = "left",
                   numericInput(ns("fc_threshold"), tagList("Min fold change",
                                hint_tip("Absolute fold-change cut-off for calling a gene up/down.")),
                                1.5, 1, Inf, width = "100%")
            ),
            column(4, align = "left",
                   numericInput(ns("padj_threshold"), tagList("Max adjusted p",
                                hint_tip("Benjamini-Hochberg adjusted p-value cut-off.")),
                                0.05, 0, Inf, width = "100%")
            )
          ),
          fluidRow(
            column(12,
                   div(id = ns("compare_container")),
                   fluidRow(
                     column(2, align = "center",
                            actionButton(ns("add_compare"), "Add Comparison",
                                         title = "Add a contrast-vs-control row; one per comparison.")),
                     column(2, align = "center",
                            actionButton(ns("apply"), "Run Comparisons", class = "btn-default",
                                         title = "Wald-test every comparison; fold changes are shrunken with ashr.")),
                     column(4, align = "center",
                            downloadButton(ns("dl_dge_xlsx"), "Download Results (xlsx)", icon = icon("table"),
                                           title = "Full Wald results (all genes) for every comparison, one sheet each."))
                   )
            )
          )
      )
    ),
    fluidRow(
      uiOutput(ns("sig_tables_ui"))
    ),
    uiOutput(ns("volcano_box_ui"))
  )
}


rnaseq_dge_server <- function(id, dds) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      # ----- dynamic comparison UI -----
      counter <- reactiveVal(integer(0))
      next_id <- reactiveVal(1L)

      observeEvent(input$add_compare, {
        req(dds())
        conditions <- unique(colData(dds())$condition)

        # ids of existing comparison UIs; next_id never gets recycled so a
        # stale remove-button value can never delete a freshly added row
        n_now <- next_id()
        next_id(n_now + 1L)
        counter(c(counter(), n_now))

        insertUI(selector = paste0("#", ns("compare_container")), where = "beforeEnd", immediate = TRUE,
                 ui = fluidRow(
                   div(id = ns(paste0("compare_insertUI_", n_now)),
                       column(5,
                              selectInput(ns(paste0("compare_contrast_", n_now)), "Contrast (of interest)",
                                          conditions, selected = conditions[2])
                       ),
                       column(5,
                              selectInput(ns(paste0("compare_control_", n_now)), "Control (reference)",
                                          conditions, selected = conditions[1])
                       ),
                       column(2,
                              actionButton(ns(paste0("remove_compare_", n_now)), "Remove",
                                           title = "Remove this comparison row."))
                   )
                 )
        )
      })

      observe({
        req(length(counter()) > 0)
        lapply(counter(), function(i) {
          remove_btn <- input[[paste0("remove_compare_", i)]]
          if (isTruthy(remove_btn)) {
            removeUI(selector = paste0("#", ns(paste0("compare_insertUI_", i))), immediate = TRUE)
            counter(counter()[-which(counter() == i)])
          }
        })
      })

      # ----- Wald test results -----
      results <- eventReactive(input$apply, {
        req(dds())
        if (length(counter()) == 0) {
          showNotification("Add at least one comparison first", type = "warning", duration = 6)
          return(NULL)
        }

        compare_list <- purrr::map(counter(), function(i) {
          c(input[[paste0("compare_contrast_", i)]], input[[paste0("compare_control_", i)]])
        })
        bad <- purrr::map_lgl(compare_list, ~ is.na(.x[1]) || is.na(.x[2]) || .x[1] == .x[2])
        if (any(bad)) {
          showNotification("Contrast and control must be different conditions; ignoring those comparison(s).",
                           type = "warning", duration = 8)
          compare_list <- compare_list[!bad]
        }
        if (length(compare_list) == 0) {
          showNotification("No valid comparisons selected", type = "warning", duration = 6)
          return(NULL)
        }

        withProgress(message = "Wald test ...", value = 0.1, {
          out <- purrr::imap(compare_list, function(cmp, i) {
            incProgress(1 / length(compare_list), detail = paste0(cmp[1], " vs ", cmp[2]))
            rnaseq_wald_result(dds(), cmp[1], cmp[2])
          })
        })
        names(out) <- purrr::map_chr(compare_list, paste0, collapse = "_vs_")
        out
      })

      # ----- filter significant -----
      results_sig <- reactive({
        req(results())
        purrr::map(results(), ~ dplyr::filter(.x,
                                              baseMean > input$baseMean_threshold,
                                              abs(log2FoldChange) > log2(input$fc_threshold),
                                              padj < input$padj_threshold))
      })

      # ----- full result workbook: one sheet per comparison -----
      output$dl_dge_xlsx <- downloadHandler(
        filename = function() paste0('RNAseq_', Sys.Date(), '_DGE_results.xlsx'),
        content = function(file) {
          res_all <- results()
          req(length(res_all) > 0)
          write_xlsx_sheets(res_all, file)
        }
      )

      # ----- result tables -----
      output$sig_tables_ui <- renderUI({
        req(results_sig())
        tabpanel_list <- lapply(seq_along(results_sig()), function(i) {
          tabPanel(title = names(results_sig())[i],
                   rnaseq_dt(results_sig()[[i]]))
        })
        do.call(shinydashboard::tabBox, args = c(id = ns("sig_table_tabbox"), width = 12, tabpanel_list))
      })

      # ----- volcano -----
      volcano_p <- reactive({
        req(results(), input$sig_table_tabbox)
        # feed the FULL result table so non-significant genes render as the
        # gray background cloud; rnaseq_plot_volcano applies the thresholds
        df <- results()[[input$sig_table_tabbox]]
        rnaseq_plot_volcano(df,
                            plot_title = input$sig_table_tabbox,
                            baseMean_thred = input$baseMean_threshold,
                            fc_thred = input$fc_threshold,
                            padj_thred = input$padj_threshold)
      })

      output$sig_volcano <- renderPlot({
        req(volcano_p())
        volcano_p()
      }, res = 96)

      output$dl_volcano <- downloadHandler(
        filename = function() {
          cmp <- gsub("[^[:alnum:].-]", "_", input$sig_table_tabbox)
          paste0(format(Sys.Date(), "%Y-%m-%d"), "_", cmp, "_volcano.pdf")
        },
        content = function(file) {
          req(volcano_p())
          export_plot_pdf(volcano_p(), file, width = 7, height = 6)
        })

      # ----- volcano box (hidden until comparison results exist) -----
      output$volcano_box_ui <- renderUI({
        req(results())
        fluidRow(
          box(width = 6, status = "primary", solidHeader = TRUE,
              div(style = "position:relative",
                  tags$div(style = "position:absolute;top:5px;right:5px;z-index:10;",
                           downloadButton(ns("dl_volcano"), label = NULL, icon = icon("download"), class = "btn-xs", title = "Download volcano plot (PDF)")),
                  plotOutput(ns("sig_volcano"), height = "400px")
              )
          )
        )
      })

      return(list(results = results, results_sig = results_sig))
    }
  )
}


# =============================================================================
# Sub-module 4: Enrichment (ORA + GSEA)
# =============================================================================

rnaseq_ora_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(status = "primary", solidHeader = TRUE, width = 12,
        title = tagList("Enrichment",
                        hint_tip("ORA tests the significant genes (from the thresholds set in step 3) for over-represented gene sets; GSEA ranks all genes by the Wald statistic. Both run against MSigDB.")),
        fluidRow(
          column(2, align = "center",
                 radioButtons(ns("species"), "Species", selected = "hs", inline = TRUE,
                              choiceNames = c("hs", "mm"), choiceValues = c("hs", "mm"))
          ),
          column(2, align = "center",
                 radioButtons(ns("database"), tagList("Database",
                              hint_tip("HALLMARK: 50 hallmark sets; GOBP: Gene Ontology biological process; KEGG: KEGG pathways.")),
                              selected = "HALLMARK", inline = TRUE,
                              choiceNames = c("HALLMARK", "GOBP", "KEGG"),
                              choiceValues = c("HALLMARK", "GOBP", "KEGG"))
          ),
          column(2,
                 numericInput(ns("min_setSize"), tagList("Min set size",
                              hint_tip("Gene sets smaller than this are skipped.")),
                              10, min = 1, max = 100)
          ),
          column(2,
                 numericInput(ns("max_setSize"), tagList("Max set size",
                              hint_tip("Gene sets larger than this are skipped.")),
                              500, min = 50, max = 1000)
          ),
          column(1,
                 numericInput(ns("ora_pval_cutoff"), "ORA p cutoff", 0.05, min = 0, max = 0.2)
          ),
          column(1, align = "center", br(),
                 actionButton(ns("ora_apply"), "Run ORA", class = "btn-default", width = "100%",
                              title = "Over-representation analysis on the significant genes of every comparison."))
          ,
          column(1,
                 numericInput(ns("gsea_pval_cutoff"), "GSEA p cutoff", 0.05, min = 0, max = 1)
          ),
          column(1, align = "center", br(),
                 actionButton(ns("gsea_apply"), "Run GSEA", class = "btn-default", width = "100%",
                              title = "Gene-set enrichment analysis on all genes ranked by Wald statistic."))
        )
    ),
    uiOutput(ns("enrichment_tabs"))
  )
}


rnaseq_ora_server <- function(id, results, results_sig, mat) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      # ----- shared MSigDB geneset (loaded at most once per species/database) -----
      genesets <- reactiveVal(NULL)
      genesets_key <- reactiveVal("")

      observeEvent(c(input$ora_apply, input$gsea_apply), ignoreInit = TRUE, {
        req(input$species, input$database)
        key <- paste(input$species, input$database, sep = "/")
        if (!identical(genesets_key(), key)) {
          withProgress(message = "Loading MSigDB genesets ...", value = 0.3, {
            genesets(rnaseq_msigdb_geneset(input$species, input$database))
          })
          genesets_key(key)
        }
      })

      # ----- ORA (runs immediately on click, cached per parameters/data) -----
      ora_list <- reactiveVal(NULL)
      ora_cache <- reactiveValues(cache = list())

      observeEvent(input$ora_apply, ignoreInit = TRUE, {
        req(genesets())
        res_sig <- tryCatch(results_sig(), error = function(e) NULL)
        if (is.null(res_sig) || length(res_sig) == 0) {
          showNotification("Run comparisons in the 'Differential Expression' tab first",
                           duration = 5, type = "error")
          return()
        }
        key <- paste(input$species, input$database, input$ora_pval_cutoff,
                     input$min_setSize, input$max_setSize, rlang::hash(res_sig), sep = "|")
        cached <- ora_cache$cache[[key]]
        if (!is.null(cached)) {
          ora_list(cached)
          showNotification("ORA result loaded from cache", duration = 1, type = "message")
          return()
        }
        # ORA background: every gene tested in the Wald step (i.e. retained in
        # the full results), not the enricher default of all genesets' genes
        res_all <- tryCatch(results(), error = function(e) NULL)
        universe <- if (!is.null(res_all)) {
          unique(unlist(purrr::map(res_all, ~ .x$gene_name), use.names = FALSE))
        } else character(0)
        if (length(universe) == 0) {
          showNotification("No tested genes found for the ORA background", duration = 5, type = "error")
          return()
        }
        withProgress(message = "ORA enrichment ...", value = 0.1, {
          ora_list(rnaseq_run_ora(res_sig, universe, genesets(),
                                  input$ora_pval_cutoff, input$min_setSize, input$max_setSize))
        })
        ora_cache$cache[[key]] <- ora_list()
      })

      output$ora_tables_ui <- renderUI({
        req(ora_list())
        tabpanel_list <- lapply(seq_along(ora_list()), function(i) {
          tabPanel(title = names(ora_list())[i], rnaseq_dt(ora_list()[[i]]))
        })
        do.call(shinydashboard::tabBox, args = c(id = ns("ora_tabbox"), width = 12, tabpanel_list))
      })

      ora_dot_p <- reactive({
        req(ora_list(), input$ora_tabbox)
        validate(need(nrow(ora_list()[[input$ora_tabbox]]) > 0, "No significantly enriched terms"))
        rnaseq_plot_ora_dot(ora_list()[[input$ora_tabbox]])
      })

      output$ora_dotplot <- renderPlot({
        req(ora_dot_p())
        ora_dot_p()
      }, res = 96)

      output$dl_ora_dot <- downloadHandler(
        filename = function() {
          cmp <- gsub("[^[:alnum:].-]", "_", input$ora_tabbox)
          paste0(format(Sys.Date(), "%Y-%m-%d"), "_", cmp, "_ORA_dotplot.pdf")
        },
        content = function(file) {
          req(ora_dot_p())
          export_plot_pdf(ora_dot_p(), file, width = 8, height = 6)
        })

      # ----- ORA result workbook (all comparisons) -----
      output$dl_ora_xlsx <- downloadHandler(
        filename = function() paste0('RNAseq_', Sys.Date(), '_ORA_results.xlsx'),
        content = function(file) {
          write_xlsx_sheets(ora_list(), file)
        })

      # ----- pathway genes heatmap (DEGs) -----
      output$ora_heatmap_ui <- renderUI({
        req(ora_list(), input$ora_tabbox)
        validate(need(nrow(ora_list()[[input$ora_tabbox]]) > 0, "No significantly enriched terms"))
        choices <- unique(ora_list()[[input$ora_tabbox]]$ID)
        selectizeInput(ns("ora_pathway"), "Pathway",
                       choices = choices, selected = choices[1],
                       options = list(placeholder = "Type to search...", maxOptions = 10000))
      })

      pathway_genes <- reactive({
        req(mat(), genesets(), input$ora_pathway)
        genesets() %>%
          dplyr::filter(gs_name == input$ora_pathway) %>%
          pull(gene_symbol) %>%
          intersect(rownames(mat()))
      })

      ora_heat_p <- reactive({
        req(pathway_genes(), results_sig(), input$ora_tabbox)
        degs <- intersect(pathway_genes(), results_sig()[[input$ora_tabbox]]$gene_name)
        validate(need(length(degs) > 0, "No DEGs found for this pathway"))
        rnaseq_plot_heatmap(mat()[degs, , drop = FALSE],
                            cluster_rows = FALSE,
                            show_rownames = length(degs) <= 50,
                            row_label_size = 7)
      })

      output$ora_heatmap_deg <- renderPlot(ora_heat_p())

      output$dl_ora_heat <- downloadHandler(
        filename = function() {
          cmp <- gsub("[^[:alnum:].-]", "_", input$ora_tabbox)
          pwy <- gsub("[^[:alnum:].-]", "_", input$ora_pathway)
          paste0(format(Sys.Date(), "%Y-%m-%d"), "_", cmp, "_", pwy, "_DEG_heatmap.pdf")
        },
        content = function(file) {
          req(ora_heat_p())
          ggsave(file, ora_heat_p(), width = 8, height = 6, device = pdf)
        })

      # ----- GSEA (runs immediately on click, cached per parameters/data) -----
      gsea_objs <- reactiveVal(NULL)
      gsea_cache <- reactiveValues(cache = list())

      observeEvent(input$gsea_apply, ignoreInit = TRUE, {
        req(genesets())
        res_all <- tryCatch(results(), error = function(e) NULL)
        if (is.null(res_all) || length(res_all) == 0) {
          showNotification("Run comparisons in the 'Differential Expression' tab first",
                           duration = 5, type = "error")
          return()
        }
        key <- paste(input$species, input$database, input$gsea_pval_cutoff,
                     input$min_setSize, input$max_setSize, rlang::hash(res_all), sep = "|")
        cached <- gsea_cache$cache[[key]]
        if (!is.null(cached)) {
          gsea_objs(cached)
          showNotification("GSEA result loaded from cache", duration = 1, type = "message")
          return()
        }
        withProgress(message = "Running GSEA ...", value = 0.1, {
          gsea_objs(rnaseq_run_gsea(rnaseq_gsea_ranked(res_all), genesets(),
                                    pval_cutoff = input$gsea_pval_cutoff,
                                    min_setSize = input$min_setSize, max_setSize = input$max_setSize))
        })
        gsea_cache$cache[[key]] <- gsea_objs()
      })

      gsea_list <- reactive({
        req(gsea_objs())
        purrr::map(gsea_objs(), ~ as_tibble(.x@result))
      })

      output$gsea_tables_ui <- renderUI({
        req(gsea_list())
        tabpanel_list <- lapply(seq_along(gsea_list()), function(i) {
          tabPanel(title = names(gsea_list())[i], rnaseq_dt(gsea_list()[[i]]))
        })
        do.call(shinydashboard::tabBox, args = c(id = ns("gsea_tabbox"), width = 12, tabpanel_list))
      })

      gsea_dot_p <- reactive({
        req(gsea_list(), input$gsea_tabbox)
        validate(need(nrow(gsea_list()[[input$gsea_tabbox]]) > 0, "No significant GSEA terms"))
        rnaseq_plot_gsea_dot(gsea_list()[[input$gsea_tabbox]])
      })

      output$gsea_dotplot <- renderPlot({
        req(gsea_dot_p())
        gsea_dot_p()
      }, res = 96)

      output$dl_gsea_dot <- downloadHandler(
        filename = function() {
          cmp <- gsub("[^[:alnum:].-]", "_", input$gsea_tabbox)
          paste0(format(Sys.Date(), "%Y-%m-%d"), "_", cmp, "_GSEA_dotplot.pdf")
        },
        content = function(file) {
          req(gsea_dot_p())
          export_plot_pdf(gsea_dot_p(), file, width = 8, height = 6)
        })

      # ----- GSEA result workbook (all comparisons) -----
      output$dl_gsea_xlsx <- downloadHandler(
        filename = function() paste0('RNAseq_', Sys.Date(), '_GSEA_results.xlsx'),
        content = function(file) {
          write_xlsx_sheets(gsea_list(), file)
        })

      # pathway selector for running score plot (follows the comparison result tab)
      output$gsea_pathway_ui <- renderUI({
        req(gsea_list(), input$gsea_tabbox)
        validate(need(length(gsea_list()) > 0, "No GSEA results (check that comparisons were run)"))
        validate(need(nrow(gsea_list()[[input$gsea_tabbox]]) > 0, "No significant GSEA terms"))
        selectizeInput(ns("gsea_pathway"), "Pathway",
                    choices = gsea_list()[[input$gsea_tabbox]]$Description,
                    selected = gsea_list()[[input$gsea_tabbox]]$Description[1],
                    options = list(placeholder = "Type to search...", maxOptions = 10000))
      })

      gsea_gsea_p <- reactive({
        req(gsea_objs(), input$gsea_tabbox, input$gsea_pathway)
        obj <- gsea_objs()[[input$gsea_tabbox]]
        validate(need(!is.null(obj) && nrow(obj@result) > 0, "No significant GSEA terms"))
        # gseaplot2 expects an ID or index, not the (wrapped) description string
        gene_set_id <- obj@result$ID[match(input$gsea_pathway, obj@result$Description)]
        if (is.na(gene_set_id)) gene_set_id <- input$gsea_pathway
        enrichplot::gseaplot2(obj, geneSetID = gene_set_id,
                              title = input$gsea_pathway)
      })

      output$gsea_gseaplot <- renderPlot(gsea_gsea_p())

      output$dl_gsea_plot <- downloadHandler(
        filename = function() {
          cmp <- gsub("[^[:alnum:].-]", "_", input$gsea_tabbox)
          pwy <- gsub("[^[:alnum:].-]", "_", input$gsea_pathway)
          paste0(format(Sys.Date(), "%Y-%m-%d"), "_", cmp, "_", pwy, "_GSEA_running_score.pdf")
        },
        content = function(file) {
          req(gsea_gsea_p())
          export_plot_pdf(gsea_gsea_p(), file, width = 8, height = 7)
        })

      # ----- enrichment result tabs (hidden until ORA/GSEA results exist) -----
      output$enrichment_tabs <- renderUI({
        tab_list <- list()
        if (!is.null(ora_list())) {
          tab_list$ORA <- tabPanel("ORA",
            fluidRow(uiOutput(ns("ora_tables_ui"))),
            fluidRow(
              box(width = 6, status = "primary", solidHeader = TRUE,
                  div(style = "position:relative",
                      tags$div(style = "position:absolute;top:5px;right:5px;z-index:10;",
                               downloadButton(ns("dl_ora_dot"), label = NULL, icon = icon("download"), class = "btn-xs", title = "Download ORA dot plot (PDF)")),
                      plotOutput(ns("ora_dotplot"), height = "490px"))),
              box(width = 6, status = "primary", solidHeader = TRUE,
                  uiOutput(ns("ora_heatmap_ui")),
                  h5("DEGs in the selected pathway"),
                  div(style = "position:relative",
                      tags$div(style = "position:absolute;top:5px;right:5px;z-index:10;",
                               downloadButton(ns("dl_ora_heat"), label = NULL, icon = icon("download"), class = "btn-xs", title = "Download pathway DEG heatmap (PDF)")),
                      plotOutput(ns("ora_heatmap_deg"), height = "400px")))
            ),
            fluidRow(
              column(4, downloadButton(ns("dl_ora_xlsx"), "Download ORA Results (xlsx)", icon = icon("table")))
            )
          )
        }
        if (!is.null(gsea_objs())) {
          tab_list$GSEA <- tabPanel("GSEA",
            fluidRow(uiOutput(ns("gsea_tables_ui"))),
            fluidRow(
              box(width = 6, status = "primary", solidHeader = TRUE,
                  div(style = "position:relative",
                      tags$div(style = "position:absolute;top:5px;right:5px;z-index:10;",
                               downloadButton(ns("dl_gsea_dot"), label = NULL, icon = icon("download"), class = "btn-xs", title = "Download GSEA dot plot (PDF)")),
                      plotOutput(ns("gsea_dotplot"), height = "490px"))),
              box(width = 6, status = "primary", solidHeader = TRUE,
                  uiOutput(ns("gsea_pathway_ui")),
                  div(style = "position:relative",
                      tags$div(style = "position:absolute;top:5px;right:5px;z-index:10;",
                               downloadButton(ns("dl_gsea_plot"), label = NULL, icon = icon("download"), class = "btn-xs", title = "Download GSEA running-score plot (PDF)")),
                      plotOutput(ns("gsea_gseaplot"), height = "430px")))
            ),
            fluidRow(
              column(4, downloadButton(ns("dl_gsea_xlsx"), "Download GSEA Results (xlsx)", icon = icon("table")))
            )
          )
        }
        if (length(tab_list) == 0) return(NULL)
        do.call(tabsetPanel, unname(tab_list))
      })
    }
  )
}


# =============================================================================
# Outer Module
# =============================================================================

rnaseqUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$head(tags$style(htmltools::HTML("
      .rnaseq-steptabs .nav-tabs-custom { background: #3c8dbc; box-shadow: none; }
      .rnaseq-steptabs .nav-tabs-custom > .nav-tabs { border-bottom-color: transparent; }
      .rnaseq-steptabs .nav-tabs-custom > .nav-tabs > li { border-top: none; margin-bottom: 0; }
      .rnaseq-steptabs .nav-tabs-custom > .nav-tabs > li > a,
      .rnaseq-steptabs .nav-tabs-custom > .nav-tabs > li > a:hover { color: #fff; background: transparent; }
      .rnaseq-steptabs .nav-tabs-custom > .nav-tabs > li:not(.active) > a:hover,
      .rnaseq-steptabs .nav-tabs-custom > .nav-tabs > li:not(.active) > a:focus {
        background: rgba(255,255,255,0.15); color: #fff;
      }
      .rnaseq-steptabs .nav-tabs-custom > .nav-tabs > li.active { border-top: none; }
      .rnaseq-steptabs .nav-tabs-custom > .nav-tabs > li.active > a,
      .rnaseq-steptabs .nav-tabs-custom > .nav-tabs > li.active:hover > a {
        background-color: #fff; color: #3c8dbc; border-top-color: transparent;
      }
    "))),
    fluidRow(
      useShinyjs(),
      div(class = "rnaseq-steptabs",
        do.call(shinydashboard::tabBox,
                args = c(width = 12,
                         list(
                           tabPanel("1. Upload Data", rnaseq_upload_ui(ns("upload"))),
                           tabPanel("2. QC & DESeq2", rnaseq_qc_ui(ns("qc"))),
                           tabPanel("3. Differential Expression", rnaseq_dge_ui(ns("dge"))),
                           tabPanel("4. Enrichment", rnaseq_ora_ui(ns("ora")))
                         )))
      )
    )
  )
}


rnaseqServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      data <- rnaseq_upload_server("upload")
      qc <- rnaseq_qc_server("qc", data$counts, data$coldata)
      res <- rnaseq_dge_server("dge", qc$dds)
      rnaseq_ora_server("ora", res$results, res$results_sig, qc$mat)
    }
  )
}

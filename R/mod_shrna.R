# shinyBioTools - shRNA Module

library(httr)
library(rvest)
library(Biostrings)

# =============================================================================
# splashRNA Helper Functions
# =============================================================================

#' Query splashRNA Database for shRNA Design
#'
#' @description Queries the splashRNA web server (http://splashrna.mskcc.org)
#'   to retrieve optimized shRNA sequences for a given Entrez Gene ID.
#'
#' @param id Character. A single Entrez Gene ID.
#' @param n Integer. Number of shRNA predictions to retrieve (1-6).
#' @param anno Character. Annotation prefix for result IDs.
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item ID - shRNA identifier (e.g., "Gene#1")
#'     \item Antisense - Antisense sequence for shRNA
#'     \item Score - Quality score from splashRNA
#'   }
#'
#' @examples
#' \dontrun{
#' splashRNA("7157", n = 3, anno = "TP53")
#' }
#'
#' @export
splashRNA <- function(id, n = 3, anno = "Gene") {
  n <- max(1, min(6, n))

  resp <- tryCatch(
    POST(
      url = "http://splashrna.mskcc.org/show_results",
      body = list(fasta = paste0("> entrezID\n", id),
                  n_predictions = as.character(n), removeRE = "on",
                  apa = "on", academic = "on", email = "user@example.com"),
      add_headers(Origin = "http://splashrna.mskcc.org",
                  Referer = "http://splashrna.mskcc.org/"),
      encode = "form",
      timeout(30)
    ),
    error = function(e) stop("splashRNA request failed: ", conditionMessage(e))
  )

  if (http_error(resp)) {
    stop("splashRNA returned HTTP status ", resp$status_code,
         ". The service may be unavailable.")
  }

  html <- read_html(resp)
  base_xpath <- "/html/body/div[2]/div[2]/div/table/tbody/tr["

  cell_text <- function(row, col) {
    node <- html_node(html, xpath = paste0(base_xpath, row, "]/td[", col, "]"))
    if (is.null(node)) NA_character_ else html_text(node, trim = TRUE)
  }

  results <- lapply(seq_len(n), function(i) {
    data.frame(
      ID = paste0(anno, "#", i),
      Antisense = cell_text(i, 2),
      Score = cell_text(i, 3),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, results)
  # drop predictions whose table cells were missing
  out[!is.na(out$Antisense), , drop = FALSE]
}

#' Batch Query splashRNA Database
#'
#' @description Queries the splashRNA web server for multiple genes sequentially
#'   with a delay between requests to avoid overloading the server.
#'
#' @param IDs Character vector of Entrez Gene IDs.
#' @param N Integer. Number of shRNA predictions per gene (1-6).
#' @param Anno Character vector of annotation prefixes for each gene.
#'   If NULL, auto-generates "Gene_1", "Gene_2", etc.
#' @param on_gene Optional callback \code{function(i, total)} invoked after each
#'   gene query (e.g., for progress bars).
#'
#' @return A data.frame with the same structure as \code{\link{splashRNA}},
#'   containing results for all genes combined. Failed queries return empty
#'   rows with warnings.
#'
#' @examples
#' \dontrun{
#' splashRNA_batch(c("7157", "672"), N = 3, Anno = c("TP53", "BRCA1"))
#' }
#'
#' @export
splashRNA_batch <- function(IDs, N = 3, Anno = NULL, on_gene = NULL) {
  if (is.null(Anno)) Anno <- paste0("Gene_", seq_along(IDs))

  results <- lapply(seq_along(IDs), function(i) {
    message("Processing gene ", i, "/", length(IDs), ": ", IDs[i])
    if (i < length(IDs)) Sys.sleep(2)
    res <- tryCatch(
      splashRNA(IDs[i], N, Anno[i]),
      error = function(e) {
        warning("Failed: ", IDs[i], " - ", e$message)
        data.frame(ID = character(0), Antisense = character(0),
                   Score = character(0), stringsAsFactors = FALSE)
      }
    )
    if (is.function(on_gene)) on_gene(i, length(IDs))
    res
  })
  do.call(rbind, results)
}

#' Design shRNA Primers for miR-30 Based Vectors
#'
#' @description Generates forward and reverse primers for cloning shRNA
#'   sequences into miR-30 based expression vectors using overlapping PCR.
#'
#' @param splashRNA_result Data frame with columns:
#'   \itemize{
#'     \item ID - Identifier for each shRNA
#'     \item Antisense - Antisense sequence from splashRNA
#'     \item Score - (Optional) Quality score
#'   }
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item ID - Primer identifier (e.g., "Gene#1-F", "Gene#1-R")
#'     \item Primer - Primer sequence
#'   }
#'
#' @examples
#' \dontrun{
#' result <- data.frame(ID = "TP53#1", Antisense = "ATGC...", stringsAsFactors = FALSE)
#' shRNA_primer(result)
#' }
#'
#' @export
shRNA_primer <- function(splashRNA_result) {
  loop <- DNAString("TAGTGAAGCCACAGATGTA")
  ovlp5 <- DNAString("TGCTGTTGACAGTGAGCG")
  ovlp3 <- DNAString("TGCCTACTGCCTCGGACT")

  primers <- lapply(seq_len(nrow(splashRNA_result)), function(i) {
    anti_seq <- splashRNA_result$Antisense[i]
    if (is.na(anti_seq) || !grepl("^[ATCG]+$", anti_seq) || nchar(anti_seq) != 22) {
      warning("Skipping ", splashRNA_result$ID[i],
              ": antisense sequence is not a valid 22-mer (got ",
              if (is.na(anti_seq)) "NA" else paste0(nchar(anti_seq), " nt"), ")")
      return(NULL)
    }
    anti <- DNAString(anti_seq)
    mm_base <- switch(as.character(anti[22]),
                      "A" = "C", "G" = "A", "C" = "A", "T" = "C")
    sense <- c(DNAString(mm_base), reverseComplement(anti[-22]))

    primer_F <- as.character(paste0(loop, anti, ovlp3))
    primer_R <- as.character(reverseComplement(
      DNAString(paste0(ovlp5, sense, loop))))

    id <- splashRNA_result$ID[i]
    data.frame(ID = c(paste0(id, "-F"), paste0(id, "-R")),
               Primer = c(primer_F, primer_R), stringsAsFactors = FALSE)
  })
  do.call(rbind, primers)
}

# =============================================================================
# UI Function
# =============================================================================

shRNAUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(title = "Settings", width = 3, status = "primary", solidHeader = TRUE,
      radioButtons(ns("species"), "Species", inline = TRUE,
                   choices = c("Human" = "hs", "Mouse" = "mm")),
      radioButtons(ns("type"), "Input Type",
                   choices = c("Symbol" = "SYMBOL", "Ensembl" = "ENSEMBL",
                               "Entrez" = "ENTREZID")),
      numericInput(ns("number"), "shRNA per gene", 3, 1, 6),
      textAreaInput(ns("input"), "Gene List",
                    placeholder = "One gene per line", height = "150px"),
      helpText("ID conversion runs locally; primer design queries splashRNA online.",
               hint_tip("Gene IDs are converted with the local org.Hs/org.Mm annotation database. 'Design Primers' then queries the splashRNA web service (splashrna.mskcc.org) once per gene, a few seconds apart - this needs an internet connection and can take a minute for long gene lists. Primers are for miR-30-based vectors (22-mer antisense + loop + overlap arms).")),
      actionButton(ns("submit"), "Convert IDs", icon = icon("sync"),
                   class = "btn-default",
                   title = "Map the gene list to SYMBOL / ENTREZID / ENSEMBL."),
      actionButton(ns("design"), "Design Primers", icon = icon("magic"),
                   class = "btn-default",
                   title = "Fetch the top shRNA predictions per gene and build miR-30 cloning primers (online)."),
      downloadButton(ns("download"), "Results (xlsx)", icon = icon("table"),
                     title = "Converted gene IDs, antisense sequences and primer pairs, one sheet each.")
    ),
    column(9,
      box(title = "Converted Gene IDs", width = NULL, status = "primary", solidHeader = TRUE, collapsible = TRUE,
        DT::dataTableOutput(ns("gene_table"))
      ),
      box(title = "shRNA Antisense Sequences", width = NULL, status = "primary", solidHeader = TRUE, collapsible = TRUE,
        DT::dataTableOutput(ns("antisense_table"))
      ),
      box(title = "Designed Primers", width = NULL, status = "primary", solidHeader = TRUE, collapsible = TRUE,
        DT::dataTableOutput(ns("primer_table"))
      )
    )
  )
}

shRNAServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    rv <- reactiveValues(gene_info = NULL, antisense = NULL, primer = NULL)
    
    # Convert gene IDs
    observeEvent(input$submit, {
      req(input$input)
      genes <- trimws(strsplit(input$input, "\n")[[1]])
      genes <- genes[genes != ""]

      if (length(genes) == 0) {
        showNotification("Enter at least one gene", type = "warning")
        return()
      }

      rv$gene_info <- id_convert(genes, input$type, input$species)
    })
    
    # Design shRNA primers
    observeEvent(input$design, {
      req(rv$gene_info)
      
      valid <- rv$gene_info[!is.na(rv$gene_info$ENTREZID), ]
      entrez_ids <- unique(valid$ENTREZID)
      
      if (length(entrez_ids) == 0) {
        showNotification("No valid Entrez IDs found", type = "warning")
        return()
      }
      
      symbols <- valid[!duplicated(valid$ENTREZID), c("SYMBOL", "ENTREZID")]

      withProgress(message = "Querying splashRNA...", value = 0, {
        rv$antisense <- splashRNA_batch(
          entrez_ids, input$number,
          symbols$SYMBOL[match(entrez_ids, symbols$ENTREZID)],
          on_gene = function(i, total) incProgress(1 / total)
        )

        rv$primer <- if (!is.null(rv$antisense) && nrow(rv$antisense) > 0) {
          tryCatch(shRNA_primer(rv$antisense),
                   error = function(e) {
                     showNotification(paste("Primer design failed:", e$message),
                                      type = "error", duration = 10)
                     NULL
                   })
        } else {
          NULL
        }

        if (is.null(rv$antisense) || nrow(rv$antisense) == 0) {
          showNotification("No shRNA predictions were retrieved from splashRNA",
                           type = "warning", duration = 10)
        }
      })
    })
    
    # Outputs
    output$gene_table <- DT::renderDataTable({
      req(rv$gene_info)
      DT::datatable(rv$gene_info[, c("query", "SYMBOL", "ENTREZID", "ENSEMBL")],
                    rownames = FALSE, options = dt_table_options())
    })

    output$antisense_table <- DT::renderDataTable({
      req(rv$antisense)
      DT::datatable(rv$antisense, rownames = FALSE, options = dt_table_options())
    })

    output$primer_table <- DT::renderDataTable({
      req(rv$primer)
      DT::datatable(rv$primer, rownames = FALSE, options = dt_table_options()) %>%
        DT::formatStyle(columns = "Primer", fontFamily = "Consolas, monospace")
    })

    output$download <- downloadHandler(
      filename = function() paste0("shRNA_", Sys.Date(), "_results.xlsx"),
      content = function(file) {
        sheets <- list(
          "Gene IDs" = rv$gene_info,
          "Antisense" = rv$antisense,
          "Primers" = rv$primer
        )
        write_xlsx_sheets(sheets, file)
      }
    )
  })
}

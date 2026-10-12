# shinyLab - SynergyFinder Module
# Drug synergy analysis using the synergyfinder Bioconductor package
#
# Uploads a formatted CSV, calculates ZIP/HSA/Bliss/Loewe synergy plus
# sensitivity, previews results as static charts, and exports
# organized data (.xlsx) and vector plots (.zip of PDFs).

library(shiny)
library(shinydashboard)
library(tidyverse)
library(magrittr)
library(patchwork)
library(synergyfinder)
library(writexl)
library(DT)
library(zip)


# =============================================================================
# Synergy Calculation Function
# =============================================================================

run_synergyfinder <- function(data,
                              type = c('viability', 'inhibition'),
                              save_dir = '.') {

  type <- type[1]
  if (!(type %in% c('viability', 'inhibition'))) {
    stop("type must be 'viability' or 'inhibition'")
  }

  # Remove NA rows
  data <- drop_na(data, block_id)

  # Reshape data
  data <- ReshapeData(
    data = data,
    data_type = type,
    impute = TRUE,
    impute_method = NULL,
    noise = TRUE,
    seed = 1
  )

  # Analyze synergy
  res <- CalculateSynergy(
    data = data,
    method = c("ZIP", "HSA", "Bliss", "Loewe"),
    Emin = NA,
    Emax = NA,
    correct_baseline = "non"
  ) %>%
    CalculateSensitivity()

  # Save results
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  blocks <- unique(res$drug_pairs$block_id)

  for (block in blocks) {
    # Dose-response (draws on the open device - returns a recorded plot)
    pdf(file.path(save_dir, str_glue('DoseResponse_{block}.pdf')), width = 12,
        height = 6)
    tryCatch(
      PlotDoseResponse(data = res, block_ids = block, drugs = c(1, 2),
                       save_file = FALSE),
      error = function(e) NULL
    )
    dev.off()

    # Synergy maps, one page per method (display = FALSE returns ggplot objects)
    drug1 <- res$drug_pairs %>% filter(block_id == block) %>% pull(drug1) %>% .[1]
    drug2 <- res$drug_pairs %>% filter(block_id == block) %>% pull(drug2) %>% .[1]

    p1 <- map(c("ZIP", "HSA", "Bliss", "Loewe"), ~ PlotSynergy(
      data = res,
      type = "2D",
      method = .x,
      block_ids = block,
      drugs = c(1, 2),
      dynamic = FALSE,
      display = FALSE,
      save_file = FALSE
    ))
    pdf(file.path(save_dir, str_glue('{drug1}_{drug2}_Synergy_{block}.pdf')),
        width = 6, height = 5)
    invisible(map(p1, print))
    dev.off()

    # 2D synergy surfaces (static plots draw to the open device directly)
    pdf(file.path(save_dir, str_glue('{drug1}_{drug2}_Surface_{block}.pdf')),
        width = 6, height = 6)
    tryCatch(
      invisible(pmap(list(
        c("ZIP_synergy", "HSA_synergy", "Bliss_synergy", "Loewe_synergy"),
        c("ZIP", "HSA", "Bliss", "Loewe")
      ), function(value, method) {
        Plot2DrugSurface(
          data = res,
          plot_block = block,
          drugs = c(1, 2),
          plot_value = value,
          summary_statistic = 'mean',
          dynamic = FALSE,
          interpolate_len = 3
        )
      })),
      error = function(e) NULL
    )
    dev.off()

    # Synergy barometer (skip when IC50 fitting failed for either drug)
    ic50_1 <- res$drug_pairs %>% filter(block_id == block) %>%
      pull(ic50_1) %>% as.numeric()
    conc_1 <- filter(res$response, block_id == block) %>%
      pull(conc1) %>% unique() %>% as.numeric()
    c_1 <- conc_1[which(min_rank(abs(conc_1 - ic50_1)) == 1)]

    ic50_2 <- res$drug_pairs %>% filter(block_id == block) %>% pull(ic50_2) %>% as.numeric()
    conc_2 <- filter(res$response, block_id == block) %>%
      pull(conc2) %>% unique() %>% as.numeric()
    c_2 <- conc_2[which(min_rank(abs(conc_2 - ic50_2)) == 1)]

    if (!is.na(ic50_1) && !is.na(ic50_2) &&
        length(c_1) == 1 && length(c_2) == 1) {
      pdf(file.path(save_dir, str_glue('{drug1}_{drug2}_Barometer_{block}.pdf')),
          width = 6, height = 6)
      tryCatch({
        p3 <- PlotBarometer(
          data = res,
          plot_block = block,
          plot_concs = c(c_1, c_2),
          needle_text_offset = -2.5
        )
        if (inherits(p3, "ggplot")) print(p3)
      }, error = function(e) NULL)
      dev.off()
    }
  }

  # Sensitivity and synergy plot (draws to the open device directly)
  pdf(file.path(save_dir, 'synergy_sensitivity.pdf'), width = 16, height = 4)
  tryCatch(
    invisible(map(c("ZIP", "HSA", "Bliss", "Loewe"), ~ PlotSensitivitySynergy(
      data = res,
      plot_synergy = .x,
      point_size = 3,
      label_size = 8,
      show_labels = TRUE,
      dynamic = FALSE
    ))),
    error = function(e) NULL
  )
  dev.off()

  writexl::write_xlsx(
    list(
      "Synergy Scores" = res$synergy_scores,
      "Drug Pairs & Sensitivity" = res$drug_pairs
    ),
    file.path(save_dir, 'result_synergyfinder.xlsx')
  )

  return(res)
}

# Static preview wrappers. Plot2DrugSurface and PlotSensitivitySynergy draw
# onto the open graphics device (renderPlot captures that), PlotSynergy
# returns a ggplot. NULL is returned when a chart fails so the app can show a
# placeholder instead of crashing.

preview_plot <- function(expr) {
  tryCatch(expr, error = function(e) {
    showNotification(paste("Preview unavailable:", conditionMessage(e)),
                     type = "warning", duration = 6)
    NULL
  })
}

synergy_surface_plot <- function(res, block) {
  preview_plot(Plot2DrugSurface(
    data = res, plot_block = block, drugs = c(1, 2),
    plot_value = "response", summary_statistic = "mean",
    dynamic = FALSE, interpolate_len = 5
  ))
  invisible(NULL)
}

synergy_method_plot <- function(res, block, method) {
  p <- preview_plot({
    p <- PlotSynergy(
      data = res, type = "heatmap", method = method, block_ids = block,
      drugs = c(1, 2), dynamic = FALSE, display = FALSE, save_file = FALSE
    )
    if (is.list(p)) p <- p[[1]]
    p
  })
  if (!is.null(p)) print(p)
  invisible(NULL)
}

synergy_sensitivity_plot <- function(res, method) {
  preview_plot(PlotSensitivitySynergy(
    data = res, plot_synergy = method, point_size = 3, label_size = 8,
    show_labels = TRUE, dynamic = FALSE
  ))
  invisible(NULL)
}


# =============================================================================
# UI Function
# =============================================================================

synergyUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(title = "Settings", width = 3, status = "primary", solidHeader = TRUE,
      fileInput(ns("input_file"), "Upload dose-response CSV", accept = ".csv"),
      helpText("CSV with block_id, drug1, drug2, conc1, conc2, response.",
               hint_tip("Expected columns: block_id, drug1, drug2, conc_unit1, conc_unit2, conc1, conc2, response. Each block needs a 0-dose row for both drugs (the models require the untreated baseline); missing dose combinations are imputed. Response is % viability (or % inhibition - pick the matching Effect type).")),
      radioButtons(ns("type"), "Effect type",
                   choices = c("viability", "inhibition"),
                   selected = "viability", inline = FALSE),
      div(style = "display: flex; justify-content: space-between; align-items: center; gap: 4px; flex-wrap: nowrap;",
        actionButton(ns("submit"), "Analyze", icon = icon("play"), class = "btn-default btn-sm",
                     title = "Reshape the data, fit the ZIP / HSA / Bliss / Loewe models and compute sensitivity (IC50, RI, CSS) per block."),
        div(style = "display: flex; gap: 4px; flex-wrap: nowrap;",
          downloadButton(ns("download_data"), "Data", icon = icon("table"), class = "btn-default btn-sm",
                         title = "Synergy scores per dose, drug-pair fits (IC50 / RI / CSS), tidy response and the raw input as one .xlsx, one sheet each."),
          downloadButton(ns("download_plots"), "Plots", icon = icon("file-pdf"), class = "btn-default btn-sm",
                         title = "All vector PDFs generated by synergyfinder: dose response, synergy maps, 2D surfaces, barometers and sensitivity plots, zipped.")))
    ),
    box(title = "Results", width = 9, status = "primary", solidHeader = TRUE,
      uiOutput(ns("results_ui"))
    )
  )
}


# =============================================================================
# Server Function
# =============================================================================

synergyServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    rv <- reactiveValues(tmpdir = NULL, res = NULL, input_data = NULL)

    # =========================================================================
    # Calculate synergy on submit
    # =========================================================================
    observeEvent(input$submit, {
      req(input$input_file)

      raw <- tryCatch(read_csv(input$input_file$datapath, show_col_types = FALSE),
                      error = function(e) {
                        showNotification(paste("Could not read CSV:", e$message),
                                         type = "error", duration = 8)
                        NULL
                      })
      req(raw)

      # Fresh per-submit directory; remove the previous one
      if (!is.null(rv$tmpdir)) unlink(rv$tmpdir, recursive = TRUE)
      tmpdir <- tempfile(pattern = "synergy_")
      dir.create(tmpdir, recursive = TRUE)
      rv$tmpdir <- tmpdir

      res <- tryCatch(
        withProgress(
          run_synergyfinder(raw, type = input$type, save_dir = tmpdir),
          message = 'Calculating synergy ...', value = 0.4
        ),
        error = function(e) {
          showNotification(paste("Synergy calculation failed:", e$message),
                           type = "error", duration = 10)
          NULL
        }
      )
      req(res)

      rv$res <- res
      rv$input_data <- raw
      showNotification("Synergy calculation finished.", type = "message", duration = 5)
    })

    session$onSessionEnded(function() {
      tmp <- isolate(rv$tmpdir)
      if (!is.null(tmp)) unlink(tmp, recursive = TRUE)
    })

    res <- reactive(req(rv$res))

    # =========================================================================
    # Results UI
    # =========================================================================
    output$results_ui <- renderUI({
      req(res())
      blocks <- sort(unique(res()$drug_pairs$block_id))

      tabsetPanel(type = "tabs",
        tabPanel("Plots", br(),
          fluidRow(
            column(3, selectInput(ns("block"), tagList("Block",
                                    hint_tip("A block is one drug-pair x dose-matrix experiment; multi-block CSVs produce one entry each.")),
                                  choices = blocks)),
            column(3, selectInput(ns("method"), tagList("Method",
                                    hint_tip("Reference model for the synergy score: ZIP, HSA, Bliss or Loewe. Scores > 0 suggest synergy, < 0 antagonism.")),
                                  choices = c("ZIP", "HSA", "Bliss", "Loewe"),
                                  selected = "ZIP"))
          ),
          fluidRow(
            column(6, h5("Dose-response surface"), plotOutput(ns("dose_plot"), height = "420px")),
            column(6, h5("Synergy heatmap"), plotOutput(ns("synergy_heatmap"), height = "420px"))
          ),
          h5("Sensitivity vs synergy (all blocks)"),
          plotOutput(ns("sensitivity_plot"), height = "380px")
        ),
        tabPanel("Synergy Scores", br(),
          h5("Per-dose synergy scores (ZIP / HSA / Bliss / Loewe)"),
          DT::dataTableOutput(ns("synergy_tbl"))
        ),
        tabPanel("Drug Pairs & Sensitivity", br(),
          h5("Model fits, IC50 / RI / CSS per drug pair"),
          DT::dataTableOutput(ns("drug_pairs_tbl"))
        ),
        tabPanel("Files", br(),
          h5("PDF files included in the plots download"),
          uiOutput(ns("file_list"))
        )
      )
    })

    output$dose_plot <- renderPlot({ synergy_surface_plot(res(), input$block) }, res = 96)
    output$synergy_heatmap <- renderPlot({ synergy_method_plot(res(), input$block, input$method) }, res = 96)
    output$sensitivity_plot <- renderPlot({ synergy_sensitivity_plot(res(), input$method) }, res = 96)

    output$synergy_tbl <- DT::renderDataTable({
      d <- res()$synergy_scores
      DT::datatable(d, rownames = FALSE, options = dt_table_options()) %>%
        DT::formatRound(columns = c("conc1", "conc2", "ZIP_synergy", "HSA_synergy",
                                    "Bliss_synergy", "Loewe_synergy"), digits = 3)
    })

    output$drug_pairs_tbl <- DT::renderDataTable({
      d <- res()$drug_pairs %>% mutate(across(where(is.numeric), ~ round(.x, 4)))
      DT::datatable(d, rownames = FALSE, options = dt_table_options())
    })

    output$file_list <- renderUI({
      req(rv$tmpdir)
      files <- list.files(rv$tmpdir, pattern = "\\.pdf$")
      if (length(files) == 0) return(helpText("No plot files generated."))
      tags$ul(lapply(sort(files), tags$li))
    })

    # =========================================================================
    # Downloads
    # =========================================================================
    output$download_data <- downloadHandler(
      filename = function() paste0('Synergy_', Sys.Date(), '_data.xlsx'),
      content = function(file) {
        sheets <- list(
          "Synergy Scores" = res()$synergy_scores,
          "Drug Pairs & Sensitivity" = res()$drug_pairs,
          "Response (Tidy)" = res()$response,
          "Input Data" = rv$input_data
        )
        write_xlsx_sheets(sheets, file)
      }
    )

    output$download_plots <- downloadHandler(
      filename = function() paste0('Synergy_', Sys.Date(), '_plots.zip'),
      content = function(fname) {
        req(rv$tmpdir)
        pdf_files <- list.files(rv$tmpdir, pattern = "\\.pdf$", full.names = FALSE)
        req(length(pdf_files) > 0)
        zip::zip(zipfile = fname, files = pdf_files, root = rv$tmpdir)
      },
      contentType = "application/zip"
    )
  })
}

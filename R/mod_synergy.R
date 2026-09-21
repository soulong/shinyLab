# shinyLab - SynergyFinder Module
# Drug synergy analysis using synergyfinder Bioconductor package

library(shiny)
library(shinydashboard)
library(tidyverse)
library(magrittr)
library(patchwork)
library(synergyfinder)
library(writexl)


# =============================================================================
# UI Function
# =============================================================================

synergyUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(title = "Settings", width = 3, status = "primary", solidHeader = TRUE,
      fileInput(ns("input_file"), "Upload formatted CSV file", accept = ".csv"),
      radioButtons(ns("type"), "Effect type", 
                   choices = c("viability", "inhibition"), 
                   selected = "viability", inline = FALSE),
      hr(),
      actionButton(ns("submit"), "Submit", icon = icon("play"), class = "btn-default"),
      br(), br(),
      downloadButton(ns("download"), "Download", icon = icon("file-download"))
    ),
    box(title = "Results", width = 9, status = "primary", solidHeader = TRUE,
      verbatimTextOutput(ns("res_show"))
    )
  )
}


# =============================================================================
# Server Function
# =============================================================================

synergyServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    rv <- reactiveValues(tmpdir = NULL)

    # =========================================================================
    # Calculate synergy on submit
    # =========================================================================
    data <- eventReactive(input$submit, {
      req(input$input_file)

      # Fresh per-submit directory; remove the previous one
      if (!is.null(rv$tmpdir)) unlink(rv$tmpdir, recursive = TRUE)
      tmpdir <- tempfile(pattern = "synergy_")
      dir.create(tmpdir, recursive = TRUE)
      rv$tmpdir <- tmpdir

      withProgress(
        read_csv(input$input_file$datapath) %>%
          run_synergyfinder(type = input$type, save_dir = tmpdir),
        message = 'Calculating synergy ...', value = 0.4
      )

      res_files <- list.files(tmpdir, pattern = '(\\.pdf)|(\\.xlsx)',
                             full.names = FALSE)

      return(list(tmpdir = tmpdir, res_files = res_files))
    })

    session$onSessionEnded(function() {
      tmp <- isolate(rv$tmpdir)
      if (!is.null(tmp)) unlink(tmp, recursive = TRUE)
    })

    # =========================================================================
    # Display results
    # =========================================================================
    output$res_show <- renderPrint({
      d <- data()
      cat("Synergy calculation finished.\n\nOutput files:\n")
      print(d$res_files)
    })

    # =========================================================================
    # Download handler
    # =========================================================================
    output$download <- downloadHandler(
      filename = 'synergy_result.zip',
      content = function(fname) {
        req(data())
        zip::zip(zipfile = fname, files = data()$res_files, root = data()$tmpdir)
      },
      contentType = "application/zip"
    )
  })
}


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

  writexl::write_xlsx(res, file.path(save_dir, 'result_synergyfinder.xlsx'))

  blocks <- unique(res$drug_pairs$block_id)

  for (block in blocks) {
    # Dose-response matrix (synergyfinder 3.2 draws on call - capture into PDF)
    pdf(file.path(save_dir, str_glue('DoseResponse_{block}.pdf')), width = 12,
        height = 6)
    p_dr <- tryCatch(
      PlotDoseResponse(data = res, block_ids = block, drugs = c(1, 2),
                       save_file = FALSE),
      error = function(e) NULL
    )
    dev.off()

    # Plot synergy (display = FALSE returns ggplot objects we can save)
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

    # Plot 2D surface (static surface plots draw to the open device directly)
    pdf(file.path(save_dir, str_glue('{drug1}_{drug2}_Surface_{block}.pdf')),
        width = 6, height = 6)
    invisible(pmap(list(
      c("ZIP_synergy", "HSA_synergy", "Bliss_synergy", "Loewe_synergy"),
      c("ZIP", "HSA", "Bliss", "Loewe")
    ), function(value, method) {
      Plot2DrugSurface(
        data = res,
        plot_block = block,
        plot_value = value,
        summary_statistic = 'mean',
        dynamic = FALSE,
        interpolate_len = 3
      )
    }))
    dev.off()

    # Synergy barometer (skip when IC50 fitting failed for either drug)
    ic50_1 <- res$drug_pairs %>% filter(block_id == block) %>%
      pull(ic50_1) %>% as.numeric()
    conc_1 <- filter(res$response, block_id == block) %>%
      pull(conc1) %>% unique()
    c_1 <- conc_1[which(min_rank(abs(conc_1 - ic50_1)) == 1)]

    ic50_2 <- res$drug_pairs %>% filter(block_id == block) %>% pull(ic50_2) %>% as.numeric()
    conc_2 <- filter(res$response, block_id == block) %>%
      pull(conc2) %>% unique() %>% as.numeric()
    c_2 <- conc_2[which(min_rank(abs(conc_2 - ic50_2)) == 1)]

    if (!is.na(ic50_1) && !is.na(ic50_2) &&
        length(c_1) == 1 && length(c_2) == 1) {
      pdf(file.path(save_dir, str_glue('{drug1}_{drug2}_Barometer_{block}.pdf')),
          width = 6, height = 6)
      p3 <- PlotBarometer(
        data = res,
        plot_block = block,
        plot_concs = c(c_1, c_2),
        needle_text_offset = -2.5
      )
      if (inherits(p3, "ggplot")) print(p3)
      dev.off()
    }
  }

  # Sensitivity and synergy plot (3.2 draws to the open device directly)
  pdf(file.path(save_dir, 'synergy_sensitivity.pdf'), width = 16, height = 4)
  invisible(map(c("ZIP", "HSA", "Bliss", "Loewe"), ~ PlotSensitivitySynergy(
    data = res,
    plot_synergy = .x,
    point_size = 3,
    label_size = 8,
    show_labels = TRUE,
    dynamic = FALSE
  )))
  dev.off()

  return(res)
}

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
      actionButton(ns("submit"), "Submit", icon = icon("play"), class = "btn-primary"),
      br(), br(),
      downloadButton(ns("download"), "Download", icon = icon("file-download"))
    ),
    box(title = "Results", width = 9, status = "primary",
      tabsetPanel(
        tabPanel("Output",
          verbatimTextOutput(ns("res_show"))
        ),
        tabPanel("Parameters",
          fluidRow(
            column(6,
              sliderInput(ns("w_plot"), "Plot Width", 300, 1200, 700, 50),
              sliderInput(ns("h_plot"), "Plot Height", 200, 800, 400, 50)
            ),
            column(6,
              sliderInput(ns("sz_label"), "Label Size", 8, 20, 12, 1),
              sliderInput(ns("ang_x"), "X Axis Angle", 0, 90, 45, 15)
            )
          )
        )
      )
    )
  )
}


# =============================================================================
# Server Function
# =============================================================================

synergyServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    rv <- reactiveValues(data = NULL, tmpdir = NULL, res_files = NULL)
    
    # =========================================================================
    # Calculate synergy on submit
    # =========================================================================
    data <- eventReactive(input$submit, {
      req(input$input_file)
      
      # Create temp directory for results
      tmpdir <- str_glue("{tempdir()}_synergy_download")
      if (dir.exists(tmpdir)) {
        files <- list.files(tmpdir, all.files = TRUE, full.names = TRUE, 
                           recursive = TRUE, include.dirs = TRUE)
        unlink(files, recursive = TRUE)
      }
      
      withProgress(
        read_csv(input$input_file$datapath) %>%
          run_synergyfinder(type = input$type, save_dir = tmpdir),
        message = 'Calculating synergy ...', value = 0.4
      )
      
      res_files <- list.files(tmpdir, pattern = '(\\.pdf)|(\\.xlsx)', 
                             full.names = FALSE)
      
      rv$tmpdir <- tmpdir
      rv$res_files <- res_files
      
      return(list(tmpdir = tmpdir, res_files = res_files))
    })
    
    # =========================================================================
    # Display results
    # =========================================================================
    output$res_show <- renderPrint({
      data()
    })
    
    # =========================================================================
    # Download handler
    # =========================================================================
    output$download <- downloadHandler(
      filename = 'synergy_result.zip',
      content = function(fname) {
        req(data())
        old_wd <- getwd()
        setwd(data()$tmpdir)
        on.exit(setwd(old_wd))
        zip::zip(zipfile = fname, files = data()$res_files)
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
                              save_file = 'result_synergyfinder.xlsx',
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
    correct_baseline = "none"
  ) %>%
    CalculateSensitivity()
  
  # Save results
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  if (!is.null(save_file)) {
    writexl::write_xlsx(res, file.path(save_dir, save_file))
    
    # Plot dose response
    suppressWarnings(
      PlotDoseResponse(
        data = res,
        block_ids = unique(res$drug_pairs$block_id),
        drugs = c(1, 2),
        save_file = TRUE,
        file_type = "pdf",
        save_dir = save_dir
      )
    )
    
    # Generate plots for each block
    for (block in unique(res$drug_pairs$block_id)) {
      # Plot synergy
      p1 <- map(c("ZIP", "HSA", "Bliss", "Loewe"),
                ~ PlotSynergy(
                  data = res,
                  type = "2D",
                  method = .x,
                  block_ids = block,
                  drugs = c(1, 2),
                  dynamic = FALSE,
                  save_file = FALSE
                ))
      
      drug1 <- res$drug_pairs %>% filter(block_id == block) %>% pull(drug1) %>% .[1]
      drug2 <- res$drug_pairs %>% filter(block_id == block) %>% pull(drug2) %>% .[1]
      
      pdf(file.path(save_dir, str_glue('{drug1}_{drug2}_Synergy_{block}.pdf')), 
          width = 6, height = 5)
      print(p1)
      dev.off()
      
      # Plot 2D surface
      p2 <- map(c("ZIP_synergy", "HSA_synergy", "Bliss_synergy", "Loewe_synergy"),
                ~ Plot2DrugSurface(
                  data = res,
                  plot_block = block,
                  plot_value = .x,
                  summary_statistic = 'mean',
                  dynamic = FALSE,
                  interpolate_len = 3
                ))
      
      pdf(file.path(save_dir, str_glue('{drug1}_{drug2}_Surface_{block}.pdf')), 
          width = 6, height = 6)
      print(p2)
      dev.off()
      
      # Synergy barometer
      ic50_1 <- res$drug_pairs %>% filter(block_id == block) %>% 
        pull(ic50_1) %>% as.numeric()
      conc_1 <- filter(res$response, block_id == block) %>% 
        pull(conc1) %>% unique()
      c_1 <- conc_1[which(min_rank(abs(conc_1 - ic50_1)) == 1)]
      
      ic50_2 <- res$drug_pairs %>% filter(block_id == block) %>% pull(ic50_2)
      conc_2 <- filter(res$response, block_id == block) %>% 
        pull(conc2) %>% unique() %>% as.numeric()
      c_2 <- conc_2[which(min_rank(abs(conc_2 - ic50_2)) == 1)]
      
      p3 <- PlotBarometer(
        data = res,
        plot_block = block,
        plot_concs = c(c_1, c_2),
        needle_text_offset = -2.5
      )
      
      ggsave(file.path(save_dir, str_glue('{drug1}_{drug2}_Barometer_{block}.pdf')), 
             p3, height = 6, width = 6)
    }
    
    # Sensitivity and synergy plot
    p <- map(c("ZIP", "HSA", "Bliss", "Loewe"),
             ~ PlotSensitivitySynergy(
               data = res,
               plot_synergy = .x,
               point_size = 3,
               label_size = 8,
               show_labels = TRUE,
               dynamic = FALSE
             ))
    
    ggsave(file.path(save_dir, 'synergy_sensitivity.pdf'), 
           wrap_plots(p, nrow = 1),
           height = 4, width = 16)
  }
  
  return(res)
}

# shinyLab - DSF Module
# Differential Scanning Fluorimetry analysis module

library(patchwork)
library(mgcv)
library(gratia)
# library(rcdk)
library(rio)
library(DT)
library(zip)
library(janitor)
options(rio.import.class='tbl')

# =============================================================================
# Plot Helper Functions
# =============================================================================

build_scatter_plot <- function(df, x, y, color_col, facet_vars, 
                                xlim_min = NULL, xlim_max = NULL,
                                hline_y = NULL, ncol = 4) {
  use_color <- color_col != "" && color_col %in% names(df)
  
  if (use_color) {
    p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], color = .data[[color_col]])) +
      geom_point()
  } else {
    p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
      geom_point()
  }
  
  if (!is.null(hline_y)) {
    p <- p + geom_hline(yintercept = hline_y, linetype = "dashed", color = "gray50")
  }
  
  if (length(facet_vars) > 0 && !all(facet_vars == "")) {
    p <- p + facet_wrap(as.formula(paste("~", paste(facet_vars, collapse = " + "))), 
                        scales = "fixed", ncol = ncol)
  }
  
  p <- p + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  if (!is.null(xlim_min) || !is.null(xlim_max)) {
    p <- p + coord_cartesian(xlim = c(
      if (is.null(xlim_min)) NA_real_ else xlim_min,
      if (is.null(xlim_max)) NA_real_ else xlim_max
    ))
  }
  
  p
}


build_curve_plot <- function(data, ref_ligand_val, color_by, facet_vars,
                              is_derivative = FALSE, x_limits = NULL, ncol = 4) {
  has_ref <- !is.null(ref_ligand_val) && ref_ligand_val != "" && 
             ref_ligand_val %in% data$ligand
  is_valid_color <- !is.null(color_by) && color_by != "" && 
                    color_by %in% names(data) && color_by != "ligand"
  is_color_numeric <- if (is_valid_color) is.numeric(data[[color_by]]) else FALSE
  
  geom_fn <- if (is_derivative) geom_line else geom_path
  y_var <- if (is_derivative) "derivative" else "fluorescence"
  
  p <- ggplot()
  
  if (has_ref) {
    p <- p + geom_fn(
      data = filter(data, ligand == ref_ligand_val),
      aes(temperature, .data[[y_var]], group = well),
      color = 'grey50'
    )
  }
  
  filter_data <- if (has_ref) filter(data, ligand != ref_ligand_val) else data
  
  if (is_valid_color) {
    p <- p + geom_fn(
      data = filter_data,
      aes(temperature, .data[[y_var]], group = well, color = .data[[color_by]])
    )
    if (is_color_numeric) {
      p <- p + scale_color_viridis_c(option = 'C', na.value = 'darkred')
    } else {
      p <- p + scale_color_viridis_d(option = 'C', na.value = 'darkred')
    }
  } else {
    p <- p + geom_fn(
      data = filter_data,
      aes(temperature, .data[[y_var]], group = well),
      color = 'darkred'
    )
  }
  
  if (length(facet_vars) > 0 && !all(facet_vars == "")) {
    p <- p + facet_wrap(as.formula(paste("~", paste(facet_vars, collapse = " + "))), 
                        ncol = ncol)
  }
  
  p <- p + theme_bw() + theme(legend.position = 'top')
  
  if (!is.null(x_limits)) {
    p <- p + coord_cartesian(xlim = x_limits)
  }
  
  p
}


compute_derivatives <- function(mc_tidy) {
  mc_tidy %>%
    filter(!is.na(ligand), !is.na(target)) %>%
    nest(.by = c(plate, target, ligand, well)) %>%
    mutate(pred = purrr::map(data, \(d) {
      tryCatch({
        mod <- mgcv::gam(fluorescence ~ s(temperature), data = d)
        res <- gratia::derivatives(mod, n = 100, order = 1)
        tibble(temperature = res$temperature,
               derivative = res$.derivative)
      }, error = function(e) NULL)
    })) %>%
    unnest(pred) %>%
    filter(!is.na(derivative))
}

# =============================================================================
# UI Function
# =============================================================================

dsfUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(width = 3,
      box(title = "File Upload", width = 12, status = "primary", solidHeader = TRUE,
        fileInput(ns("file_analysis"), "Analysis Data Files (multiple)", 
                  accept = c(".txt", ".tsv"), multiple = TRUE),
        checkboxInput(ns("use_custom_metadata"), "Use custom metadata file", value = TRUE),
        conditionalPanel(
          condition = "input.use_custom_metadata == true",
          ns = ns,
          fileInput(ns("file_metadata"), "Plate meta (.xlsx)", accept = ".xlsx")
        ),
        radioButtons(ns("tm_type"), "Tm Calculation",
                     choices = c("Tm (Derivative)" = "tm_d", "Tm (Boltzmann)" = "tm_b"),
                     selected = "tm_d", inline = FALSE),
        hr(),
        actionButton(ns("submit"), "Submit", icon = icon("play"), class = "btn-primary"),
        br(), br(),
        uiOutput(ns("ref_selectors")),
        hr(),
        downloadButton(ns("download_data"), "Download Data", icon = icon("table")),
        br(), br(),
        downloadButton(ns("download_plots"), "Download Plots", icon = icon("file-pdf"))
      ),
      box(title = "Plot Settings", width = 12, status = "primary", solidHeader = TRUE,
        collapsible = TRUE,
        selectInput(ns("color_var"), "Color By", choices = "", selected = ""),
        selectInput(ns("facet_vars"), "Facet By (hold Ctrl for multiple)", 
                    choices = "", selected = "", multiple = TRUE),
        hr(),
        sliderInput(ns("plot_width"), "Plot Width", 4, 16, 7, 1),
        sliderInput(ns("plot_height"), "Plot Height", 4, 16, 5, 1),
        numericInput(ns("facet_ncol"), "Facet Columns", value = 4, min = 1, max = 12),
        hr(),
        numericInput(ns("xlim_min"), "X Min", value = NULL),
        numericInput(ns("xlim_max"), "X Max", value = NULL)
      )
    ),
    column(width = 9,
      box(title = "Results", width = 12, status = "primary",
        tabsetPanel(
          tabPanel("Data Preview",
            DT::dataTableOutput(ns("data_preview"))
          ),
          tabPanel("TM Scatter",
            uiOutput(ns("plot_tm_wrapper"))
          ),
          tabPanel("Delta TM Scatter",
            uiOutput(ns("plot_dtm_wrapper"))
          ),
          tabPanel("Raw Curves",
            uiOutput(ns("plot_raw_wrapper"))
          ),
          tabPanel("Derivative Curves",
            uiOutput(ns("plot_deri_wrapper"))
          ),
          tabPanel("Ligand Details",
            uiOutput(ns("ligand_selector")),
            uiOutput(ns("plot_ligand_detail_wrapper"))
          )
        )
      )
    )
  )
}


# =============================================================================
# Server Function
# =============================================================================

dsfServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rv <- reactiveValues(
      data = NULL,
      df_wide = NULL,
      mc_tidy = NULL,
      meta = NULL,
      has_ref_ligand = FALSE,
      ref_ligand_val = "",
      error_message = NULL
    )
    
    # =========================================================================
    # Process data on submit
    # =========================================================================
    data <- eventReactive(input$submit, {
      req(input$file_analysis)

      tryCatch({
        withProgress(message = 'Processing data...', value = 0, {
          
          incProgress(0.0, message = "Reading analysis results...")
          
          all_files <- input$file_analysis$datapath %>% 
            set_names(., nm = input$file_analysis$name)
          
          # Read analysis results
          res_file <- names(all_files) %>% 
            str_detect("AnalysisResults") %>% 
            which() %>% all_files[.]
          if(length(res_file) == 0) {
            stop("No AnalysisResults file found in the uploaded files")
          }
          if(length(res_file) > 1) {
            stop("Multiple AnalysisResults files detected. Please upload only one.")
          }
          res <- import(res_file) %>%
            janitor::clean_names() %>% 
            janitor::remove_empty('cols') %>% 
            dplyr::rename(plate = experiment_file_name) %>%
            mutate(tm = .data[[input$tm_type]],
                   well = transform_well_style(well),
                   plate = str_remove(plate, "\\.eds$") )
          
          incProgress(0.2, message = "Reading metadata...")
          # Conditionally read metadata
          if (input$use_custom_metadata) {
            req(input$file_metadata)
            meta_files <- input$file_metadata$datapath %>% 
              set_names(., str_remove(input$file_metadata$name, "\\.(eds\\.)?xlsx$"))
            meta <- meta_files %>% 
              map(read_metadata) %>% 
              list_rbind(., names_to="plate")
            if(!all(c("target", "ligand") %in% colnames(meta))) {
              stop("Custom metadata must contain 'target' and 'ligand' columns")
            }
            # inner_join: keep only (plate, well) present in BOTH res and meta
            merged <- res %>% 
              inner_join(meta, by = c("plate", "well"))
            
          } else {
            if (!all(c("target", "ligand") %in% colnames(res))) {
              stop("AnalysisResults must contain 'target' and 'ligand' columns when not using custom metadata")
            }
            merged <- res %>%
              dplyr::select(plate, well, target, ligand, tm,
                     any_of(c("conc", "smiles")))
            meta <- merged %>%
              dplyr::select(any_of(c("plate", "well", "target", "ligand", "conc", "smiles"))) %>%
              distinct()
          }

          # Explicit filter: keep only rows with non-NA ligand and target
          merged <- merged %>%
            filter(!is.na(ligand), !is.na(target)) %>%
            {if("conc" %in% colnames(.)) mutate(., conc = as.numeric(conc)) else .} %>% 
            relocate(tm, .after = last_col())

          incProgress(0.2, message = "Processing melting curve...")
          # Read raw data files with stricter plate name matching
          match_mc_files <- c()
          for(x in unique(res$plate)) {
            pattern <- paste0("RawData_", x, "\\.eds\\.txt$")
            match_idx <- which(str_detect(names(all_files), pattern))
            if (length(match_idx) == 0) {
              # Fallback: match plate name before .eds.txt
              match_idx <- which(str_detect(names(all_files), 
                paste0(x, "\\.eds\\.txt$")))
            }
            if(length(match_idx) >= 1) {
              matched <- all_files[match_idx[1]]
              names(matched) <- x
              match_mc_files <- c(match_mc_files, matched)
            }
          }
          
          if (length(match_mc_files) == 0) {
            stop("No matching raw data files found for any plate")
          }
          
          mc <- match_mc_files %>% 
            map(\(x) import(x)) %>%
            list_rbind(names_to = 'plate') %>%
            janitor::clean_names() %>%
            dplyr::select(!c(well)) %>%
            dplyr::rename(well = well_position) %>%
            mutate(well = transform_well_style(well)) %>% 
            mutate(temperature = as.numeric(temperature),
                   fluorescence = as.numeric(fluorescence))
          
          incProgress(0.25, detail = "Merging with raw data...")
          mc_tidy <- mc %>%
            inner_join(merged %>% dplyr::select(plate, well, target, ligand) %>% distinct(),
                      by = c("plate", "well"))
          
          incProgress(0.3, detail = "Computing derivatives...")
          derivatives <- compute_derivatives(mc_tidy)
          
          # Auto-detect reference values
          ligand_choices <- sort(unique(meta$ligand))
          ref_ligand_default <- find_best_match(ligand_choices, c("blank", "dmso", "pbs"))
          
          # Auto-compute d_tm with detected reference
          has_ref_ligand <- ref_ligand_default != "" && ref_ligand_default %in% merged$ligand
          
          if (has_ref_ligand) {
            ref <- merged %>%
              filter(ligand == ref_ligand_default) %>%
              reframe(tm_ref = median(tm, na.rm = TRUE), .by = c(plate, target))
            
            merged_with_dtm <- merged %>%
              left_join(ref, by = c("plate", "target")) %>%
              mutate(d_tm = tm - tm_ref) %>%
              dplyr::select(!tm_ref)
          } else {
            merged_with_dtm <- merged %>% mutate(d_tm = NA_real_)
          }
          
          incProgress(0.35, detail = "Done!")
          
          rv$meta <- meta
          rv$mc_tidy <- mc_tidy
          rv$derivatives <- derivatives
          rv$auto_ref_ligand <- ref_ligand_default
          rv$ref_ligand_val <- ref_ligand_default
          rv$has_ref_ligand <- has_ref_ligand
          rv$cached_df_with_dtm <- merged_with_dtm
          
          return(list(
            merged = merged_with_dtm,
            mc_tidy = mc_tidy, 
            meta = meta,
            derivatives = derivatives,
            ligand_choices = ligand_choices
          ))
        })
      }, error = function(e) {
        showNotification(paste("Error processing data:", e$message), 
                         type = "error", duration = 10)
        rv$error_message <- e$message
        return(NULL)
      })
    })
    
    # =========================================================================
    # Dynamic ref selectors (shown after submit)
    # =========================================================================
    output$ref_selectors <- renderUI({
      req(data())
      tagList(
        selectInput(ns("ref_ligand_input"), "Reference Ligand",
                    choices = data()$ligand_choices,
                    selected = rv$auto_ref_ligand)
      )
    })
    
    # =========================================================================
    # Calculate d_tm reactively when ref changes
    # =========================================================================
    df_with_dtm <- reactive({
      req(data(), input$ref_ligand_input)
      
      merged <- data()$merged
      ref_ligand_val <- input$ref_ligand_input
      
      if (is.null(ref_ligand_val) || ref_ligand_val == "") {
        rv$has_ref_ligand <- FALSE
        rv$ref_ligand_val <- ""
        return(merged %>% mutate(d_tm = NA_real_))
      }
      
      if (!ref_ligand_val %in% merged$ligand) {
        showNotification(paste("Reference ligand", ref_ligand_val, "not found in data"), 
                         type = "warning", duration = 5)
        rv$has_ref_ligand <- FALSE
        rv$ref_ligand_val <- ""
        return(merged %>% mutate(d_tm = NA_real_))
      }
      
      # Use cached result if ref hasn't changed from auto-detected default
      if (!is.null(rv$cached_df_with_dtm) && 
          ref_ligand_val == rv$auto_ref_ligand &&
          rv$ref_ligand_val == ref_ligand_val) {
        rv$has_ref_ligand <- TRUE
        return(rv$cached_df_with_dtm)
      }
      
      # Recompute with new reference
      ref <- merged %>%
        filter(ligand == ref_ligand_val) %>%
        reframe(tm_ref = median(tm, na.rm = TRUE), .by = c(plate, target))
      
      df <- merged %>%
        left_join(ref, by = c("plate", "target")) %>%
        mutate(d_tm = tm - tm_ref) %>%
        dplyr::select(!tm_ref)
      
      rv$has_ref_ligand <- TRUE
      rv$ref_ligand_val <- ref_ligand_val
      return(df)
    })
    
    # Create wide format with d_tm
    df_wide <- reactive({
      req(df_with_dtm())
      
      df <- df_with_dtm()
      meta <- data()$meta
      
      df %>%
        pivot_wider(id_cols = any_of(c("plate", setdiff(colnames(meta), c("well", "target")))),
                    names_from = target, values_from = d_tm,
                    names_prefix = 'dTM_',
                    values_fn = median) %>%
        mutate(across(where(is.numeric), \(x) round(x, 3)))
    })
    
    # =========================================================================
    # Update select inputs after data is loaded
    # =========================================================================
    observeEvent(data(), {
      req(data())
      df <- data()$merged
      
      col_choices <- colnames(df)
      
      updateSelectInput(session, "color_var", choices = c("", col_choices), 
                        selected = "target")
      updateSelectInput(session, "facet_vars", choices = col_choices, 
                        selected = "target")
    })
    
    # =========================================================================
    # Plot wrappers (dynamic width/height)
    # =========================================================================
    output$plot_tm_wrapper <- renderUI({
      list(input$plot_width, input$plot_height)
      div(
        style = sprintf("width: %sin;", input$plot_width),
        plotOutput(ns("plot_tm"), height = sprintf("%sin", input$plot_height))
      )
    })
    
    output$plot_dtm_wrapper <- renderUI({
      list(input$plot_width, input$plot_height)
      div(
        style = sprintf("width: %sin;", input$plot_width),
        plotOutput(ns("plot_dtm"), height = sprintf("%sin", input$plot_height))
      )
    })
    
    output$plot_raw_wrapper <- renderUI({
      list(input$plot_width, input$plot_height)
      div(
        style = sprintf("width: %sin;", input$plot_width),
        plotOutput(ns("plot_raw"), height = sprintf("%sin", input$plot_height))
      )
    })
    
    output$plot_deri_wrapper <- renderUI({
      list(input$plot_width, input$plot_height)
      div(
        style = sprintf("width: %sin;", input$plot_width),
        plotOutput(ns("plot_deri"), height = sprintf("%sin", input$plot_height))
      )
    })
    
    output$plot_ligand_detail_wrapper <- renderUI({
      list(input$plot_width, input$plot_height)
      div(
        style = sprintf("width: %sin;", input$plot_width),
        plotOutput(ns("plot_ligand_detail"), height = sprintf("%sin", input$plot_height + 1))
      )
    })
    
    # =========================================================================
    # TM Scatter Plot
    # =========================================================================
    output$plot_tm <- renderPlot({
      req(data(), df_with_dtm())
      
      # Explicitly register all plot parameter reactive dependencies
      list(input$color_var, input$facet_vars, input$xlim_min, input$xlim_max, 
           input$facet_ncol, input$plot_width, input$plot_height)
      
      df <- df_with_dtm() %>%
        mutate(tm = ifelse(tm > 100, 100, ifelse(tm < 0, 0, tm)))
      
      build_scatter_plot(df, "ligand", "tm", input$color_var, input$facet_vars,
                          input$xlim_min, input$xlim_max, ncol = input$facet_ncol)
    })
    
    # =========================================================================
    # Delta TM Scatter Plot
    # =========================================================================
    output$plot_dtm <- renderPlot({
      req(data(), df_with_dtm())
      
      list(input$color_var, input$facet_vars, input$xlim_min, input$xlim_max, 
           input$facet_ncol, input$plot_width, input$plot_height)
      
      if (!rv$has_ref_ligand) {
        plot(NULL, xlim = c(0, 1), ylim = c(0, 1), 
             xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n')
        text(0.5, 0.5, "Please set Reference Ligand\nto view Delta TM plots", 
             cex = 1.5, col = "red")
        return(NULL)
      }
      
      df <- df_with_dtm() %>%
        mutate(d_tm = ifelse(d_tm > 10, 10, ifelse(d_tm < -10, -10, d_tm)))
      
      build_scatter_plot(df, "ligand", "d_tm", input$color_var, input$facet_vars,
                          input$xlim_min, input$xlim_max, hline_y = 0, ncol = input$facet_ncol)
    })
    
    # =========================================================================
    # Raw Curves Plot
    # =========================================================================
    output$plot_raw <- renderPlot({
      req(data())
      
      list(input$color_var, input$facet_vars, input$facet_ncol, input$plot_width, input$plot_height)
      
      build_curve_plot(data()$mc_tidy, rv$ref_ligand_val, 
                        input$color_var, input$facet_vars, ncol = input$facet_ncol)
    })
    
    # =========================================================================
    # Derivative Curves Plot
    # =========================================================================
    output$plot_deri <- renderPlot({
      req(data())
      
      list(input$color_var, input$facet_vars, input$facet_ncol, input$plot_width, input$plot_height)
      
      deri <- data()$derivatives
      validate(need(nrow(deri) > 0, "No derivatives could be computed"))
      
      build_curve_plot(deri, rv$ref_ligand_val, 
                        input$color_var, input$facet_vars, 
                        is_derivative = TRUE, x_limits = c(25, 55), ncol = input$facet_ncol)
    })
    
    # =========================================================================
    # Ligand Detail Plot
    # =========================================================================
    output$ligand_selector <- renderUI({
      req(data())
      selectInput(ns("selected_ligand"), "Select Ligand", 
                  choices = data()$ligand_choices, 
                  selected = data()$ligand_choices[1])
    })
    
    output$plot_ligand_detail <- renderPlot({
      req(data(), input$selected_ligand)
      
      list(input$color_var, input$facet_ncol, input$plot_width, input$plot_height)
      
      roi <- input$selected_ligand
      ref_ligand_val <- rv$ref_ligand_val
      has_ref <- !is.null(ref_ligand_val) && ref_ligand_val != "" && 
                 ref_ligand_val %in% data()$mc_tidy$ligand
      
      # Subset raw data
      if (has_ref) {
        df_sub <- data()$mc_tidy %>%
          filter(ligand %in% c(ref_ligand_val, roi))
        deri_sub <- data()$derivatives %>%
          filter(ligand %in% c(ref_ligand_val, roi))
      } else {
        df_sub <- data()$mc_tidy %>%
          filter(ligand == roi)
        deri_sub <- data()$derivatives %>%
          filter(ligand == roi)
      }
      
      multi_plate <- n_distinct(df_sub$plate) > 1
      facet_vars <- if (multi_plate) vars(target, plate) else vars(target)
      
      p1 <- build_curve_plot(df_sub, ref_ligand_val,
                             input$color_var, NULL) +
        facet_wrap(facet_vars, ncol = input$facet_ncol) +
        labs(title = roi) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      suppressMessages(p1 <- p1 + guides(color = guide_legend('')))
      
      p2 <- build_curve_plot(deri_sub, ref_ligand_val,
                             input$color_var, NULL,
                             is_derivative = TRUE, x_limits = c(25, 55)) +
        facet_wrap(facet_vars, ncol = input$facet_ncol) +
        labs(title = roi) +
        theme(legend.position = 'none',
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      
      # Chemical structure (if SMILES available)
      meta <- data()$meta
      if (roi %in% meta$ligand && 'smiles' %in% colnames(meta) && requireNamespace("rcdk", quietly = TRUE)) {
        tryCatch({
          mol <- meta %>%
            filter(ligand == roi) %>%
            pull(smiles) %>% unique() %>%
            rcdk::parse.smiles() %>% .[[1]]
          grob <- grid::rasterGrob(view.image.2d(mol))
          p3 <- ggplot() + annotation_custom(grob) +
            labs(title = roi) + theme_void()
          print(p1 + p2 + p3 + plot_layout(widths = c(1, 1, 0.2)))
        }, error = function(e) {
          print(p1 + p2 + plot_spacer() + plot_layout(widths = c(1, 1, 0.2)))
        })
      } else {
        print(p1 + p2 + plot_spacer() + plot_layout(widths = c(1, 1, 0.2)))
      }
    })
    
    # =========================================================================
    # Data Preview
    # =========================================================================
    output$data_preview <- DT::renderDataTable({
      req(data())
      df_with_dtm()
    }, options = list(pageLength = 20, scrollX = TRUE))
    
    # =========================================================================
    # Download Data
    # =========================================================================
    output$download_data <- downloadHandler(
      filename = function() {
        paste0('DSF_', Sys.Date(), '_merged.xlsx')
      },
      content = function(file) {
        req(data())
        writexl::write_xlsx(list(
          "Merged Data" = df_with_dtm(),
          "Wide Format" = df_wide(),
          "Metadata" = data()$meta,
          "Raw Fluorescence" = data()$mc_tidy,
          "1st Derivative" = data()$derivatives
        ), file)
      }
    )
    
    # =========================================================================
    # Download Plots
    # =========================================================================
    output$download_plots <- downloadHandler(
      filename = function() {
        paste0('DSF_', Sys.Date(), '_plots.zip')
      },
      content = function(file) {
        req(data())
        
        plot_dir <- tempfile(pattern = "dsf_plots_")
        dir.create(plot_dir, recursive = TRUE)
        
        tryCatch({
          df <- df_with_dtm()
          mc_tidy <- data()$mc_tidy
          deri <- data()$derivatives
          meta <- data()$meta
          ref_ligand_val <- rv$ref_ligand_val
          color_col <- input$color_var
          facet_vars <- input$facet_vars
          
          withProgress(message = 'Generating plots...', value = 0, {
            
            incProgress(0.15, detail = "TM plot...")
            p_tm <- df %>%
              mutate(tm = ifelse(tm > 100, 100, ifelse(tm < 0, 0, tm))) %>%
              {build_scatter_plot(., "ligand", "tm", color_col, facet_vars, ncol = input$facet_ncol)}
            ggsave(file.path(plot_dir, 'TM_scatter.pdf'), p_tm, 
                   width = input$plot_width, height = input$plot_height)
            
            incProgress(0.15, detail = "Delta TM plot...")
            if (rv$has_ref_ligand) {
              p_dtm <- df %>%
                mutate(d_tm = ifelse(d_tm > 10, 10, ifelse(d_tm < -10, -10, d_tm))) %>%
                {build_scatter_plot(., "ligand", "d_tm", color_col, facet_vars, 
                                     hline_y = 0, ncol = input$facet_ncol)}
              ggsave(file.path(plot_dir, 'Delta_TM_scatter.pdf'), p_dtm, 
                     width = input$plot_width, height = input$plot_height)
            }
            
            incProgress(0.15, detail = "Raw curves...")
            p_raw <- build_curve_plot(mc_tidy, ref_ligand_val, color_col, facet_vars, ncol = input$facet_ncol)
            ggsave(file.path(plot_dir, 'Raw_curves.pdf'), p_raw, 
                   width = input$plot_width, height = input$plot_height)
            
            incProgress(0.15, detail = "Derivative curves...")
            if (nrow(deri) > 0) {
              p_deri <- build_curve_plot(deri, ref_ligand_val, color_col, facet_vars,
                                          is_derivative = TRUE, x_limits = c(25, 55), ncol = input$facet_ncol)
              ggsave(file.path(plot_dir, 'Derivative_curves.pdf'), p_deri, 
                     width = input$plot_width, height = input$plot_height)
            }
            
            incProgress(0.2, detail = "Ligand details...")
            rois <- unique(meta$ligand)
            for (roi in rois[1:min(20, length(rois))]) {
              if (!is.null(ref_ligand_val) && ref_ligand_val != "" && ref_ligand_val %in% mc_tidy$ligand) {
                df_sub <- mc_tidy %>% filter(ligand %in% c(ref_ligand_val, roi))
                deri_sub <- deri %>% filter(ligand %in% c(ref_ligand_val, roi))
              } else {
                df_sub <- mc_tidy %>% filter(ligand == roi)
                deri_sub <- deri %>% filter(ligand == roi)
              }
              
              multi_plate <- n_distinct(df_sub$plate) > 1
              fvars <- if (multi_plate) vars(target, plate) else vars(target)
              
              p1 <- build_curve_plot(df_sub, ref_ligand_val, color_col, NULL) +
                facet_wrap(fvars, ncol = input$facet_ncol) +
                labs(title = roi) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
              suppressMessages(p1 <- p1 + guides(color = guide_legend('')))
              
              p2 <- build_curve_plot(deri_sub, ref_ligand_val, color_col, NULL,
                                      is_derivative = TRUE, x_limits = c(25, 55)) +
                facet_wrap(fvars, ncol = input$facet_ncol) +
                labs(title = roi) +
                theme(legend.position = 'none',
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
              
              if (roi %in% meta$ligand && 'smiles' %in% colnames(meta) && requireNamespace("rcdk", quietly = TRUE)) {
                tryCatch({
                  mol <- meta %>% filter(ligand == roi) %>% pull(smiles) %>% unique() %>%
                    rcdk::parse.smiles() %>% .[[1]]
                  grob <- grid::rasterGrob(view.image.2d(mol))
                  p3 <- ggplot() + annotation_custom(grob) + labs(title = roi) + theme_void()
                  pdf(file.path(plot_dir, str_glue('ligand_{roi}.pdf')), width = 15, height = 4)
                  print(p1 + p2 + p3 + plot_layout(widths = c(1, 1, 0.2)))
                  dev.off()
                }, error = function(e) {
                  pdf(file.path(plot_dir, str_glue('ligand_{roi}.pdf')), width = 15, height = 4)
                  print(p1 + p2 + plot_spacer() + plot_layout(widths = c(1, 1, 0.2)))
                  dev.off()
                })
              } else {
                pdf(file.path(plot_dir, str_glue('ligand_{roi}.pdf')), width = 15, height = 4)
                print(p1 + p2 + plot_spacer() + plot_layout(widths = c(1, 1, 0.2)))
                dev.off()
              }
            }
            
            incProgress(0.1, detail = "Creating zip...")
            pdf_files <- list.files(plot_dir, pattern = "\\.pdf$", full.names = FALSE)
            
            if (length(pdf_files) == 0) {
              showNotification("No plots were generated", type = "error", duration = 10)
              return(NULL)
            }
            
            old_wd <- getwd()
            setwd(plot_dir)
            on.exit(setwd(old_wd), add = TRUE)
            zip::zip(zipfile = file, files = pdf_files)
          })
          
        }, error = function(e) {
          showNotification(paste("Error creating zip:", e$message), 
                           type = "error", duration = 10)
        })
      },
      contentType = "application/zip"
    )
  })
}
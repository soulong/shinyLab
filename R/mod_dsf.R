# shinyLab - DSF Module
# Differential Scanning Fluorimetry analysis module

library(patchwork)
library(mgcv)
library(gratia)
library(rcdk)
library(rio)
library(DT)
library(zip)
library(janitor)
options(rio.import.class='tbl')

# =============================================================================
# UI Function
# =============================================================================

dsfUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(width = 3,
      box(title = "File Upload", width = 12, status = "primary", solidHeader = TRUE,
        fileInput(ns("file_analysis"), "Analysis Results (.txt)", 
                  accept = c(".txt", ".tsv")),
        fileInput(ns("file_rawdata"), "Raw Data Files (multiple)", 
                  accept = c(".txt", ".tsv"), multiple = TRUE),
        checkboxInput(ns("use_custom_metadata"), "Use custom metadata file", value = FALSE),
        conditionalPanel(
          condition = "input.use_custom_metadata == true",
          ns = ns,
          fileInput(ns("file_metadata"), "Plate Info (.xlsx)", accept = ".xlsx")
        ),
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
            plotOutput(ns("plot_tm"), height = "600px")
          ),
          tabPanel("Delta TM Scatter",
            plotOutput(ns("plot_dtm"), height = "600px")
          ),
          tabPanel("Raw Curves",
            plotOutput(ns("plot_raw"), height = "600px")
          ),
          tabPanel("Derivative Curves",
            plotOutput(ns("plot_deri"), height = "600px")
          ),
          tabPanel("Ligand Details",
            uiOutput(ns("ligand_selector")),
            plotOutput(ns("plot_ligand_detail"), height = "500px")
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
      info = NULL,
      has_ref_ligand = FALSE,
      ref_ligand_val = "",
      error_message = NULL
    )
    
    # =========================================================================
    # Process data on submit
    # =========================================================================
    data <- eventReactive(input$submit, {
      req(input$file_analysis)
      req(input$file_rawdata)
      
      tryCatch({
        withProgress(message = 'Processing data...', value = 0, {
          
          incProgress(0.2, detail = "Reading analysis results...")
          
          # Read analysis results
          raw <- import(input$file_analysis$datapath) %>%
            janitor::clean_names()
          
          incProgress(0.2, detail = "Reading metadata...")
          
          # Conditionally read metadata
          if (input$use_custom_metadata) {
            req(input$file_metadata)
            info <- read_metadata(input$file_metadata$datapath) %>%
              janitor::remove_constant()
            
            # Merge data
            merged <- raw %>%
              janitor::remove_empty('cols') %>%
              mutate(well = transform_well_style(well)) %>%
              select(plate = experiment_file_name, well, tm = tm_d) %>%
              mutate(plate = map_chr(plate, \(x) str_replace_all(x, '.*_Admin_', '') %>%
                                       str_replace_all('.eds', ''))) %>%
              right_join(info, .)
          } else {
            # Validate required columns
            if (!all(c("target", "ligand") %in% colnames(raw))) {
              stop("AnalysisResults must contain 'target' and 'ligand' columns when not using custom metadata")
            }
            
            # Extract relevant columns from raw data
            merged <- raw %>%
              janitor::remove_empty('cols') %>%
              mutate(well = transform_well_style(well)) %>%
              select(plate = experiment_file_name,
                     well,
                     target,
                     ligand,
                     tm = tm_d,
                     any_of(c("conc", "smiles"))) %>%
              mutate(plate = map_chr(plate, \(x) str_replace_all(x, '.*_Admin_', '') %>%
                                       str_replace_all('.eds', '')))
            
            # Create info dataframe
            info <- merged %>%
              select(any_of(c("plate", "well", "target", "ligand", "conc", "smiles"))) %>%
              distinct()
          }
          
          merged <- merged %>%
            filter(!is.na(ligand), !is.na(target)) %>%
            relocate(plate, well, target, .before = 1) %>%
            relocate(tm, .after = last_col())
          
          # Add conc as numeric if exists
          if ("conc" %in% colnames(merged)) {
            merged <- merged %>% mutate(conc = as.numeric(conc))
          }
          
          incProgress(0.2, detail = "Processing raw data...")
          
          # Read raw data files
          mc1 <- input$file_rawdata$datapath %>%
            set_names(nm = map_chr(input$file_rawdata$name, \(x) str_replace_all(x, '.*_Admin_', '') %>%
                                                          str_replace_all('.eds.txt', ''))) %>%
            map(\(x) import(x)) %>%
            list_rbind(names_to = 'plate') %>%
            janitor::clean_names() %>%
            select(!c(well)) %>%
            rename(well = well_position) %>%
            mutate(well = transform_well_style(well))
          
          incProgress(0.2, detail = "Merging raw data...")
          
          # Tidy raw data
          mc_tidy <- mc1 %>%
            left_join(info) %>%
            mutate(temperature = as.numeric(temperature),
                   fluorescence = as.numeric(fluorescence)) %>%
            reframe(fluorescence = mean(fluorescence), 
                    .by = c(plate, well, target, ligand, temperature)) %>%
            filter(!is.na(ligand), !is.na(target))
          
          # Create wide format (without d_tm, will be calculated later)
          df_wide_base <- merged %>%
            select(-any_of("tm"))
          
          incProgress(0.2, detail = "Done!")
          
          rv$info <- info
          rv$mc_tidy <- mc_tidy
          
          return(list(
            df_base = merged,  # Base data without d_tm
            df_wide_base = df_wide_base,
            mc_tidy = mc_tidy, 
            info = info,
            ligand_choices = sort(unique(info$ligand)),
            target_choices = sort(unique(info$target))
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
      
      ligand_choices <- data()$ligand_choices
      target_choices <- data()$target_choices
      
      # Auto-select ref_ligand
      ref_ligand_default <- find_best_match(ligand_choices, c("blank", "dmso", "pbs"))
      
      # Auto-select ref_target
      ref_target_default <- find_best_match(target_choices, c("blank", "empty", "noprotein", "no_protein"))
      
      tagList(
        selectInput(ns("ref_ligand_input"), "Reference Ligand",
                    choices = ligand_choices,
                    selected = ref_ligand_default),
        selectInput(ns("ref_target_input"), "Reference Target",
                    choices = target_choices,
                    selected = ref_target_default)
      )
    })
    
    # =========================================================================
    # Calculate d_tm reactively when ref changes
    # =========================================================================
    df_with_dtm <- reactive({
      req(data())
      
      df_base <- data()$df_base
      ref_ligand_val <- input$ref_ligand_input
      
      # Validate
      if (is.null(ref_ligand_val) || ref_ligand_val == "") {
        rv$has_ref_ligand <- FALSE
        rv$ref_ligand_val <- ""
        return(df_base %>% mutate(d_tm = NA_real_))
      }
      
      if (!ref_ligand_val %in% df_base$ligand) {
        showNotification(paste("Reference ligand", ref_ligand_val, "not found in data"), 
                         type = "warning", duration = 5)
        rv$has_ref_ligand <- FALSE
        rv$ref_ligand_val <- ""
        return(df_base %>% mutate(d_tm = NA_real_))
      }
      
      # Calculate Delta TM
      ref <- df_base %>%
        filter(ligand == ref_ligand_val) %>%
        reframe(tm_ref = median(tm), .by = c(target))
      
      df <- df_base %>%
        left_join(ref, by = "target") %>%
        mutate(d_tm = tm - tm_ref) %>%
        select(!tm_ref)
      
      rv$has_ref_ligand <- TRUE
      rv$ref_ligand_val <- ref_ligand_val
      
      return(df)
    })
    
    # Create wide format with d_tm
    df_wide <- reactive({
      req(df_with_dtm())
      
      df <- df_with_dtm()
      info <- data()$info
      
      df_wide <- df %>%
        pivot_wider(id_cols = any_of(setdiff(colnames(info), 'well')),
                    names_from = target, values_from = d_tm,
                    names_prefix = 'dTM_',
                    values_fn = median) %>%
        mutate(across(where(is.numeric), \(x) round(x, 3)))
      
      return(df_wide)
    })
    
    # =========================================================================
    # Update select inputs after data is loaded
    # =========================================================================
    observeEvent(data(), {
      req(data())
      df <- data()$df_base
      
      col_choices <- colnames(df)
      
      updateSelectInput(session, "color_var", choices = c("", col_choices), 
                        selected = "target")
      updateSelectInput(session, "facet_vars", choices = col_choices, 
                        selected = "target")
    })
    
    # =========================================================================
    # TM Scatter Plot
    # =========================================================================
    output$plot_tm <- renderPlot({
      req(data())
      req(df_with_dtm())
      
      withProgress(message = 'Rendering plot...', value = 0, {
        incProgress(0.3, detail = "Preparing data...")
        
        df <- df_with_dtm()
        color_col <- input$color_var
        use_color <- color_col != "" && color_col %in% names(df)
        
        p <- df %>%
          mutate(tm = ifelse(tm > 100, 100, ifelse(tm < 0, 0, tm)))
        
        if (use_color) {
          p <- p %>%
            ggplot(aes(x = ligand, y = tm, color = .data[[color_col]])) +
            geom_point(show.legend = TRUE)
        } else {
          p <- p %>%
            ggplot(aes(x = ligand, y = tm)) +
            geom_point(show.legend = FALSE)
        }
        
        if (length(input$facet_vars) > 0 && !all(input$facet_vars == "")) {
          p <- p + facet_wrap(as.formula(paste("~", paste(input$facet_vars, collapse = " + "))), 
                              scales = "fixed")
        }
        
        p <- p + theme_bw() + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        
        if (!is.null(input$xlim_min) && !is.null(input$xlim_max)) {
          p <- p + coord_cartesian(xlim = c(input$xlim_min, input$xlim_max))
        }
        
        incProgress(0.7, detail = "Finalizing...")
        
        print(p)
      })
    })
    
    # =========================================================================
    # Delta TM Scatter Plot
    # =========================================================================
    output$plot_dtm <- renderPlot({
      req(data())
      req(df_with_dtm())
      
      if (!rv$has_ref_ligand) {
        plot(NULL, xlim = c(0, 1), ylim = c(0, 1), 
             xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n')
        text(0.5, 0.5, "Please set Reference Ligand\nto view Delta TM plots", 
             cex = 1.5, col = "red")
        return(NULL)
      }
      
      withProgress(message = 'Rendering plot...', value = 0, {
        incProgress(0.3, detail = "Preparing data...")
        
        df <- df_with_dtm()
        color_col <- input$color_var
        use_color <- color_col != "" && color_col %in% names(df)
        
        p <- df %>%
          mutate(d_tm = ifelse(d_tm > 10, 10, ifelse(d_tm < -10, -10, d_tm)))
        
        if (use_color) {
          p <- p %>%
            ggplot(aes(x = ligand, y = d_tm, color = .data[[color_col]])) +
            geom_point(show.legend = TRUE)
        } else {
          p <- p %>%
            ggplot(aes(x = ligand, y = d_tm)) +
            geom_point(show.legend = FALSE)
        }
        
        p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
        
        if (length(input$facet_vars) > 0 && !all(input$facet_vars == "")) {
          p <- p + facet_wrap(as.formula(paste("~", paste(input$facet_vars, collapse = " + "))), 
                              scales = "fixed")
        }
        
        p <- p + theme_bw() + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        
        if (!is.null(input$xlim_min) && !is.null(input$xlim_max)) {
          p <- p + coord_cartesian(xlim = c(input$xlim_min, input$xlim_max))
        }
        
        incProgress(0.7, detail = "Finalizing...")
        
        print(p)
      })
    })
    
    # =========================================================================
    # Raw Curves Plot
    # =========================================================================
    output$plot_raw <- renderPlot({
      req(data())
      
      withProgress(message = 'Rendering plot...', value = 0, {
        incProgress(0.3, detail = "Preparing data...")
        
        mc_tidy <- data()$mc_tidy
        
        ref_ligand_val <- rv$ref_ligand_val
        has_ref <- !is.null(ref_ligand_val) && ref_ligand_val != "" && ref_ligand_val %in% mc_tidy$ligand
        
        color_by <- input$color_var
        is_valid_color <- !is.null(color_by) && color_by != "" && 
                          color_by %in% names(mc_tidy) && color_by != "ligand"
        
        is_color_numeric <- FALSE
        if (is_valid_color) {
          is_color_numeric <- is.numeric(mc_tidy[[color_by]])
        }
        
        p <- ggplot()
        
        if (has_ref) {
          p <- p +
            geom_path(
              data = filter(mc_tidy, ligand == ref_ligand_val),
              aes(temperature, fluorescence, group = well),
              color = 'grey50',
              show.legend = TRUE
            )
        }
        
        if (is_valid_color) {
          filter_data <- if (has_ref) filter(mc_tidy, ligand != ref_ligand_val) else mc_tidy
          p <- p +
            geom_path(
              data = filter_data,
              aes(temperature, fluorescence, group = well,
                  color = .data[[color_by]]),
              show.legend = TRUE
            )
          if (is_color_numeric) {
            p <- p + scale_color_viridis_c(option = 'C', na.value = 'darkred')
          } else {
            p <- p + scale_color_viridis_d(option = 'C', na.value = 'darkred')
          }
        } else {
          filter_data <- if (has_ref) filter(mc_tidy, ligand != ref_ligand_val) else mc_tidy
          p <- p +
            geom_path(
              data = filter_data,
              aes(temperature, fluorescence, group = well),
              color = 'darkred',
              show.legend = TRUE
            )
        }
        
        if (length(input$facet_vars) > 0 && !all(input$facet_vars == "")) {
          p <- p + facet_wrap(as.formula(paste("~", paste(input$facet_vars, collapse = " + "))), 
                              ncol = 6)
        }
        
        p <- p + theme_bw() + 
          theme(legend.position = 'top')
        
        incProgress(0.7, detail = "Finalizing...")
        
        print(p)
      })
    })
    
    # =========================================================================
    # Derivative Curves Plot
    # =========================================================================
    output$plot_deri <- renderPlot({
      req(data())
      
      withProgress(message = 'Rendering plot...', value = 0, {
        incProgress(0.2, detail = "Preparing data...")
        
        mc_tidy <- data()$mc_tidy
        
        ref_ligand_val <- rv$ref_ligand_val
        has_ref <- !is.null(ref_ligand_val) && ref_ligand_val != "" && ref_ligand_val %in% mc_tidy$ligand
        
        color_by <- input$color_var
        is_valid_color <- !is.null(color_by) && color_by != "" && 
                          color_by %in% names(mc_tidy) && color_by != "ligand"
        
        is_color_numeric <- FALSE
        if (is_valid_color) {
          is_color_numeric <- is.numeric(mc_tidy[[color_by]])
        }
        
        incProgress(0.3, detail = "Computing derivatives...")
        
        deri <- mc_tidy %>%
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
        
        incProgress(0.3, detail = "Plotting...")
        
        p <- ggplot()
        
        if (has_ref) {
          p <- p +
            geom_line(
              data = filter(deri, ligand == ref_ligand_val),
              aes(temperature, derivative, group = well),
              color = 'grey50',
              show.legend = TRUE
            )
        }
        
        if (is_valid_color) {
          filter_data <- if (has_ref) filter(deri, ligand != ref_ligand_val) else deri
          p <- p +
            geom_line(
              data = filter_data,
              aes(temperature, derivative, group = well,
                  color = .data[[color_by]]),
              show.legend = TRUE
            )
          if (is_color_numeric) {
            p <- p + scale_color_viridis_c(option = 'C', na.value = 'darkred')
          } else {
            p <- p + scale_color_viridis_d(option = 'C', na.value = 'darkred')
          }
        } else {
          filter_data <- if (has_ref) filter(deri, ligand != ref_ligand_val) else deri
          p <- p +
            geom_line(
              data = filter_data,
              aes(temperature, derivative, group = well),
              color = 'darkred',
              show.legend = TRUE
            )
        }
        
        if (length(input$facet_vars) > 0 && !all(input$facet_vars == "")) {
          p <- p + facet_wrap(as.formula(paste("~", paste(input$facet_vars, collapse = " + "))), 
                              ncol = 6)
        }
        
        p <- p + coord_cartesian(xlim = c(25, 55)) +
          theme_bw() +
          theme(legend.position = 'top')
        
        incProgress(0.2, detail = "Finalizing...")
        
        print(p)
      })
    })
    
    # =========================================================================
    # Ligand Detail Plot
    # =========================================================================
    output$ligand_selector <- renderUI({
      req(data())
      info <- data()$info
      selectInput(ns("selected_ligand"), "Select Ligand", 
                  choices = unique(info$ligand), 
                  selected = unique(info$ligand)[1])
    })
    
    output$plot_ligand_detail <- renderPlot({
      req(data())
      req(input$selected_ligand)
      
      withProgress(message = 'Rendering plot...', value = 0, {
        incProgress(0.2, detail = "Preparing data...")
        
        mc_tidy <- data()$mc_tidy
        info <- data()$info
        
        roi <- input$selected_ligand
        ref_ligand_val <- rv$ref_ligand_val
        has_ref <- !is.null(ref_ligand_val) && ref_ligand_val != "" && ref_ligand_val %in% mc_tidy$ligand
        
        if (has_ref) {
          df_sub <- mc_tidy %>%
            filter(ligand %in% c(ref_ligand_val, roi))
        } else {
          df_sub <- mc_tidy %>%
            filter(ligand == roi)
        }
        
        multi_plate <- n_distinct(df_sub$plate) > 1
        facet_vars <- if (multi_plate) vars(target, plate) else vars(target)
        
        color_by_used <- input$color_var
        is_valid_color <- !is.null(color_by_used) && color_by_used != "" && 
                          color_by_used %in% names(df_sub) && color_by_used != "ligand"
        
        is_color_numeric <- FALSE
        if (is_valid_color) {
          is_color_numeric <- is.numeric(df_sub[[color_by_used]])
        }
        
        incProgress(0.3, detail = "Plotting raw curves...")
        
        # Raw curve
        p1 <- ggplot()
        
        if (has_ref) {
          p1 <- p1 +
            geom_path(
              data = filter(df_sub, ligand == ref_ligand_val),
              aes(temperature, fluorescence, group = well),
              color = 'grey50',
              show.legend = TRUE
            )
        }
        
        if (is_valid_color) {
          filter_data <- if (has_ref) filter(df_sub, ligand != ref_ligand_val) else df_sub
          p1 <- p1 +
            geom_path(
              data = filter_data,
              aes(temperature, fluorescence, group = well,
                  color = .data[[color_by_used]]),
              show.legend = TRUE
            ) +
            (if (is_color_numeric) scale_color_viridis_c(option = 'C', na.value = 'darkred')
             else scale_color_viridis_d(option = 'C', na.value = 'darkred'))
        } else {
          filter_data <- if (has_ref) filter(df_sub, ligand != ref_ligand_val) else df_sub
          p1 <- p1 +
            geom_path(
              data = filter_data,
              aes(temperature, fluorescence, group = well),
              color = 'darkred',
              show.legend = TRUE
            )
        }
        
        p1 <- p1 +
          guides(color = guide_legend('')) +
          facet_wrap(facet_vars, ncol = 6) +
          labs(title = roi) +
          theme_bw() +
          theme(
            legend.position = 'top',
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
          )
        
        incProgress(0.3, detail = "Plotting derivative curves...")
        
        # Derivative curve
        deri <- df_sub %>%
          nest(.by = c(plate, target, ligand, well)) %>%
          mutate(pred = purrr::map(data, \(d) {
            tryCatch({
              mod <- mgcv::gam(fluorescence ~ s(temperature), data = d)
              res <- gratia::derivatives(mod, n = 100, order = 1)
              tibble(temperature = res$temperature,
                     derivative = res$.derivative)
            }, error = function(e) NULL)
          })) %>%
          unnest(pred)
        
        p2 <- ggplot()
        
        if (has_ref) {
          p2 <- p2 +
            geom_line(
              data = filter(deri, ligand == ref_ligand_val),
              aes(temperature, derivative, group = well),
              color = 'grey50',
              show.legend = TRUE
            )
        }
        
        if (is_valid_color) {
          filter_data <- if (has_ref) filter(deri, ligand != ref_ligand_val) else deri
          p2 <- p2 +
            geom_line(
              data = filter_data,
              aes(temperature, derivative, group = well,
                  color = .data[[color_by_used]]),
              show.legend = TRUE
            ) +
            (if (is_color_numeric) scale_color_viridis_c(option = 'C', na.value = 'darkred')
             else scale_color_viridis_d(option = 'C', na.value = 'darkred'))
        } else {
          filter_data <- if (has_ref) filter(deri, ligand != ref_ligand_val) else deri
          p2 <- p2 +
            geom_line(
              data = filter_data,
              aes(temperature, derivative, group = well),
              color = 'darkred',
              show.legend = TRUE
            )
        }
        
        p2 <- p2 +
          facet_wrap(facet_vars, ncol = 6) +
          coord_cartesian(xlim = c(25, 55)) +
          labs(title = roi) +
          theme_bw() +
          theme(
            legend.position = 'none',
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
          )
        
        incProgress(0.2, detail = "Finalizing...")
        
        # Chemical structure
        if (roi %in% info$ligand && 'smiles' %in% colnames(info)) {
          try_res <- try({
            mol <- info %>%
              filter(ligand == roi) %>%
              pull(smiles) %>%
              unique() %>%
              rcdk::parse.smiles() %>%
              .[[1]]
            grob <- grid::rasterGrob(view.image.2d(mol))
            p3 <- ggplot() +
              annotation_custom(grob) +
              labs(title = roi) +
              theme_void()
            print(p1 + p2 + p3 + plot_layout(widths = c(1, 1, 0.2)))
          }, silent = TRUE)
          
          if (inherits(try_res, "try-error")) {
            print(p1 + p2 + plot_spacer() + plot_layout(widths = c(1, 1, 0.2)))
          }
        } else {
          print(p1 + p2 + plot_spacer() + plot_layout(widths = c(1, 1, 0.2)))
        }
      })
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
          "Metadata" = data()$info
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
        
        # Create temporary directory
        plot_dir <- tempfile(pattern = "dsf_plots_")
        dir.create(plot_dir, recursive = TRUE)
        
        tryCatch({
          df <- df_with_dtm()
          mc_tidy <- data()$mc_tidy
          info <- data()$info
          ref_ligand_val <- rv$ref_ligand_val
          has_ref <- !is.null(ref_ligand_val) && ref_ligand_val != "" && ref_ligand_val %in% mc_tidy$ligand
          
          withProgress(message = 'Generating plots...', value = 0, {
            
            incProgress(0.15, detail = "TM plot...")
            # TM plot
            color_col <- input$color_var
            use_color <- color_col != "" && color_col %in% names(df)
            
            p_tm <- df %>%
              mutate(tm = ifelse(tm > 100, 100, ifelse(tm < 0, 0, tm)))
            
            if (use_color) {
              p_tm <- p_tm %>%
                ggplot(aes(x = ligand, y = tm, color = .data[[color_col]])) +
                geom_point(show.legend = TRUE)
            } else {
              p_tm <- p_tm %>%
                ggplot(aes(x = ligand, y = tm)) +
                geom_point(show.legend = FALSE)
            }
            
            if (length(input$facet_vars) > 0 && !all(input$facet_vars == "")) {
              p_tm <- p_tm + facet_wrap(as.formula(paste("~", paste(input$facet_vars, collapse = " + "))), 
                                         scales = "fixed")
            }
            
            p_tm <- p_tm + theme_bw() + 
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
            
            ggsave(file.path(plot_dir, 'TM_scatter.pdf'), p_tm, 
                   width = input$plot_width, height = input$plot_height)
            
            incProgress(0.15, detail = "Delta TM plot...")
            # Delta TM plot
            if (rv$has_ref_ligand) {
              p_dtm <- df %>%
                mutate(d_tm = ifelse(d_tm > 10, 10, ifelse(d_tm < -10, -10, d_tm)))
              
              if (use_color) {
                p_dtm <- p_dtm %>%
                  ggplot(aes(x = ligand, y = d_tm, color = .data[[color_col]])) +
                  geom_point(show.legend = TRUE)
              } else {
                p_dtm <- p_dtm %>%
                  ggplot(aes(x = ligand, y = d_tm)) +
                  geom_point(show.legend = FALSE)
              }
              
              p_dtm <- p_dtm + geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
              
              if (length(input$facet_vars) > 0 && !all(input$facet_vars == "")) {
                p_dtm <- p_dtm + facet_wrap(as.formula(paste("~", paste(input$facet_vars, collapse = " + "))), 
                                            scales = "fixed")
              }
              
              p_dtm <- p_dtm + theme_bw() + 
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
              
              ggsave(file.path(plot_dir, 'Delta_TM_scatter.pdf'), p_dtm, 
                     width = input$plot_width, height = input$plot_height)
            }
            
            incProgress(0.15, detail = "Raw curves...")
            # Raw curves
            is_valid_color_raw <- !is.null(color_col) && color_col != "" && 
                                  color_col %in% names(mc_tidy) && color_col != "ligand"
            is_color_numeric_raw <- if(is_valid_color_raw) is.numeric(mc_tidy[[color_col]]) else FALSE
            
            p_raw <- ggplot()
            
            if (has_ref) {
              p_raw <- p_raw +
                geom_path(data = filter(mc_tidy, ligand == ref_ligand_val),
                          aes(temperature, fluorescence, group = well), color = 'grey50')
            }
            
            if (is_valid_color_raw) {
              filter_data <- if (has_ref) filter(mc_tidy, ligand != ref_ligand_val) else mc_tidy
              p_raw <- p_raw + geom_path(data = filter_data,
                                         aes(temperature, fluorescence, group = well, color = .data[[color_col]])) +
                (if(is_color_numeric_raw) scale_color_viridis_c(option = 'C', na.value = 'darkred')
                 else scale_color_viridis_d(option = 'C', na.value = 'darkred'))
            } else {
              filter_data <- if (has_ref) filter(mc_tidy, ligand != ref_ligand_val) else mc_tidy
              p_raw <- p_raw + geom_path(data = filter_data,
                                         aes(temperature, fluorescence, group = well), color = 'darkred')
            }
            
            if (length(input$facet_vars) > 0 && !all(input$facet_vars == "")) {
              p_raw <- p_raw + facet_wrap(as.formula(paste("~", paste(input$facet_vars, collapse = " + "))), 
                                          ncol = 6)
            }
            
            p_raw <- p_raw + theme_bw() + theme(legend.position = 'top')
            
            ggsave(file.path(plot_dir, 'Raw_curves.pdf'), p_raw, 
                   width = input$plot_width, height = input$plot_height)
            
            incProgress(0.15, detail = "Derivative curves...")
            # Derivative curves
            deri <- mc_tidy %>%
              nest(.by = c(plate, target, ligand, well)) %>%
              mutate(pred = purrr::map(data, \(d) {
                tryCatch({
                  mod <- mgcv::gam(fluorescence ~ s(temperature), data = d)
                  res <- gratia::derivatives(mod, n = 100, order = 1)
                  tibble(temperature = res$temperature, derivative = res$.derivative)
                }, error = function(e) NULL)
              })) %>%
              unnest(pred) %>%
              filter(!is.na(derivative))
            
            is_valid_color_deri <- !is.null(color_col) && color_col != "" && 
                                   color_col %in% names(deri) && color_col != "ligand"
            is_color_numeric_deri <- if(is_valid_color_deri) is.numeric(deri[[color_col]]) else FALSE
            
            p_deri <- ggplot()
            
            if (has_ref) {
              p_deri <- p_deri +
                geom_line(data = filter(deri, ligand == ref_ligand_val),
                          aes(temperature, derivative, group = well), color = 'grey50')
            }
            
            if (is_valid_color_deri) {
              filter_data <- if (has_ref) filter(deri, ligand != ref_ligand_val) else deri
              p_deri <- p_deri + geom_line(data = filter_data,
                                           aes(temperature, derivative, group = well, color = .data[[color_col]])) +
                (if(is_color_numeric_deri) scale_color_viridis_c(option = 'C', na.value = 'darkred')
                 else scale_color_viridis_d(option = 'C', na.value = 'darkred'))
            } else {
              filter_data <- if (has_ref) filter(deri, ligand != ref_ligand_val) else deri
              p_deri <- p_deri + geom_line(data = filter_data,
                                           aes(temperature, derivative, group = well), color = 'darkred')
            }
            
            if (length(input$facet_vars) > 0 && !all(input$facet_vars == "")) {
              p_deri <- p_deri + facet_wrap(as.formula(paste("~", paste(input$facet_vars, collapse = " + "))), 
                                            ncol = 6)
            }
            
            p_deri <- p_deri + coord_cartesian(xlim = c(25, 55)) +
              theme_bw() + theme(legend.position = 'top')
            
            ggsave(file.path(plot_dir, 'Derivative_curves.pdf'), p_deri, 
                   width = input$plot_width, height = input$plot_height)
            
            incProgress(0.2, detail = "Ligand details...")
            # Individual ligand plots
            rois <- unique(info$ligand)
            
            for (roi in rois[1:min(20, length(rois))]) {
              if (has_ref) {
                df_sub <- mc_tidy %>%
                  filter(ligand %in% c(ref_ligand_val, roi))
              } else {
                df_sub <- mc_tidy %>%
                  filter(ligand == roi)
              }
              
              multi_plate <- n_distinct(df_sub$plate) > 1
              facet_vars <- if (multi_plate) vars(target, plate) else vars(target)
              
              color_by_used <- input$color_var
              is_valid_color <- !is.null(color_by_used) && color_by_used != "" && 
                                color_by_used %in% names(df_sub) && color_by_used != "ligand"
              
              is_color_numeric <- FALSE
              if (is_valid_color) {
                is_color_numeric <- is.numeric(df_sub[[color_by_used]])
              }
              
              p1 <- ggplot()
              
              if (has_ref) {
                p1 <- p1 +
                  geom_path(
                    data = filter(df_sub, ligand == ref_ligand_val),
                    aes(temperature, fluorescence, group = well),
                    color = 'grey50'
                  )
              }
              
              if (is_valid_color) {
                filter_data <- if (has_ref) filter(df_sub, ligand != ref_ligand_val) else df_sub
                p1 <- p1 +
                  geom_path(
                    data = filter_data,
                    aes(temperature, fluorescence, group = well,
                        color = .data[[color_by_used]])
                  ) +
                  (if (is_color_numeric) scale_color_viridis_c(option = 'C', na.value = 'darkred')
                   else scale_color_viridis_d(option = 'C', na.value = 'darkred'))
              } else {
                filter_data <- if (has_ref) filter(df_sub, ligand != ref_ligand_val) else df_sub
                p1 <- p1 +
                  geom_path(
                    data = filter_data,
                    aes(temperature, fluorescence, group = well),
                    color = 'darkred'
                  )
              }
              
              p1 <- p1 +
                guides(color = guide_legend('')) +
                facet_wrap(facet_vars, ncol = 6) +
                labs(title = roi) +
                theme_bw() +
                theme(legend.position = 'top',
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
              
              deri <- df_sub %>%
                nest(.by = c(plate, target, ligand, well)) %>%
                mutate(pred = purrr::map(data, \(d) {
                  tryCatch({
                    mod <- mgcv::gam(fluorescence ~ s(temperature), data = d)
                    res <- gratia::derivatives(mod, n = 100, order = 1)
                    tibble(temperature = res$temperature,
                           derivative = res$.derivative)
                  }, error = function(e) NULL)
                })) %>%
                unnest(pred)
              
              p2 <- ggplot()
              
              if (has_ref) {
                p2 <- p2 +
                  geom_line(
                    data = filter(deri, ligand == ref_ligand_val),
                    aes(temperature, derivative, group = well),
                    color = 'grey50'
                  )
              }
              
              if (is_valid_color) {
                filter_data <- if (has_ref) filter(deri, ligand != ref_ligand_val) else deri
                p2 <- p2 +
                  geom_line(
                    data = filter_data,
                    aes(temperature, derivative, group = well,
                        color = .data[[color_by_used]])
                  ) +
                  (if (is_color_numeric) scale_color_viridis_c(option = 'C', na.value = 'darkred')
                   else scale_color_viridis_d(option = 'C', na.value = 'darkred'))
              } else {
                filter_data <- if (has_ref) filter(deri, ligand != ref_ligand_val) else deri
                p2 <- p2 +
                  geom_line(
                    data = filter_data,
                    aes(temperature, derivative, group = well),
                    color = 'darkred'
                  )
              }
              
              p2 <- p2 +
                facet_wrap(facet_vars, ncol = 6) +
                coord_cartesian(xlim = c(25, 55)) +
                labs(title = roi) +
                theme_bw() +
                theme(legend.position = 'none',
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
              
              if (roi %in% info$ligand && 'smiles' %in% colnames(info)) {
                try_res <- try({
                  mol <- info %>%
                    filter(ligand == roi) %>%
                    pull(smiles) %>%
                    unique() %>%
                    rcdk::parse.smiles() %>%
                    .[[1]]
                  grob <- grid::rasterGrob(view.image.2d(mol))
                  p3 <- ggplot() +
                    annotation_custom(grob) +
                    labs(title = roi) +
                    theme_void()
                  pdf(file.path(plot_dir, str_glue('ligand_{roi}.pdf')), width = 15, height = 4)
                  print(p1 + p2 + p3 + plot_layout(widths = c(1, 1, 0.2)))
                  dev.off()
                }, silent = TRUE)
                
                if (inherits(try_res, "try-error")) {
                  pdf(file.path(plot_dir, str_glue('ligand_{roi}.pdf')), width = 15, height = 4)
                  print(p1 + p2 + plot_spacer() + plot_layout(widths = c(1, 1, 0.2)))
                  dev.off()
                }
              } else {
                pdf(file.path(plot_dir, str_glue('ligand_{roi}.pdf')), width = 15, height = 4)
                print(p1 + p2 + plot_spacer() + plot_layout(widths = c(1, 1, 0.2)))
                dev.off()
              }
            }
            
            incProgress(0.1, detail = "Creating zip...")
            
            # Get PDF files
            pdf_files <- list.files(plot_dir, pattern = "\\.pdf$", full.names = FALSE)
            
            if (length(pdf_files) == 0) {
              showNotification("No plots were generated", type = "error", duration = 10)
              return(NULL)
            }
            
            # Create zip
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
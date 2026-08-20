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
  if (nrow(mc_tidy) == 0) {
    return(tibble(plate = character(), target = character(),
                  ligand = character(), well = character(),
                  temperature = numeric(), derivative = numeric()))
  }
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
    unnest(pred, keep_empty = TRUE) %>%
    dplyr::select(!any_of(c("data", "pred"))) %>%
    { if (!"derivative" %in% colnames(.)) mutate(., derivative = NA_real_) else . } %>%
    filter(!is.na(derivative))
}


# -----------------------------------------------------------------------------
# Analysis Results Reader
# -----------------------------------------------------------------------------
# QuantStudio DSF exports append a second "Hits" summary table below the main
# per-well results table, which causes fread() to abort when the trailing
# table has fewer columns. Read only the first (main) table block.
read_analysis_results <- function(file) {
  lines <- readLines(file, warn = FALSE)
  hdr <- grep("^#[\t ]", lines)
  if (length(hdr) < 1) return(NULL)
  first <- hdr[1]
  last <- if (length(hdr) > 1) hdr[2] - 1 else length(lines)
  data.table::fread(text = paste(lines[first:last], collapse = "\n"),
                    data.table = FALSE) |> tibble::as_tibble()
}

# RawData exports append a second "Derivative" table below the fluorescence
# data. Read only the fluorescence table.
read_fluorescence_table <- function(file) {
  lines <- readLines(file, warn = FALSE)
  hdr <- grep("^Well\tWell Position\tReading\tTemperature\tFluorescence", lines)
  if (length(hdr) < 1) return(NULL)
  next_hdr <- grep("^Well\tWell Position\tReading\tTemperature\tDerivative", lines)
  last <- if (length(next_hdr) > 0) next_hdr[1] - 1 else length(lines)
  data.table::fread(text = paste(lines[hdr[1]:last], collapse = "\n"),
                    data.table = FALSE) |> tibble::as_tibble()
}

# =============================================================================
# UI Function
# =============================================================================

dsfUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(width = 3,
      box(title = "File Upload", width = NULL, status = "primary", solidHeader = TRUE,
        fileInput(ns("file_analysis"), "Analysis Data Files (multiple)", 
                  accept = c(".txt", ".tsv"), multiple = TRUE),
        helpText("Upload the QuantStudio DSF export files: the '...AnalysisResults.txt' file plus every '...RawData_<plate>.eds.txt' file. Plate names in the results must match the RawData file names."),
        checkboxInput(ns("use_custom_metadata"), "Use custom metadata file", value = FALSE),
        conditionalPanel(
          condition = "input.use_custom_metadata == true",
          ns = ns,
          fileInput(ns("file_metadata"), "Plate meta (.xlsx)", accept = ".xlsx"),
          helpText("Optional. Name the file like its raw-data file (e.g. '...RawData_<plate>.eds.xlsx' for '...RawData_<plate>.eds.txt') so annotations are assigned to the correct plate (preferred for multi-plate). Or name it like the analysis-results file for a single plate. If it can't be matched by name, annotations are applied to all wells. Each sheet name becomes a well-annotation column ('target', 'ligand' required; optionally 'conc', 'smiles'); fill each sheet with 384-well values in range B2:Y17."),
          helpText("If no file is provided, each well is treated as its own protein/ligand combination.")
        ),
        radioButtons(ns("tm_type"), "Tm Calculation",
                     choices = c("Tm (Derivative)" = "tm_d", "Tm (Boltzmann)" = "tm_b"),
                     selected = "tm_d", inline = FALSE),
        hr(),
        actionButton(ns("submit"), "Submit", icon = icon("play"), class = "btn-default"),
        br(), br(),
        uiOutput(ns("ref_selectors")),
        hr(),
        downloadButton(ns("download_data"), "Download Data", icon = icon("table")),
        br(), br(),
        downloadButton(ns("download_plots"), "Download Plots", icon = icon("file-pdf"))
      ),
      box(title = "Plot Settings", width = NULL, status = "primary", solidHeader = TRUE,
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
      box(title = "Results", width = NULL, status = "primary", solidHeader = TRUE,
        uiOutput(ns("results_tabs"))
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
      used_fallback = FALSE,
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
          res <- read_analysis_results(res_file) %>%
            janitor::clean_names() %>% 
            janitor::remove_empty('cols') %>% 
            dplyr::rename(plate = experiment_file_name) %>%
            mutate(tm = .data[[input$tm_type]],
                   well = transform_well_style(well),
                   plate = str_remove(plate, "\\.eds$") )
          
          incProgress(0.2, message = "Reading metadata...")
          # Metadata is optional. If a custom file is provided, join it;
          # otherwise treat each well as its own protein/ligand combination.
          meta <- NULL
          used_fallback <- FALSE
          if (input$use_custom_metadata) {
            if (is.null(input$file_metadata)) {
              showNotification("No plate metadata file provided - each well will be treated as its own protein/ligand combination", 
                               type = "warning", duration = 8)
              used_fallback <- TRUE
} else {
              meta_files <- input$file_metadata$datapath %>% 
                set_names(., input$file_metadata$name)
              meta <- meta_files %>% 
                map(read_metadata) %>% 
                list_rbind(., names_to="plate")
              if(!all(c("target", "ligand") %in% colnames(meta))) {
                stop("Custom metadata must contain 'target' and 'ligand' columns")
              }
              
              # Align the metadata's plate name with the data plate by matching the
              # metadata file name to the data files:
              #   - PRIMARY:   match a raw-data file (e.g. '...RawData_<plate>.eds.xlsx'
              #               for '...RawData_<plate>.eds.txt') -> assign that plate.
              #   - SECONDARY: match the analysis-results file (single plate) -> assign
              #               the unique res$plate.
              #   - otherwise: leave plate unchanged -> the join below falls back to
              #               a well-only match (annotations shared across all plates).
              norm_stem <- function(nm) {
                tolower(nm) %>%
                  str_remove("\\.eds\\.txt$") %>% str_remove("\\.txt$") %>%
                  str_remove("\\.eds\\.xlsx$") %>% str_remove("\\.xlsx$") %>%
                  str_remove("\\.eds$")
              }
              raw_names     <- names(all_files)[str_detect(names(all_files), "\\.eds\\.txt$")]
              analysis_name <- names(res_file)
              raw_file_to_plate <- function(raw_name) {
                for (p in unique(res$plate)) {
                  if (str_detect(raw_name, fixed(p))) return(p)
                }
                NA_character_
              }
              
              meta_target_plate <- setNames(character(length(input$file_metadata$name)),
                                            input$file_metadata$name)
              for (m in input$file_metadata$name) {
                ms <- norm_stem(m)
                hit_raw <- raw_names[str_detect(norm_stem(raw_names), fixed(ms))]
                if (length(hit_raw) >= 1) {
                  meta_target_plate[m] <- raw_file_to_plate(hit_raw[1])
                } else if (norm_stem(analysis_name) == ms) {
                  ups <- unique(res$plate)
                  meta_target_plate[m] <- if (length(ups) == 1) ups else NA_character_
                } else {
                  meta_target_plate[m] <- NA_character_
                }
              }
              stem_to_plate <- setNames(meta_target_plate, norm_stem(names(meta_target_plate)))
              meta <- meta %>% mutate(
                plate = dplyr::coalesce(unname(stem_to_plate[ norm_stem(plate) ]), plate)
              )
              
              # Inform the user which matching path was used.
              for (m in input$file_metadata$name) {
                tp <- meta_target_plate[m]
                if (!is.na(tp)) {
                  showNotification(paste0("Metadata '", m, "' matched by file name and assigned to plate '", tp, "'."),
                                   type = "message", duration = 8)
                } else {
                  showNotification(paste0("Metadata '", m, "' not matched by file name; ",
                                   "annotations applied to all wells (by well only)."),
                                   type = "warning", duration = 8)
                }
              }
            }
          } else {
            used_fallback <- TRUE
          }
          
          if (used_fallback) {
            # Keep all raw analysis-result columns (matching the metadata mode).
            # Use the exported Protein/Ligand annotations when present; fall back
            # to the well label only where a value is missing/empty, so that any
            # partial annotation carried by the QuantStudio results export is kept.
            merged <- res
            if ("protein" %in% colnames(merged)) {
              merged <- merged %>% mutate(target = if_else(is.na(protein) | trimws(protein) == "", well, protein))
            } else {
              merged <- merged %>% mutate(target = well)
            }
            if ("ligand" %in% colnames(merged)) {
              merged <- merged %>% mutate(ligand = if_else(is.na(ligand) | trimws(ligand) == "", well, ligand))
            } else {
              merged <- merged %>% mutate(ligand = well)
            }
            meta <- merged %>%
              dplyr::select(any_of(c("plate", "well", "target", "ligand", "conc", "smiles"))) %>%
              distinct()
          } else {
            merged <- res %>% 
              inner_join(meta, by = c("plate", "well"))
            if (nrow(merged) == 0) {
              showNotification("Metadata plate name did not match the experiment plate; matched annotations by well only.",
                               type = "warning", duration = 8)
              # The metadata xlsx filename rarely matches the experiment plate
              # name, so fall back to joining annotations by well only. This
              # applies the well-layout annotations to every plate.
              merged <- res %>% 
                inner_join(meta %>% select(-any_of("plate")), by = "well")
            }
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
            map(read_fluorescence_table) %>%
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
          ref_ligand_default <- if (used_fallback) {
            ""  # per-well fallback: no meaningful blank/DMSO reference
          } else {
            find_best_match(ligand_choices, c("blank", "dmso", "pbs"))
          }
          
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
          rv$used_fallback <- used_fallback
          
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
      if (rv$used_fallback) return(NULL)
      choices <- data()$ligand_choices
      sel <- if (is.null(rv$auto_ref_ligand) || rv$auto_ref_ligand == "" ||
                 !rv$auto_ref_ligand %in% choices) "" else rv$auto_ref_ligand
      tagList(
        selectInput(ns("ref_ligand_input"), "Reference Ligand",
                    choices = c("(None)" = "", choices),
                    selected = sel)
      )
    })
    
    # =========================================================================
    # Results tabs (Delta TM omitted in per-well fallback mode)
    # =========================================================================
    output$results_tabs <- renderUI({
      req(data())
      
      if (rv$used_fallback) {
        tabsetPanel(
          tabPanel("Data Preview",
            DT::dataTableOutput(ns("data_preview"))
          ),
          tabPanel("TM Scatter",
            uiOutput(ns("plot_tm_wrapper"))
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
      } else {
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
      }
    })
    
    # =========================================================================
    # Calculate d_tm reactively when ref changes
    # =========================================================================
    df_with_dtm <- reactive({
      req(data())
      
      if (rv$used_fallback) {
        rv$has_ref_ligand <- FALSE
        rv$ref_ligand_val <- ""
        return(data()$merged %>% dplyr::select(!any_of("d_tm")))
      }
      
      req(input$ref_ligand_input)
      
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
      req(df_with_dtm(), !rv$used_fallback)
      
      df <- df_with_dtm()
      meta <- data()$meta
      
      df %>%
        pivot_wider(id_cols = unique(any_of(c("plate", setdiff(colnames(meta), c("well", "target"))))),
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
      if (nrow(deri) == 0) validate("No derivatives could be computed")
      
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
        sheets <- list(
          "Merged Data" = df_with_dtm(),
          "Metadata" = data()$meta,
          "Raw Fluorescence" = data()$mc_tidy,
          "1st Derivative" = data()$derivatives
        )
        if (!rv$used_fallback) {
          sheets[["Wide Format"]] <- df_wide()
        }
        writexl::write_xlsx(sheets, file)
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
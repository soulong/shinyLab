# shinyLab - DSF Module
# Differential Scanning Fluorimetry analysis module
#
# Layout: the left column holds the upload / reference ligand / download
# controls; the right side has two tabs - "Data" (result tables) and "Plots"
# (interactive plotly curves, one row per group showing raw fluorescence and
# its derivative side by side).

library(patchwork)
library(mgcv)
library(gratia)
# library(rcdk)
library(DT)
library(zip)
library(janitor)
library(plotly)

# =============================================================================
# Reader / Helper Functions
# =============================================================================

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

# Find the Tm column requested by the user ("tm_d" / "tm_b") in a cleaned
# AnalysisResults table. Exports name this column differently between versions
# ("Tm (Derivative)" -> tm_derivative, "Tm D" -> tm_d, ...), so match flexibly
# instead of requiring one exact name. The result carries a "matched"
# attribute: "exact", "fuzzy" or "fallback".
resolve_tm_column <- function(df, choice) {
  cols <- colnames(df)
  exact <- if (identical(choice, "tm_d")) c("tm_d", "tm_derivative") else c("tm_b", "tm_boltzmann")
  hit <- exact[exact %in% cols]
  if (length(hit) > 0) return(structure(hit[1], matched = "exact"))

  pat <- if (identical(choice, "tm_d")) "deriv" else "boltz"
  fuzzy <- cols[str_detect(cols, fixed("tm", ignore_case = TRUE)) &
                  str_detect(cols, regex(pat, ignore_case = TRUE))]
  if (length(fuzzy) > 0) return(structure(fuzzy[1], matched = "fuzzy"))

  tm_like <- cols[str_detect(cols, regex("^tm", ignore_case = TRUE))]
  if (length(tm_like) > 0) return(structure(tm_like[1], matched = "fallback"))

  stop("No Tm column found in the AnalysisResults file (expected a column like 'Tm (Derivative)' or 'Tm (Boltzmann)').",
       call. = FALSE)
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

# All rows of a curve table belonging to one group (NA-safe, %in% based).
group_curve_df <- function(tbl, gvar, g) {
  if (!gvar %in% colnames(tbl)) return(tbl[0, , drop = FALSE])
  tbl[!is.na(tbl[[gvar]]) & tbl[[gvar]] %in% g, , drop = FALSE]
}

# Reference-ligand overlay rows for a group; NULL when the overlay is off,
# the group is the reference itself (grouped by ligand), or nothing matches.
ref_subset_df <- function(tbl, gvar, g, overlay_on, ref_val) {
  if (!isTRUE(overlay_on) || is.null(ref_val) || !nzchar(ref_val)) return(NULL)
  if (identical(gvar, "ligand") && identical(g, ref_val)) return(NULL)
  if (!"ligand" %in% colnames(tbl)) return(NULL)
  out <- tbl[!is.na(tbl$ligand) & tbl$ligand %in% ref_val, , drop = FALSE]
  if (nrow(out) == 0) NULL else out
}

# Build the melt-curve ggplot used for both the interactive plotly view and
# the vector-PDF export. `data` must contain temperature, well and the y_var
# column ("fluorescence" or "derivative").
build_curve_ggplot <- function(data, y_var = c("fluorescence", "derivative"),
                               color_by = "", temp_range = NULL, ref_data = NULL) {
  y_var <- match.arg(y_var)
  y_lab <- if (y_var == "derivative") "dF/dT" else "Fluorescence"

  if (nrow(data) == 0 && (is.null(ref_data) || nrow(ref_data) == 0)) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, label = "No data for this group", size = 4) +
             theme_void())
  }

  has_color <- !is.null(color_by) && !identical(color_by, "") &&
               color_by %in% colnames(data)

  p <- ggplot()

  if (!is.null(ref_data) && nrow(ref_data) > 0) {
    p <- p + geom_path(
      data = ref_data,
      mapping = aes(temperature, .data[[y_var]], group = well,
                    text = paste0(well, " [ref]<br>", temperature, " \u00b0C<br>",
                                  signif(.data[[y_var]], 4))),
      colour = "grey60", linewidth = 0.4, inherit.aes = FALSE
    )
  }

  if (nrow(data) > 0) {
    main_aes <- if (has_color) {
      aes(temperature, .data[[y_var]], group = well, colour = .data[[color_by]],
          text = paste0(well, "<br>", temperature, " \u00b0C<br>", signif(.data[[y_var]], 4)))
    } else {
      aes(temperature, .data[[y_var]], group = well,
          text = paste0(well, "<br>", temperature, " \u00b0C<br>", signif(.data[[y_var]], 4)))
    }

    p <- p + geom_path(data = data, mapping = main_aes,
                       colour = if (has_color) NULL else "darkred",
                       linewidth = 0.6, na.rm = TRUE)

    if (has_color) {
      p <- p + if (is.numeric(data[[color_by]])) {
        scale_colour_viridis_c(option = "C", na.value = "darkred", name = color_by)
      } else {
        scale_colour_viridis_d(option = "C", na.value = "darkred", name = color_by)
      }
    }
  }

  p <- p + theme_bw(base_size = 11) +
    theme(legend.position = "top") +
    labs(x = "Temperature (\u00b0C)", y = y_lab)

  if (!is.null(temp_range) && length(temp_range) == 2 && !any(is.na(temp_range))) {
    p <- p + coord_cartesian(xlim = temp_range)
  }

  p
}

# "Ligand X | 6 wells | median Tm = 52.3 C | median DTm = +2.1 C" for a row header.
group_header_label <- function(df, gvar, g) {
  n_wells <- 0L
  tm_med <- NA_real_
  dtm_med <- NA_real_
  has_dtm_col <- "d_tm" %in% colnames(df)
  if (gvar %in% colnames(df)) {
    idx <- !is.na(df[[gvar]]) & df[[gvar]] %in% g
    n_wells <- sum(idx)
    tm_med <- suppressWarnings(median(df$tm[idx], na.rm = TRUE))
    if (has_dtm_col) dtm_med <- suppressWarnings(median(df$d_tm[idx], na.rm = TRUE))
  }
  dtm_part <- if (has_dtm_col && !is.na(dtm_med)) sprintf("  |  median \u0394Tm = %+.1f \u00b0C", dtm_med) else ""
  sprintf("%s  |  %d wells  |  median Tm = %.1f \u00b0C%s", g, n_wells, tm_med, dtm_part)
}

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- substr(x, 1, 80)
  if (!nzchar(x)) "group" else x
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
        downloadButton(ns("download_data"), "Download Data (xlsx)", icon = icon("table")),
        br(), br(),
        downloadButton(ns("download_plots"), "Download Plots (PDF zip)", icon = icon("file-pdf"))
      ),
      box(title = "Plot Settings", width = NULL, status = "primary", solidHeader = TRUE,
        collapsible = TRUE,
        selectInput(ns("row_group"), "Group curves by (one row per group)",
                    choices = c("ligand", "target", "well"), selected = "ligand"),
        sliderInput(ns("rows_per_page"), "Rows per page", min = 1, max = 8, value = 3, step = 1),
        hr(),
        selectInput(ns("color_var"), "Color curves by", choices = c("(single color)" = ""), selected = ""),
        sliderInput(ns("temp_range"), "Temperature range (\u00b0C)", min = 20, max = 100,
                    value = c(20, 100), step = 1),
        checkboxInput(ns("overlay_ref"), "Overlay reference ligand curves (grey)", value = TRUE),
        hr(),
        radioButtons(ns("pdf_scope"), "PDF export scope",
                     choices = c("All groups (one page per group)" = "all",
                                 "Current page only" = "page"),
                     selected = "all")
      )
    ),
    column(width = 9,
      box(title = "Results", width = NULL, status = "primary", solidHeader = TRUE,
        uiOutput(ns("results_ui"))
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
      used_fallback = TRUE,
      auto_ref_ligand = "",
      page = 1L
    )

    # =========================================================================
    # Process data on submit
    # =========================================================================
    data <- eventReactive(input$submit, {
      if (is.null(input$file_analysis)) {
        showNotification("Please upload the QuantStudio export files first.",
                         type = "error", duration = 8)
        return(NULL)
      }

      tryCatch({
        withProgress(message = 'Processing data...', value = 0, {
          # Read all upload/config inputs once at submit time so later changes
          # to them don't silently invalidate downstream reactives.
          isolate({
            incProgress(0.05, message = "Reading analysis results...")

            all_files <- input$file_analysis$datapath %>%
              set_names(., nm = input$file_analysis$name)

            # ---- Analysis results ------------------------------------------
            res_idx <- which(str_detect(names(all_files), "AnalysisResults"))
            if (length(res_idx) == 0) {
              stop("No '*AnalysisResults*' file found among the uploaded files.")
            }
            if (length(res_idx) > 1) {
              stop("Multiple AnalysisResults files uploaded - please upload only one.")
            }
            res_name <- names(all_files)[res_idx]
            res <- read_analysis_results(all_files[[res_idx]])
            if (is.null(res)) {
              stop("Could not parse the AnalysisResults file (no '#'-prefixed header line found).")
            }
            res <- res %>%
              janitor::clean_names() %>%
              janitor::remove_empty('cols')

            if (!"well" %in% colnames(res)) {
              stop("The AnalysisResults file has no 'Well' column.")
            }

            tm_col <- resolve_tm_column(res, input$tm_type)
            if (identical(attr(tm_col, "matched"), "fallback")) {
              showNotification(paste0("Requested Tm method column not found - falling back to column '", tm_col, "'."),
                               type = "warning", duration = 8)
            }

            # Plate name: QuantStudio exports 'Experiment File Name'; fall back
            # to the results file name (e.g. '<plate>.eds_AnalysisResults.txt').
            if ("experiment_file_name" %in% colnames(res)) {
              res <- res %>% dplyr::rename(plate = experiment_file_name)
            } else {
              fallback_plate <- res_name %>%
                str_remove("_?AnalysisResults\\.txt$") %>% str_remove("\\.eds$")
              showNotification(paste0("AnalysisResults has no 'Experiment File Name' column - using '", fallback_plate, "' as the plate name."),
                               type = "warning", duration = 8)
              res <- res %>% mutate(plate = fallback_plate)
            }

            res <- res %>%
              mutate(tm = suppressWarnings(as.numeric(.data[[tm_col]])),
                     well = transform_well_style(well),
                     plate = str_remove(plate, "\\.eds$"))

            if (all(is.na(res$tm))) {
              stop(paste0("Tm column '", tm_col, "' could not be parsed as numeric."))
            }

            incProgress(0.15, message = "Reading metadata...")
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
                  list_rbind(., names_to = "plate")
                if (!all(c("target", "ligand") %in% colnames(meta))) {
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
                analysis_name <- res_name
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
              # Custom metadata provides target/ligand(/conc/smiles); drop the
              # export's own annotation columns first so the join can't collide
              # into suffixed columns (ligand.x / ligand.y).
              merged <- res %>%
                dplyr::select(-any_of(c("protein", "target", "ligand", "conc", "smiles"))) %>%
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
              { if ("conc" %in% colnames(.)) mutate(., conc = suppressWarnings(as.numeric(conc))) else . } %>%
              relocate(tm, .after = last_col())

            incProgress(0.2, message = "Processing melting curve...")
            # Read raw data files with stricter plate name matching
            match_mc_files <- c()
            for (x in unique(res$plate)) {
              pattern <- paste0("RawData_", x, "\\.eds\\.txt$")
              match_idx <- which(str_detect(names(all_files), pattern))
              if (length(match_idx) == 0) {
                # Fallback: match plate name before .eds.txt
                match_idx <- which(str_detect(names(all_files),
                                              paste0(x, "\\.eds\\.txt$")))
              }
              if (length(match_idx) >= 1) {
                matched <- all_files[match_idx[1]]
                names(matched) <- x
                match_mc_files <- c(match_mc_files, matched)
              }
            }

            missing_plates <- setdiff(unique(res$plate), names(match_mc_files))
            if (length(missing_plates) > 0) {
              showNotification(paste0("No raw data file found for plate(s): ",
                                      paste(missing_plates, collapse = ", "),
                                      " - those plates are skipped in the curve plots."),
                               type = "warning", duration = 10)
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

            ligand_choices <- sort(unique(meta$ligand))

            # Auto-detect the reference ligand only when it really looks like a
            # blank/DMSO/PBS control; Delta Tm against an arbitrary ligand would
            # be meaningless.
            auto_ref <- ""
            if (!used_fallback) {
              cand <- find_best_match(ligand_choices, c("blank", "dmso", "pbs"))
              if (length(cand) == 1 && grepl("blank|dmso|pbs", tolower(cand))) {
                auto_ref <- cand
              } else {
                showNotification("No blank/DMSO/PBS-like ligand found - \u0394Tm stays empty until a Reference Ligand is selected.",
                                 type = "warning", duration = 8)
              }
            }

            incProgress(0.25, detail = "Merging with raw data...")
            annot <- merged %>%
              dplyr::select(-any_of(c("reading", "temperature", "fluorescence"))) %>%
              distinct()
            mc_tidy <- mc %>%
              inner_join(annot, by = c("plate", "well"))

            incProgress(0.3, detail = "Computing derivatives...")
            derivatives <- compute_derivatives(mc_tidy) %>%
              left_join(annot, by = c("plate", "well", "target", "ligand"))

            incProgress(0.35, detail = "Done!")

            rv$used_fallback <- used_fallback
            rv$auto_ref_ligand <- auto_ref
            rv$page <- 1L

            showNotification(
              sprintf("Processed %d wells on %d plate(s); %d wells have curve data.",
                      nrow(merged), dplyr::n_distinct(merged$plate),
                      nrow(dplyr::distinct(mc_tidy, plate, well))),
              type = "message", duration = 8)

            list(merged = merged,
                 mc_tidy = mc_tidy,
                 meta = meta,
                 derivatives = derivatives,
                 ligand_choices = ligand_choices)
          })
        })
      }, error = function(e) {
        showNotification(paste("Error processing data:", e$message),
                         type = "error", duration = 10)
        return(NULL)
      })
    })

    # A reference ligand is actually usable (set, exists, not in fallback mode).
    has_ref <- reactive({
      req(data())
      if (rv$used_fallback) return(FALSE)
      r <- input$ref_ligand_input
      !is.null(r) && nzchar(r) && r %in% data()$merged$ligand
    })

    # =========================================================================
    # Dynamic ref selector (shown after submit)
    # =========================================================================
    output$ref_selectors <- renderUI({
      req(data())
      if (rv$used_fallback) return(NULL)
      choices <- data()$ligand_choices
      sel <- rv$auto_ref_ligand
      if (is.null(sel) || !sel %in% choices) sel <- ""
      tagList(
        selectInput(ns("ref_ligand_input"), "Reference Ligand",
                    choices = c("(None)" = "", choices),
                    selected = sel)
      )
    })

    # =========================================================================
    # Per-well results with dTm (recomputed when the reference ligand changes)
    # =========================================================================
    df_with_dtm <- reactive({
      req(data())
      merged <- data()$merged %>% dplyr::select(-any_of("d_tm"))

      if (rv$used_fallback) return(merged)

      ref_ligand_val <- input$ref_ligand_input
      if (is.null(ref_ligand_val) || ref_ligand_val == "" ||
          !ref_ligand_val %in% merged$ligand) {
        return(merged %>% mutate(d_tm = NA_real_))
      }

      ref <- merged %>%
        filter(ligand == ref_ligand_val) %>%
        reframe(tm_ref = median(tm, na.rm = TRUE), .by = c(plate, target))

      merged %>%
        left_join(ref, by = c("plate", "target")) %>%
        mutate(d_tm = tm - tm_ref) %>%
        dplyr::select(!tm_ref)
    })

    # Median Tm / dTm per plate x target x ligand
    tm_summary <- reactive({
      req(data())
      df <- df_with_dtm()
      grp <- intersect(c("plate", "target", "ligand"), colnames(df))
      df %>%
        group_by(across(all_of(grp))) %>%
        summarise(
          n_wells = dplyr::n(),
          tm_median = round(suppressWarnings(median(tm, na.rm = TRUE)), 2),
          d_tm_median = if ("d_tm" %in% colnames(df)) {
            round(suppressWarnings(median(d_tm, na.rm = TRUE)), 2)
          } else NA_real_,
          .groups = "drop"
        ) %>%
        mutate(across(any_of(c("tm_median", "d_tm_median")),
                      ~ dplyr::if_else(is.nan(.x), NA_real_, .x))) %>%
        arrange(across(all_of(grp)))
    })

    # Wide format: one row per plate x ligand, one dTM_ column per target
    df_wide <- reactive({
      req(data(), !rv$used_fallback)
      df <- df_with_dtm()
      meta_cols <- setdiff(colnames(data()$meta), c("well", "target"))
      id_cols <- intersect(c("plate", meta_cols), colnames(df))
      req(length(id_cols) > 0, "d_tm" %in% colnames(df))

      df %>%
        pivot_wider(id_cols = all_of(id_cols), names_from = target,
                    values_from = d_tm, names_prefix = 'dTM_',
                    values_fn = ~ median(.x, na.rm = TRUE)) %>%
        mutate(across(where(is.numeric), \(x) round(x, 3)))
    })

    # =========================================================================
    # Plot row pagination
    # =========================================================================
    groups <- reactive({
      req(data(), input$row_group)
      gv <- input$row_group
      if (!gv %in% colnames(data()$mc_tidy)) return(character(0))
      sort(unique(stats::na.omit(data()$mc_tidy[[gv]])))
    })

    n_pages <- reactive({
      req(data(), input$rows_per_page)
      max(1L, ceiling(length(groups()) / input$rows_per_page))
    })

    # Keep the current page inside the valid range
    observe({
      req(data())
      pg <- if (is.null(rv$page)) 1L else as.integer(rv$page)
      np <- n_pages()
      if (pg < 1L) rv$page <- 1L else if (pg > np) rv$page <- np
    })

    visible_groups <- reactive({
      g <- groups()
      req(length(g) > 0, input$rows_per_page)
      pg <- min(max(if (is.null(rv$page)) 1L else as.integer(rv$page), 1L), n_pages())
      start <- (pg - 1L) * input$rows_per_page + 1L
      g[start:pmin(start + input$rows_per_page - 1L, length(g))]
    })

    observeEvent(input$prev_page, { rv$page <- max(1L, as.integer(rv$page) - 1L) })
    observeEvent(input$next_page, { rv$page <- min(n_pages(), as.integer(rv$page) + 1L) })
    observeEvent(input$row_group, { rv$page <- 1L })
    observeEvent(input$jump_group, {
      i <- match(input$jump_group, groups())
      if (!is.na(i)) rv$page <- ceiling(i / input$rows_per_page)
    })

    # =========================================================================
    # Update controls after data is loaded
    # =========================================================================
    observeEvent(data(), {
      d <- data()
      rv$page <- 1L

      tr <- suppressWarnings(range(d$mc_tidy$temperature, na.rm = TRUE))
      if (all(is.finite(tr))) {
        updateSliderInput(session, "temp_range",
                          min = floor(tr[1]), max = ceiling(tr[2]),
                          value = c(floor(tr[1]), ceiling(tr[2])))
      }

      color_choices <- c("(single color)" = "",
                         intersect(c("plate", "well", "target", "ligand", "conc"),
                                   colnames(d$mc_tidy)))
      sel <- if ("target" %in% names(color_choices)) "target" else ""
      updateSelectInput(session, "color_var", choices = color_choices, selected = sel)
    })

    # =========================================================================
    # Right side: Data tab + Plots tab
    # =========================================================================
    output$results_ui <- renderUI({
      req(data())

      wide_tab <- if (!rv$used_fallback) {
        tabPanel("Wide (Ligand \u00d7 Target)", br(), DT::dataTableOutput(ns("wide_tbl")))
      } else NULL

      g <- groups()
      tabsetPanel(type = "tabs",
        tabPanel("Data", br(),
          tabsetPanel(type = "tabs",
            tabPanel("Per-well Results", br(), DT::dataTableOutput(ns("data_preview"))),
            tabPanel("Tm Summary", br(), DT::dataTableOutput(ns("tm_summary_tbl"))),
            wide_tab
          )
        ),
        tabPanel("Plots", br(),
          fluidRow(
            column(width = 1, actionButton(ns("prev_page"), "", icon = icon("chevron-left"))),
            column(width = 2, div(style = "text-align: center; padding-top: 10px;",
                                  textOutput(ns("page_indicator")))),
            column(width = 1, actionButton(ns("next_page"), "", icon = icon("chevron-right"))),
            column(width = 5, selectInput(ns("jump_group"), "Jump to group:",
                                          choices = g, selected = if (length(g)) g[1] else NULL))
          ),
          uiOutput(ns("plot_rows"))
        )
      )
    })

    output$page_indicator <- renderText({
      req(data())
      n <- length(groups())
      if (n == 0) return("No groups")
      sprintf("Page %d / %d  (%d groups)", rv$page, n_pages(), n)
    })

    output$plot_rows <- renderUI({
      req(data())
      vg <- visible_groups()
      gvar <- input$row_group
      if (length(vg) == 0) {
        return(helpText("No curves available for the selected grouping."))
      }

      df_ann <- df_with_dtm()

      tagList(lapply(seq_along(vg), function(i) {
        g <- vg[[i]]
        tagList(
          div(style = "border-top: 2px solid #eee; margin-top: 14px; padding-top: 6px;",
              h5(group_header_label(df_ann, gvar, g), style = "margin: 0 0 2px 0;")),
          fluidRow(
            column(width = 6, plotlyOutput(ns(paste0("raw_", i)), height = "380px")),
            column(width = 6, plotlyOutput(ns(paste0("deri_", i)), height = "380px"))
          )
        )
      }))
    })

    # One plotly chart per visible group / y variable. plotly builds these
    # from the same ggplot objects as the PDF export, so both views match.
    curve_plotly_for <- function(g, y_var) {
      req(data())
      gvar <- input$row_group
      src <- if (y_var == "derivative") data()$derivatives else data()$mc_tidy
      df <- group_curve_df(src, gvar, g)
      ref_df <- ref_subset_df(src, gvar, g,
                              isTRUE(input$overlay_ref) && has_ref(),
                              if (has_ref()) input$ref_ligand_input else NULL)
      p <- build_curve_ggplot(df, y_var, color_by = input$color_var,
                              temp_range = input$temp_range, ref_data = ref_df)
      ggplotly(p, tooltip = "text")
    }

    observe({
      req(data())
      vg <- visible_groups()
      for (i in seq_along(vg)) {
        local({
          ii <- i
          gg <- vg[[ii]]
          output[[paste0("raw_", ii)]] <- renderPlotly({ curve_plotly_for(gg, "fluorescence") })
          output[[paste0("deri_", ii)]] <- renderPlotly({ curve_plotly_for(gg, "derivative") })
        })
      }
    })

    # =========================================================================
    # Data tab tables
    # =========================================================================
    output$data_preview <- DT::renderDataTable({
      req(data())
      df_with_dtm() %>%
        mutate(across(any_of(c("tm", "d_tm", "conc")), ~ round(.x, 3)))
    }, options = list(pageLength = 25, scrollX = TRUE))

    output$tm_summary_tbl <- DT::renderDataTable({
      req(data())
      tm_summary()
    }, options = list(pageLength = 25, scrollX = TRUE))

    output$wide_tbl <- DT::renderDataTable({
      req(data(), !rv$used_fallback)
      df_wide()
    }, options = list(pageLength = 25, scrollX = TRUE))

    # =========================================================================
    # Download Data (raw + analysis + plot data)
    # =========================================================================
    output$download_data <- downloadHandler(
      filename = function() {
        paste0('DSF_', Sys.Date(), '_data.xlsx')
      },
      content = function(file) {
        req(data())
        sheets <- list(
          "Per-well Results" = df_with_dtm(),
          "Tm Summary" = tm_summary(),
          "Metadata" = data()$meta,
          "Raw Fluorescence" = data()$mc_tidy,
          "1st Derivative" = data()$derivatives
        )
        if (!rv$used_fallback) {
          sheets[["Wide (Ligand x Target)"]] <- df_wide()
        }
        writexl::write_xlsx(sheets, file)
      }
    )

    # =========================================================================
    # Download Plots: editable vector PDF, one page per group
    # (raw fluorescence + derivative side by side), zipped
    # =========================================================================
    build_pdf_page <- function(g, gvar, color_by, temp_range, overlay_on, ref_val) {
      raw_df <- group_curve_df(data()$mc_tidy, gvar, g)
      deri_df <- group_curve_df(data()$derivatives, gvar, g)
      ref_raw <- ref_subset_df(data()$mc_tidy, gvar, g, overlay_on, ref_val)
      ref_deri <- ref_subset_df(data()$derivatives, gvar, g, overlay_on, ref_val)
      p1 <- build_curve_ggplot(raw_df, "fluorescence", color_by, temp_range, ref_raw)
      p2 <- build_curve_ggplot(deri_df, "derivative", color_by, temp_range, ref_deri)
      p1 + p2 + patchwork::plot_annotation(
        title = as.character(g),
        theme = theme_bw(base_size = 14) +
          theme(plot.title = element_text(face = "bold"))
      )
    }

    output$download_plots <- downloadHandler(
      filename = function() {
        paste0('DSF_', Sys.Date(), '_plots.zip')
      },
      content = function(file) {
        req(data())

        gs <- if (identical(input$pdf_scope, "page")) visible_groups() else groups()
        if (length(gs) == 0) {
          showNotification("No groups available to plot.", type = "error", duration = 8)
          return(NULL)
        }

        gvar <- input$row_group
        color_by <- input$color_var
        tr <- input$temp_range
        overlay_on <- isTRUE(input$overlay_ref) && has_ref()
        ref_val <- if (has_ref()) input$ref_ligand_input else NULL

        plot_dir <- tempfile(pattern = "dsf_plots_")
        dir.create(plot_dir, recursive = TRUE)

        tryCatch({
          pages <- vector("list", length(gs))
          withProgress(message = 'Generating vector PDFs...', value = 0, {
            for (i in seq_along(gs)) {
              pages[[i]] <- build_pdf_page(gs[[i]], gvar, color_by, tr, overlay_on, ref_val)
              incProgress(0.9 / length(gs),
                          detail = sprintf("%s (%d/%d)", gs[[i]], i, length(gs)))
            }
          })

          dev_fn <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

          # Combined multi-page PDF (one page per group) ...
          combined <- file.path(plot_dir, "000_all_groups.pdf")
          dev_fn(combined, width = 14, height = 6)
          tryCatch(
            for (p in pages) print(p),
            error = function(e) showNotification(paste("Some PDF pages failed:", e$message),
                                                 type = "warning", duration = 8),
            finally = grDevices::dev.off()
          )

          # ... plus one standalone PDF per group.
          for (i in seq_along(gs)) {
            ggsave(file.path(plot_dir, sprintf("%03d_%s.pdf", i, sanitize_filename(gs[[i]]))),
                   plot = pages[[i]], width = 14, height = 6, device = dev_fn)
          }

          pdf_files <- list.files(plot_dir, pattern = "\\.pdf$")
          zip::zip(zipfile = file, files = pdf_files, root = plot_dir)
        }, error = function(e) {
          showNotification(paste("Error creating zip:", e$message),
                           type = "error", duration = 10)
        })
      },
      contentType = "application/zip"
    )
  })
}

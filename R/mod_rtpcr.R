# shinyBioTools - rtPCR Module
# QuantStudio-style qPCR export: CT table + melt curves, ddCt analysis with
# user-selectable reference target / sample, static ggplot plots, and
# organized downloads (data .xlsx, plots vector-PDF .zip).

library(readxl)
library(writexl)
library(janitor)
library(zip)

rtPCRUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(3,
      box(title = "Settings", width = NULL, status = "primary", solidHeader = TRUE,
        fileInput(ns("userfile"), "Upload qPCR export (.xlsx)", accept = c(".xls", ".xlsx")),
        helpText("QuantStudio-style export.",
                 hint_tip("The reader scans sheet 'Results' (rows 25-51) for the header row that starts with 'Well'. Columns used: 2 = well, 4 = sample name, 5 = target name, 13 = CT. Optional sheet 'Melt Curve Raw' (same header row) uses columns 2 = well, 4-7 = target / temperature / fluorescence / derivative.")),
        uiOutput(ns("ref_ui")),
        checkboxInput(ns("na.do"),
                      label = tagList("Auto handle NA",
                                      hint_tip("Impute missing CT values with the group mean plus N(0, 0.5) noise so replicates stay analysable.")),
                      value = TRUE),
        radioButtons(ns("show_data"), "Display Data",
                     choices = c('ct' = 'ct', 'dCt' = 'd_ct', 'ddCt' = 'dd_ct',
                                 '2^-ddCt' = 'dd_ct_2n', '2^-ddCt (norm)' = 'dd_ct_2n_norm'),
                     selected = 'dd_ct_2n_norm'),
        div(style = "display: flex; justify-content: space-between; align-items: center; gap: 4px; flex-wrap: nowrap;",
          actionButton(ns("apply"), "Analyze", icon = icon("play"), class = "btn-default btn-sm",
                       title = "Run the ddCt analysis on the selected samples/targets with the chosen references."),
          div(style = "display: flex; gap: 4px; flex-wrap: nowrap;",
            downloadButton(ns("download_data"), "Data", icon = icon("table"), class = "btn-default btn-sm",
                           title = "Analysis result, plot data (mean +/- sd), raw CT and raw melt curves as one .xlsx, one sheet each."),
            downloadButton(ns("download_plots"), "Plots", icon = icon("file-pdf"), class = "btn-default btn-sm",
                           title = "All plots as editable vector PDFs (bar plots, melt curves), zipped.")))
      ),
      box(title = "Plot Settings", width = NULL, status = "primary",
          solidHeader = TRUE, collapsible = TRUE,
        sliderInput(ns("sz_bar"), "Label size", 8, 20, 12, 1),
        sliderInput(ns("ang_bar"), "X angle", 0, 90, 45, 15),
        helpText("Applies to the bar plots.",
                 hint_tip("Plots are static previews; the Plots download contains the same charts as editable vector PDFs."))
      )
    ),
    column(9,
      box(title = "Results", width = NULL, status = "primary", solidHeader = TRUE,
        tabsetPanel(
          tabPanel("Bar Plots", br(),
            h5("Per target", hint_tip("One panel per target; the y axis shows the selected Display Data metric; error bars are mean +/- 1 sd.")),
            plotOutput(ns("facet_plot"), height = "440px"),
            h5("All targets combined", hint_tip("All targets side by side (dodged) in one panel for direct comparison.")),
            plotOutput(ns("combine_plot"), height = "440px")
          ),
          tabPanel("Melt Curves", br(),
            h5("By target", hint_tip("Melt-curve derivative (-dF/dT); each curve is one well, colored by sample. The peak marks the melting temperature (Tm).")),
            plotOutput(ns("plot_mc_target"), height = "420px"),
            h5("By target x sample", hint_tip("One panel per target x sample combination; each curve is one well.")),
            plotOutput(ns("plot_mc_sample"), height = "560px")
          ),
          tabPanel("Data", br(), DT::dataTableOutput(ns("table_summary")))
        )
      )
    )
  )
}

rtPCRServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    rv <- reactiveValues(ct = NULL, mc = NULL, result = NULL)

    # Read uploaded file ------------------------------------------------------
    observeEvent(input$userfile, {
      req(input$userfile)

      new_path <- paste0(input$userfile$datapath, ".xlsx")
      if (!file.rename(input$userfile$datapath, new_path)) {
        showNotification("Failed to rename uploaded file", type = "error")
        return()
      }

      row_x <- NULL
      ct <- NULL
      for (i in 24:50) {
        ct_try <- tryCatch(
          read_excel(new_path, sheet = 'Results', skip = i),
          error = function(e) NULL
        )
        if (is.null(ct_try) || ncol(ct_try) == 0) break
        if (identical(colnames(ct_try)[1], 'Well')) {
          row_x <- i
          ct <- ct_try
          break
        }
      }

      if (is.null(row_x)) {
        showNotification("Could not find a valid data start row (a 'Well' header) on sheet 'Results' in the file", type = "error")
        return()
      }
      if (ncol(ct) < 13) {
        showNotification(sprintf("The 'Results' table has only %d columns; expected at least 13 (well/sample/target/CT at columns 2/4/5/13).", ncol(ct)),
                         type = "error", duration = 10)
        return()
      }

      ct <- ct %>%
        dplyr::select(2, 4, 5, 13) %>%
        set_names(c("well", "sample", "target", "ct")) %>%
        janitor::remove_empty('rows') %>%
        mutate(ct = suppressWarnings(as.numeric(ct)),
               sample = as.character(sample),
               target = as.character(target))
      rv$ct <- ct

      mc_sheet <- if ('Melt Curve Raw' %in% excel_sheets(new_path)) 'Melt Curve Raw'
                  else if ('Melt Curve Raw Data' %in% excel_sheets(new_path)) 'Melt Curve Raw Data'
                  else NA
      if (!is.na(mc_sheet)) {
        mc <- tryCatch({
          read_excel(new_path, sheet = mc_sheet, skip = row_x) %>%
            { if (ncol(.) >= 7) dplyr::select(., 2, 4:7) else . } %>%
            set_names(c("well", "target", "temperature", "fluorescence", "derivative")) %>%
            janitor::remove_empty('rows') %>%
            mutate(well = as.character(well), target = as.character(target),
                   temperature = suppressWarnings(as.numeric(temperature)),
                   fluorescence = suppressWarnings(as.numeric(fluorescence)),
                   derivative = suppressWarnings(as.numeric(derivative))) %>%
            left_join(ct[, c("well", "sample", "target")], by = c("well", "target")) %>%
            arrange(target, sample, well, temperature)
        }, error = function(e) {
          showNotification(paste("Could not read the melt curve sheet:", e$message),
                           type = "warning", duration = 8)
          NULL
        })
        rv$mc <- mc
        if (is.null(mc)) showNotification("Melt curve plots are disabled (no usable melt curve data).", type = "warning", duration = 8)
      } else {
        rv$mc <- NULL
        showNotification("No melt curve sheet found - melt curve plots are disabled.", type = "warning", duration = 8)
      }

      rv$result <- NULL
      output$ref_ui <- renderUI({
        tagList(
          selectInput(ns("ref_target"),
                      label = tagList("Reference target", hint_tip("The housekeeping target (e.g. ACTB, GAPDH) used to compute dCt = Ct - mean Ct(ref). Must stay in the Targets filter.")),
                      choices = sort(unique(ct$target)), selected = sort(unique(ct$target))[1]),
          selectInput(ns("ref_sample"),
                      label = tagList("Reference sample", hint_tip("The calibrator sample used to compute ddCt; its 2^-ddCt is set to 1 within each target. Must stay in the Samples filter.")),
                      choices = sort(unique(ct$sample)), selected = sort(unique(ct$sample))[1]),
          selectInput(ns("sample"),
                      label = tagList("Samples", hint_tip("Samples to include in the analysis and plots.")),
                      choices = sort(unique(ct$sample)), multiple = TRUE,
                      selected = sort(unique(ct$sample))),
          selectInput(ns("target"),
                      label = tagList("Targets", hint_tip("Targets (assays) to include in the analysis and plots.")),
                      choices = sort(unique(ct$target)), multiple = TRUE,
                      selected = sort(unique(ct$target)))
        )
      })
    })

    # Analyze -----------------------------------------------------------------
    observeEvent(input$apply, {
      req(rv$ct)
      if (is.null(input$sample) || is.null(input$target) ||
          length(input$sample) == 0 || length(input$target) == 0) {
        showNotification("Select at least one sample and one target", type = "warning")
        return()
      }
      req(input$ref_target, input$ref_sample)

      if (!input$ref_target %in% input$target || !input$ref_sample %in% input$sample) {
        showNotification("The reference target/sample must be included in the current selection", type = "warning")
        return()
      }

      ct_filtered <- rv$ct %>%
        dplyr::filter(sample %in% input$sample, target %in% input$target)
      if (isTRUE(input$na.do)) {
        ct_filtered <- ct_filtered %>%
          mutate(ct = ifelse(is.na(ct), rnorm(n(), mean(ct, na.rm = TRUE), sd = 0.5), ct),
                 .by = c(sample, target))
      }
      # Group means are computed AFTER imputation so the reference-target means
      # used for dCt share the same data basis as the Ct values they subtract
      ct_filtered <- ct_filtered %>%
        mutate(mean = mean(ct, na.rm = TRUE), .by = c(sample, target))

      # dCt: relative to the reference (housekeeping) target
      ref_target <- ct_filtered %>%
        dplyr::filter(target == input$ref_target) %>%
        dplyr::distinct(sample, mean) %>%
        dplyr::rename(ref_target_mean = mean)
      d_ct <- ct_filtered %>%
        left_join(ref_target, by = "sample") %>%
        mutate(d_ct = ct - ref_target_mean)

      # ddCt: relative to the reference (calibrator) sample
      ref_sample <- d_ct %>%
        dplyr::filter(sample == input$ref_sample) %>%
        reframe(ref_sample_mean = mean(d_ct, na.rm = TRUE), .by = c(target))
      dd_ct <- d_ct %>%
        left_join(ref_sample, by = "target") %>%
        mutate(dd_ct = d_ct - ref_sample_mean) %>%
        mutate(dd_ct_2n = 2^(-dd_ct),
               dd_ct_2n_mean = mean(dd_ct_2n, na.rm = TRUE), .by = c(sample, target))

      # Normalize the reference sample to 1 within each target
      ref_norm_factor <- dd_ct %>%
        dplyr::filter(sample == input$ref_sample) %>%
        reframe(norm_factor = mean(dd_ct_2n_mean, na.rm = TRUE), .by = c(target))
      dd_ct_final <- dd_ct %>%
        left_join(ref_norm_factor, by = "target") %>%
        mutate(dd_ct_2n_norm = dd_ct_2n / norm_factor) %>%
        dplyr::select(well, sample, target, ct, d_ct, dd_ct, dd_ct_2n, dd_ct_2n_norm)

      rv$result <- dd_ct_final
    })

    # Summary (plot data): mean +/- sd per sample x target of selected metric
    plot_summary <- reactive({
      req(rv$result)
      metric <- input$show_data
      rv$result %>%
        filter(!is.na(.data[[metric]])) %>%
        reframe(n = dplyr::n(), mean = mean(.data[[metric]]), sd = sd(.data[[metric]]),
                .by = c(sample, target))
    })

    # Melt curve data filtered by current selections (rv$mc stays pristine)
    mc_filtered <- reactive({
      req(rv$mc)
      d <- rv$mc
      if (length(input$sample) > 0) d <- dplyr::filter(d, sample %in% input$sample)
      if (length(input$target) > 0) d <- dplyr::filter(d, target %in% input$target)
      d
    })

    # Bar plots -----------------------------------------------------------------
    bar_plot <- function(faceted) {
      req(rv$result)
      metric <- input$show_data
      y_lab <- switch(metric,
                      ct = "Ct", d_ct = "\u0394Ct", dd_ct = "\u0394\u0394Ct",
                      dd_ct_2n = "2^(-\u0394\u0394Ct)", dd_ct_2n_norm = "2^(-\u0394\u0394Ct), normalized")
      p <- ggplot(plot_summary(), aes(sample, mean, fill = target)) +
        geom_col(width = 0.7, position = position_dodge(0.75)) +
        geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                      width = 0.4, position = position_dodge(0.75)) +
        labs(x = "", y = y_lab, fill = "Target") +
        theme_classic(base_size = input$sz_bar) +
        theme(axis.text.x = element_text(angle = input$ang_bar, vjust = 0.5, hjust = 1),
              legend.position = if (faceted) "none" else "top",
              legend.title = element_blank())
      if (faceted) p <- p + facet_wrap(vars(target), scales = "free_y")
      p
    }

    output$facet_plot <- renderPlot(bar_plot(faceted = TRUE), res = 96)
    output$combine_plot <- renderPlot(bar_plot(faceted = FALSE), res = 96)

    # Melt curve plots ----------------------------------------------------------
    mc_plot_target <- function() {
      req(mc_filtered(), nrow(mc_filtered()) > 0)
      ggplot(mc_filtered(), aes(temperature, derivative, color = sample)) +
        geom_line() +
        facet_wrap(vars(target), scales = "free_y") +
        labs(x = "Temperature (\u00b0C)", y = "-d(F)/dT", color = "Sample") +
        theme_minimal(base_size = 12)
    }

    mc_plot_sample <- function() {
      req(mc_filtered(), nrow(mc_filtered()) > 0)
      ggplot(mc_filtered(), aes(temperature, derivative, color = well)) +
        geom_line() +
        facet_grid(rows = vars(target), cols = vars(sample), scales = "free_y") +
        labs(x = "Temperature (\u00b0C)", y = "-d(F)/dT", color = "Well") +
        theme_minimal(base_size = 11) +
        theme(legend.position = "none")
    }

    output$plot_mc_target <- renderPlot(mc_plot_target(), res = 96)
    output$plot_mc_sample <- renderPlot(mc_plot_sample(), res = 96)

    # Data table ----------------------------------------------------------------
    output$table_summary <- DT::renderDataTable({
      req(rv$result)
      DT::datatable(rv$result, rownames = FALSE, options = dt_table_options())
    })

    # Downloads -----------------------------------------------------------------
    output$download_data <- downloadHandler(
      filename = function() paste0('rtPCR_', Sys.Date(), '_data.xlsx'),
      content = function(file) {
        sheets <- list(
          "Analysis Result" = rv$result,
          "Plot Data (mean-sd)" = plot_summary()
        )
        if (!is.null(rv$ct)) sheets[["Raw CT"]] <- rv$ct
        if (!is.null(rv$mc)) sheets[["Melt Curve Raw"]] <- rv$mc
        write_xlsx_sheets(sheets, file)
      }
    )

    output$download_plots <- downloadHandler(
      filename = function() paste0('rtPCR_', Sys.Date(), '_plots.zip'),
      content = function(file) {
        req(rv$result)
        metric <- input$show_data
        y_lab <- switch(metric,
                        ct = "Ct", d_ct = "dCt", dd_ct = "ddCt",
                        dd_ct_2n = "2^-ddCt", dd_ct_2n_norm = "2^-ddCt (norm)")
        plot_dir <- tempfile(pattern = "rtpcr_plots_")
        dir.create(plot_dir, recursive = TRUE)

        # The downloads show the same charts as the on-screen previews
        export_plot_pdf(bar_plot(TRUE), file.path(plot_dir, "bar_plot_faceted.pdf"),
                        width = 8, height = 6)
        export_plot_pdf(bar_plot(FALSE), file.path(plot_dir, "bar_plot_combined.pdf"),
                        width = 9, height = 6)

        if (!is.null(rv$mc) && nrow(mc_filtered()) > 0) {
          export_plot_pdf(mc_plot_target(), file.path(plot_dir, "melt_curve_by_target.pdf"),
                          width = 10, height = 6)
          export_plot_pdf(mc_plot_sample(), file.path(plot_dir, "melt_curve_by_target_sample.pdf"),
                          width = 12, height = 8)
        }

        pdf_files <- list.files(plot_dir, pattern = "\\.pdf$")
        zip::zip(zipfile = file, files = pdf_files, root = plot_dir)
        unlink(plot_dir, recursive = TRUE)
      },
      contentType = "application/zip"
    )
  })
}

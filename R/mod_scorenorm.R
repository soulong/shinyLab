# shinyLab - Score to Normal Module
# Rank-based Inverse Normal Transformation (INT) module
#
# Layout follows the shared module pattern: left "Settings" box (input +
# downloads), right "Results" box with Data / Plots / Summary tabs.

# =============================================================================
# Rank-based Inverse Normal Transformation (INT)
# Preserves original order using van der Waerden method
# =============================================================================

rank_normal_transform <- function(x) {
  n <- length(x)
  if (n == 0) return(numeric(0))

  ranks <- rank(x, ties.method = "average")
  qnorm((ranks - 0.5) / n)
}

# =============================================================================
# UI Function
# =============================================================================

scoreNormUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(title = "Input Data", width = 3, status = "primary", solidHeader = TRUE,
      p("Numbers separated by spaces, commas, tabs, or newlines:",
        hint_tip("Rank-based inverse normal transformation (van der Waerden): each value is replaced by the standard-normal quantile of its average rank, so the output is normally distributed while the original order is preserved.")),
      textAreaInput(ns("input_numbers"), NULL,
        placeholder = "1.2  5.4  3.1  2.8  10.5\n-2  0  1.5",
        height = "200px", width = "100%"),
      br(),
      div(style = "display: flex; justify-content: space-between; align-items: center; gap: 4px; flex-wrap: nowrap;",
        actionButton(ns("transform_btn"), "Transform", icon = icon("magic"), class = "btn-default btn-sm",
                     title = "Rank-based inverse normal transformation (van der Waerden): each value is replaced by the standard-normal quantile of its average rank, so the output is normally distributed while the original order is preserved."),
        div(style = "display: flex; gap: 4px; flex-wrap: nowrap;",
          downloadButton(ns("download_data"), "Data", icon = icon("table"), class = "btn-default btn-sm",
                         title = "Original vs transformed values plus summary statistics as one .xlsx, one sheet each."),
          downloadButton(ns("download_plots"), "Plots", icon = icon("file-pdf"), class = "btn-default btn-sm",
                         title = "Distribution histogram and Q-Q plot as vector PDFs, zipped.")))
    ),
    column(9,
      box(title = "Results", width = NULL, status = "primary", solidHeader = TRUE,
        tabsetPanel(
          tabPanel("Data", br(),
            DT::dataTableOutput(ns("comparison_table"))
          ),
          tabPanel("Plots", br(),
            h5("Distribution: original vs transformed (red curve = fitted normal density)"),
            plotOutput(ns("density_plot"), height = "420px"),
            h5("Normal Q-Q plot of transformed data (points on the line = normal)"),
            plotOutput(ns("qq_plot"), height = "420px")
          ),
          tabPanel("Summary", br(),
            verbatimTextOutput(ns("summary_text"))
          )
        )
      )
    )
  )
}

# =============================================================================
# Server Function
# =============================================================================

scoreNormServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    input_data <- eventReactive(input$transform_btn, {
      text <- input$input_numbers
      if (is.null(text) || trimws(text) == "") return(numeric(0))

      tokens <- unlist(strsplit(text, "[[:space:],]+"))
      tokens <- tokens[tokens != ""]
      nums <- suppressWarnings(as.numeric(tokens))
      invalid <- is.na(nums)
      if (any(invalid)) {
        showNotification(paste("Non-numeric entries ignored:",
          paste(head(tokens[invalid], 5), collapse = ", ")),
          type = "warning", duration = 5)
      }
      return(nums[!invalid])
    })

    transformed_data <- reactive({
      x <- input_data()
      if (length(x) == 0) return(numeric(0))
      rank_normal_transform(x)
    })

    data_table <- reactive({
      x <- input_data()
      y <- transformed_data()
      if (length(x) == 0) return(NULL)
      tibble(Index = seq_along(x), Original = x, Transformed = round(y, 6))
    })

    # ----- Data tab -----
    output$comparison_table <- DT::renderDataTable({
      req(data_table())
      DT::datatable(data_table(), rownames = FALSE, options = dt_table_options())
    })

    # ----- Plots -----
    output$density_plot <- renderPlot({
      x <- input_data()
      y <- transformed_data()
      if (length(x) < 2) return(NULL)

      df <- bind_rows(
        tibble(value = x, data = "Original"),
        tibble(value = y, data = "Transformed")
      )
      groups <- c("Original", "Transformed")
      vals <- list(x, y)
      dens <- bind_rows(lapply(seq_along(groups), function(i) {
        dn <- density(vals[[i]], n = 200)
        tibble(x = dn$x, y = dn$y, data = groups[i])
      }))

      ggplot(df, aes(value, fill = data)) +
        geom_histogram(aes(y = after_stat(density)),
                       position = "identity", alpha = 0.55, bins = 30, color = "white") +
        geom_line(data = dens, aes(x, y, color = data), linewidth = 1, inherit.aes = FALSE) +
        scale_fill_manual(values = c(Original = "#4e79a7", Transformed = "#e15759")) +
        scale_color_manual(values = c(Original = "#4e79a7", Transformed = "#e15759")) +
        labs(x = "Value", y = "Density", fill = NULL, color = "Density") +
        theme_bw(base_size = 12) +
        theme(legend.position = "top")
    }, res = 96)

    output$qq_plot <- renderPlot({
      y <- transformed_data()
      if (length(y) < 3) return(NULL)

      n <- length(y)
      theoretical <- qnorm(ppoints(n))
      observed <- sort(y)

      ggplot(tibble(theoretical, observed), aes(theoretical, observed)) +
        geom_point(color = "#e15759", size = 2, alpha = 0.8) +
        geom_abline(intercept = mean(y), slope = sd(y), linetype = "dashed", color = "grey40") +
        labs(x = "Theoretical normal quantiles", y = "Observed quantiles (transformed)") +
        theme_bw(base_size = 12)
    }, res = 96)

    # ----- Summary tab -----
    output$summary_text <- renderPrint({
      x <- input_data()
      y <- transformed_data()
      if (length(x) == 0) return("No data yet. Enter numbers and click 'Transform to Normal'.")
      cat("Original data:\n")
      print(summary(x))
      cat("\nShapiro-Wilk normality test (n < 5000):\n")
      if (length(x) >= 3 && length(x) <= 5000) print(shapiro.test(x)) else cat("skipped (n out of range)\n")
      cat("\nTransformed data:\n")
      print(summary(y))
      if (length(y) >= 3 && length(y) <= 5000) print(shapiro.test(y))
    })

    # ----- Downloads -----
    output$download_data <- downloadHandler(
      filename = function() paste0('ScoreNorm_', Sys.Date(), '_data.xlsx'),
      content = function(file) {
        x <- input_data()
        y <- transformed_data()
        sheets <- list(
          "Original vs Transformed" = tibble(Index = seq_along(x), Original = x, Transformed = y),
          "Summary" = tibble(
            data = rep(c("Original", "Transformed"), each = 8),
            statistic = rep(c("min", "q1", "median", "mean", "q3", "max", "sd", "n"), 2),
            value = c(min(x), quantile(x, .25), median(x), mean(x), quantile(x, .75), max(x), sd(x), length(x),
                      min(y), quantile(y, .25), median(y), mean(y), quantile(y, .75), max(y), sd(y), length(y))
          )
        )
        write_xlsx_sheets(sheets, file)
      }
    )

    output$download_plots <- downloadHandler(
      filename = function() paste0('ScoreNorm_', Sys.Date(), '_plots.zip'),
      content = function(file) {
        x <- input_data()
        y <- transformed_data()
        req(length(x) >= 2)
        plot_dir <- tempfile(pattern = "scorenorm_plots_")
        dir.create(plot_dir, recursive = TRUE)

        df <- bind_rows(
          tibble(value = x, data = "Original"),
          tibble(value = y, data = "Transformed")
        )
        groups <- c("Original", "Transformed")
        vals <- list(x, y)
        dens <- bind_rows(lapply(seq_along(groups), function(i) {
          dn <- density(vals[[i]], n = 200)
          tibble(x = dn$x, y = dn$y, data = groups[i])
        }))

        p1 <- ggplot(df, aes(value, fill = data)) +
          geom_histogram(aes(y = after_stat(density)), position = "identity",
                         alpha = 0.55, bins = 30, color = "white") +
          geom_line(data = dens, aes(x, y, color = data), linewidth = 1, inherit.aes = FALSE) +
          scale_fill_manual(values = c(Original = "#4e79a7", Transformed = "#e15759")) +
          scale_color_manual(values = c(Original = "#4e79a7", Transformed = "#e15759")) +
          labs(title = "Original vs Transformed", x = "Value", y = "Density", fill = NULL,
               color = "Density") +
          theme_bw(base_size = 12) + theme(legend.position = "top")
        export_plot_pdf(p1, file.path(plot_dir, "distribution_histogram.pdf"), 8, 5)

        n <- length(y)
        theoretical <- qnorm(ppoints(n))
        observed <- sort(y)
        p2 <- ggplot(tibble(theoretical, observed), aes(theoretical, observed)) +
          geom_point(color = "#e15759", size = 2, alpha = 0.8) +
          geom_abline(intercept = mean(y), slope = sd(y), linetype = "dashed", color = "grey40") +
          labs(title = "Q-Q Plot (Transformed)",
               x = "Theoretical normal quantiles", y = "Observed quantiles") +
          theme_bw(base_size = 12)
        export_plot_pdf(p2, file.path(plot_dir, "qq_plot_transformed.pdf"), 6, 6)

        pdf_files <- list.files(plot_dir, pattern = "\\.pdf$")
        zip::zip(zipfile = file, files = pdf_files, root = plot_dir)
        unlink(plot_dir, recursive = TRUE)
      },
      contentType = "application/zip"
    )
  })
}

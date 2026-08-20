# shinyLab - Score to Normal Module
# Rank-based Inverse Normal Transformation (INT) module

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
      p("Enter numbers separated by spaces, commas, tabs, or newlines:"),
      textAreaInput(ns("input_numbers"), NULL,
        placeholder = "1.2  5.4  3.1  2.8  10.5\n-2  0  1.5",
        height = "200px", width = "100%"),
      br(),
      actionButton(ns("transform_btn"), "Transform to Normal", icon = icon("magic"),
                   class = "btn-primary"),
      br(), br(),
      downloadButton(ns("download_data"), "Download Transformed Data (.csv)")
    ),
    column(9,
      box(title = "Original vs Transformed Data", width = NULL, status = "primary", solidHeader = TRUE,
        tableOutput(ns("comparison_table")),
        br(),
        h4("Summary Statistics"),
        verbatimTextOutput(ns("summary_orig")),
        verbatimTextOutput(ns("summary_trans")),
        br(),
        plotOutput(ns("density_plot"))
      )
    )
  )
}

# =============================================================================
# Server Function
# =============================================================================

scoreNormServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

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

    output$comparison_table <- renderTable({
      x <- input_data()
      y <- transformed_data()
      if (length(x) == 0) return(data.frame())

      data.frame(
        Index = 1:length(x),
        Original = x,
        Transformed = round(y, 4)
      )
    }, digits = 6)

    output$summary_orig <- renderPrint({
      if (length(input_data()) == 0) return("No data yet. Enter numbers and click 'Transform to Normal'.")
      summary(input_data())
    })

    output$summary_trans <- renderPrint({
      if (length(transformed_data()) == 0) return("")
      summary(transformed_data())
    })

    output$density_plot <- renderPlot({
      x <- input_data()
      y <- transformed_data()
      if (length(x) < 2) return(NULL)

      par(mfrow = c(1, 2))
      hist(x, main = "Original Data", xlab = "Value", col = "lightblue", border = "white")
      hist(y, main = "Transformed (Normal)", xlab = "Value", col = "lightcoral", border = "white")

      curve(dnorm(x, mean(y), sd(y)), add = TRUE, col = "red", lwd = 2)
    })

    output$download_data <- downloadHandler(
      filename = function() {
        paste("normal_transformed_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        df <- data.frame(
          Original = input_data(),
          Transformed = transformed_data()
        )
        write.csv(df, file, row.names = FALSE)
      }
    )
  })
}
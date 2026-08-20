# shinyBioTools - sgRNA Module

# =============================================================================
# Reverse Complement Helper
# =============================================================================

revComp <- function(seq) {
  vapply(strsplit(chartr("ATCG", "TAGC", seq), ""),
         function(x) paste(rev(x), collapse = ""), "")
}

# =============================================================================
# UI Function
# =============================================================================

sgRNAUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(title = "Settings", width = 3, status = "primary", solidHeader = TRUE,
      textAreaInput(ns("input"), "sgRNA List",
        placeholder = "Paste sgRNA sequences (one per line)\nTwo-column: Name<TAB>Sequence\nOne-column: sequence only, auto-numbered\n\nExample:\ngene1\tCACCGGAGGTCTCCTAGCA\nGATCCGATCGAACTTCGAC",
        height = "200px"),
      p(icon("info-circle"), " Default prefix/suffix are for LentiCRISPRv2. Adjust as needed.",
        style = "font-size: 12px; color: #888;"),
      hr(),
      h5("Forward Primer (5' -> 3': prefix-sgRNA-suffix)"),
      fluidRow(
        column(6, textInput(ns("f5"), "5' prefix", value = "caccG")),
        column(6, textInput(ns("f3"), "3' suffix", value = ""))
      ),
      h5("Reverse Primer (5' -> 3': prefix-revComp-suffix)"),
      fluidRow(
        column(6, textInput(ns("r5"), "5' prefix", value = "aaac")),
        column(6, textInput(ns("r3"), "3' suffix", value = "C"))
      ),
      hr(),
      actionButton(ns("design"), "Design Primers", icon = icon("magic"),
                   class = "btn-default"),
      br(), br(),
      downloadButton(ns("download"), "Download")
    ),
    column(9,
      box(title = "Designed Primers", width = NULL, status = "primary", solidHeader = TRUE,
        tableOutput(ns("primer_table"))
      )
    )
  )
}

# =============================================================================
# Server Function
# =============================================================================

sgRNAServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    rv <- reactiveValues(primers = NULL)

    observeEvent(input$design, {
      req(input$input)

      lines <- trimws(strsplit(input$input, "\n")[[1]])
      lines <- lines[lines != ""]

      if (length(lines) == 0) {
        showNotification("Enter at least one sgRNA sequence", type = "warning")
        return()
      }

      has_tab <- grepl("\t", lines)
      n_unnamed <- sum(!has_tab)
      auto_prefix <- seq_len(n_unnamed)

      prefixes <- character(length(lines))
      sequences <- character(length(lines))

      if (any(has_tab)) {
        parts <- strsplit(lines[has_tab], "\t")
        prefixes[has_tab] <- vapply(parts, `[`, "", 1)
        sequences[has_tab] <- vapply(parts, `[`, "", 2)
      }
      if (any(!has_tab)) {
        prefixes[!has_tab] <- as.character(auto_prefix)
        sequences[!has_tab] <- lines[!has_tab]
      }

      sequences <- trimws(toupper(sequences))
      prefixes <- trimws(prefixes)

      invalid <- !grepl("^[ATCG]+$", sequences)
      if (any(invalid)) {
        msg <- paste("Invalid sequences (non-ATCG) skipped:",
                     paste(sequences[invalid], collapse = ", "))
        showNotification(msg, type = "warning", duration = 10)
        prefixes <- prefixes[!invalid]
        sequences <- sequences[!invalid]
      }

      if (length(sequences) == 0) {
        showNotification("No valid sgRNA sequences found", type = "error")
        return()
      }

      f5 <- input$f5
      f3 <- input$f3
      r5 <- input$r5
      r3 <- input$r3

      rc <- revComp(sequences)

      f_primers <- paste0(f5, sequences, f3)
      r_primers <- paste0(r5, rc, r3)

      rv$primers <- data.frame(
        ID = paste0(rep(prefixes, each = 2), c("_f", "_r")),
        Primer = as.vector(rbind(f_primers, r_primers)),
        stringsAsFactors = FALSE
      )
    })

    output$primer_table <- renderTable({
      req(rv$primers)
      rv$primers
    }, striped = TRUE, hover = TRUE)

    output$download <- downloadHandler(
      filename = function() paste0("sgRNA_primers_", Sys.Date(), ".xlsx"),
      content = function(file) {
        req(rv$primers)
        writexl::write_xlsx(rv$primers, file)
      }
    )
  })
}

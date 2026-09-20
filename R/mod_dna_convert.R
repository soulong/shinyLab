# shinyBioTools - DNA Convert Module
# Batch reverse / complement / reverse-complement of a multi-line DNA list.

# =============================================================================
# Sequence Helpers
# =============================================================================

#' Clean a raw sequence string: keep letters only, uppercase, restrict to
#' IUPAC nucleotide codes (digits, whitespace and other symbols are dropped).
dna_convert_clean <- function(seq) {
  s <- toupper(gsub("[^A-Za-z]", "", seq))
  gsub("[^ACGTURYSWKMBDHVN]", "", s)
}

# Complement lookup (IUPAC ambiguity codes included). rna = TRUE returns an
# RNA alphabet complement (A <-> U) instead of DNA (A <-> T).
complement_map <- function(rna = FALSE) {
  if (rna) {
    c(A = "U", U = "A", T = "A", G = "C", C = "G",
      R = "Y", Y = "R", S = "S", W = "W", K = "M", M = "K",
      B = "V", V = "B", D = "H", H = "D", N = "N")
  } else {
    c(A = "T", T = "A", U = "A", G = "C", C = "G",
      R = "Y", Y = "R", S = "S", W = "W", K = "M", M = "K",
      B = "V", V = "B", D = "H", H = "D", N = "N")
  }
}

seq_complement <- function(seq, rna = FALSE) {
  if (length(seq) == 0) return(character(0))
  m <- complement_map(rna)
  chartr(paste(names(m), collapse = ""), paste(m, collapse = ""), seq)
}

seq_reverse <- function(seq) {
  if (length(seq) == 0) return(character(0))
  vapply(strsplit(seq, ""), function(x) paste(rev(x), collapse = ""), "")
}

# =============================================================================
# UI Function
# =============================================================================

dnaConvertUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(width = 4,
      box(title = "Input", width = NULL, status = "primary", solidHeader = TRUE,
        textAreaInput(ns("seq_input"), "DNA / RNA sequences (one per line)",
          placeholder = "ATGGCGTACGT\nAUGGCGUACGU\nGTACCGTATGGC\nACGTTN",
          height = "220px"),
        checkboxInput(ns("u2t"), "Convert U -> T (treat input as DNA)", value = TRUE),
        checkboxInput(ns("t2u"), "Convert T -> U (output as RNA)", value = FALSE),
        helpText(icon("info-circle"),
                 " One sequence per line. Characters other than IUPAC nucleotide codes ",
                 tags$code("ACGTURYSWKMBDHVN"), " are ignored."),
        hr(),
        downloadButton(ns("download_csv"), "Download CSV", icon = icon("file-csv")),
        br(), br(),
        downloadButton(ns("download_xlsx"), "Download Excel", icon = icon("file-excel"))
      )
    ),
    column(width = 8,
      box(title = "Results", width = NULL, status = "primary", solidHeader = TRUE,
        uiOutput(ns("strip_note")),
        DT::dataTableOutput(ns("result_table"))
      )
    )
  )
}

# =============================================================================
# Server Function
# =============================================================================

dnaConvertServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # U->T and T->U are mutually exclusive: enabling one disables the other.
    observeEvent(input$u2t, {
      if (isTRUE(input$u2t)) updateCheckboxInput(session, "t2u", value = FALSE)
    })
    observeEvent(input$t2u, {
      if (isTRUE(input$t2u)) updateCheckboxInput(session, "u2t", value = FALSE)
    })

    result <- reactive({
      lines <- strsplit(paste(input$seq_input, collapse = "\n"), "\r?\n")[[1]]
      lines <- trimws(lines)
      lines <- lines[nzchar(lines)]

      if (length(lines) == 0) {
        return(tibble(input = character(), reverse = character(),
                      complement = character(), `reverse complement` = character()))
      }

      cleaned <- dna_convert_clean(lines)
      rna_out <- isTRUE(input$t2u)
      if (isTRUE(input$u2t)) seqs <- chartr("U", "T", cleaned) else seqs <- cleaned
      if (rna_out) seqs <- chartr("T", "U", seqs)

      comp <- seq_complement(seqs, rna = rna_out)

      tibble(
        input = lines,
        reverse = seq_reverse(seqs),
        complement = comp,
        `reverse complement` = seq_reverse(comp)
      )
    })

    # Inline note listing lines where non-nucleotide characters were dropped
    # (digits/symbols stripped, or letters outside the IUPAC code set removed).
    output$strip_note <- renderUI({
      res <- result()
      req(nrow(res) > 0)
      raw <- res$input
      letters_only <- toupper(gsub("[^A-Za-z]", "", raw))
      cleaned <- gsub("[^ACGTURYSWKMBDHVN]", "", letters_only)
      n_bad <- sum(nchar(raw) != nchar(letters_only) |
                     nchar(letters_only) != nchar(cleaned))
      if (n_bad == 0) return(NULL)
      tags$div(style = "color: #b8860b; margin-bottom: 6px;",
               icon("exclamation-triangle"),
               sprintf("%d line(s) contained characters outside the IUPAC nucleotide code set - they were removed before conversion.",
                       n_bad))
    })

    output$result_table <- DT::renderDataTable({
      res <- result()
      DT::datatable(res, rownames = FALSE,
                    colnames = c("Input", "Reverse", "Complement", "Reverse Complement"),
                    options = list(pageLength = 25, scrollX = TRUE,
                                   dom = "ltip", searching = FALSE)) %>%
        DT::formatStyle(columns = colnames(res), fontFamily = "Consolas, monospace")
    })

    output$download_csv <- downloadHandler(
      filename = function() paste0("DNA_convert_", Sys.Date(), ".csv"),
      content = function(file) {
        req(nrow(result()) > 0)
        write.csv(result(), file, row.names = FALSE, fileEncoding = "UTF-8", quote = TRUE)
      }
    )

    output$download_xlsx <- downloadHandler(
      filename = function() paste0("DNA_convert_", Sys.Date(), ".xlsx"),
      content = function(file) {
        req(nrow(result()) > 0)
        writexl::write_xlsx(result(), file)
      }
    )
  })
}

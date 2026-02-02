# shinyBioTools - shRNA Module

shRNAUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(title = "Settings", width = 3, status = "primary", solidHeader = TRUE,
      radioButtons(ns("species"), "Species", inline = TRUE,
                   choices = c("Human" = "hs", "Mouse" = "mm")),
      radioButtons(ns("type"), "Input Type",
                   choices = c("Symbol" = "SYMBOL", "Ensembl" = "ENSEMBL",
                               "Entrez" = "ENTREZID")),
      numericInput(ns("number"), "shRNA per gene", 3, 1, 6),
      textAreaInput(ns("input"), "Gene List",
                    placeholder = "One gene per line", height = "150px"),
      hr(),
      actionButton(ns("submit"), "Convert IDs", icon = icon("sync"),
                   class = "btn-info"),
      br(), br(),
      actionButton(ns("design"), "Design Primers", icon = icon("magic"),
                   class = "btn-success"),
      br(), br(),
      downloadButton(ns("download"), "Download")
    ),
    column(9,
      box(title = "Converted Gene IDs", width = 12, collapsible = TRUE,
        tableOutput(ns("gene_table"))
      ),
      box(title = "shRNA Antisense Sequences", width = 12, collapsible = TRUE,
        tableOutput(ns("antisense_table"))
      ),
      box(title = "Designed Primers", width = 12, collapsible = TRUE,
        tableOutput(ns("primer_table"))
      )
    )
  )
}

shRNAServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    rv <- reactiveValues(gene_info = NULL, antisense = NULL, primer = NULL)
    
    # Convert gene IDs
    observeEvent(input$submit, {
      req(input$input)
      genes <- trimws(strsplit(input$input, "\n")[[1]])
      genes <- genes[genes != ""]
      
      if (length(genes) == 0) {
        showNotification("Enter at least one gene", type = "warning")
        return()
      }
      
      rv$gene_info <- id_convert(genes, input$type, input$species)
    })
    
    # Design shRNA primers
    observeEvent(input$design, {
      req(rv$gene_info)
      
      valid <- rv$gene_info[!is.na(rv$gene_info$ENTREZID), ]
      entrez_ids <- unique(valid$ENTREZID)
      
      if (length(entrez_ids) == 0) {
        showNotification("No valid Entrez IDs found", type = "warning")
        return()
      }
      
      symbols <- valid[!duplicated(valid$ENTREZID), c("SYMBOL", "ENTREZID")]
      
      withProgress(message = "Querying splashRNA...", {
        rv$antisense <- splashRNA_batch(
          entrez_ids, input$number,
          symbols$SYMBOL[match(entrez_ids, symbols$ENTREZID)]
        )
        
        if (!is.null(rv$antisense) && nrow(rv$antisense) > 0) {
          rv$primer <- shRNA_primer(rv$antisense)
        }
      })
    })
    
    # Outputs
    output$gene_table <- renderTable({
      req(rv$gene_info)
      rv$gene_info[, c("query", "SYMBOL", "ENTREZID", "ENSEMBL")]
    }, striped = TRUE, hover = TRUE)
    
    output$antisense_table <- renderTable({
      req(rv$antisense); rv$antisense
    }, striped = TRUE, hover = TRUE)
    
    output$primer_table <- renderTable({
      req(rv$primer); rv$primer
    }, striped = TRUE, hover = TRUE)
    
    output$download <- downloadHandler(
      filename = function() paste0("shRNA_", Sys.Date(), ".xlsx"),
      content = function(file) {
        sheets <- list()
        if (!is.null(rv$antisense)) sheets$Antisense <- rv$antisense
        if (!is.null(rv$primer)) sheets$Primers <- rv$primer
        if (length(sheets) > 0) write_xlsx(sheets, file)
      }
    )
  })
}

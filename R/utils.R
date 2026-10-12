# shinyBioTools - Utility Functions

# =============================================================================
# Shared Download / Plot Export Helpers
# =============================================================================

#' Write a named list of data frames to a multi-sheet Excel file
#'
#' Drops NULL / zero-row sheets, and sanitizes sheet names to Excel's 31-char
#' limit with guaranteed uniqueness. An empty `sheets` list raises an error so
#' callers never produce a corrupt file.
#'
#' @param sheets Named list of data frames (NULL entries are skipped).
#' @param file Target .xlsx path.
write_xlsx_sheets <- function(sheets, file) {
  # purrr:: explicitly: XVector (via Biostrings) masks compact in the search path
  sheets <- purrr::compact(sheets)
  sheets <- purrr::keep(sheets, ~ is.data.frame(.x) && nrow(.x) > 0)
  if (length(sheets) == 0) stop("No data to export (all sheets were empty)")

  names(sheets) <- names(sheets) %>%
    map_chr(~ if (is.null(.x) || !nzchar(.x)) "Sheet" else as.character(.x)) %>%
    substr(1, 31)
  dup <- duplicated(names(sheets))
  if (any(dup)) names(sheets)[dup] <- paste0(names(sheets)[dup], "_", seq_len(sum(dup)))

  writexl::write_xlsx(sheets, file)
}

#' Draw a plot object of any supported type onto the current device
#'
#' Handles ggplot, patchwork, gtable/grob (e.g. the ggplot2 heatmap grob),
#' recorded plots (synergyfinder), and plain lists of those.
draw_plot_object <- function(plot) {
  if (is.list(plot) && !inherits(plot, c("ggplot", "gg", "patchwork", "gtable"))) {
    return(invisible(lapply(plot, draw_plot_object)))
  }
  if (inherits(plot, c("ggplot", "patchwork", "gg"))) print(plot)
  else if (inherits(plot, "recordedplot")) grDevices::replayPlot(plot)
  else if (inherits(plot, c("gtable", "grob", "gTree"))) grid::grid.draw(plot)
  else if (is.function(plot)) plot()
  else if (!is.null(plot)) print(plot)
}

#' Save a plot object as a vector PDF (cairo when available)
#'
#' @param plot A ggplot / patchwork / gtable / grob / recorded plot / list.
#' @param file Target .pdf path.
#' @param width,height Page size in inches.
export_plot_pdf <- function(plot, file, width = 10, height = 8) {
  dev_fn <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  dev_fn(file, width = width, height = height)
  tryCatch(draw_plot_object(plot), finally = grDevices::dev.off())
  invisible(file)
}

#' Consistent DT options for data tables across modules
dt_table_options <- function(pageLength = 25) {
  list(pageLength = pageLength, lengthMenu = c(10, 25, 50, 100),
       scrollX = TRUE, searching = TRUE)
}

#' Inline hover hint
#'
#' A small info icon with a Bootstrap tooltip (delegated via app.R, so it also
#' works inside renderUI output). Use for explanations too long to show inline;
#' the native `title` attribute is the no-JS fallback.
#'
#' @param text Hint text shown on hover.
#' @param placement Tooltip position: "top", "bottom", "left" or "right".
hint_tip <- function(text, placement = "top") {
  tags$span(
    icon("info-circle"),
    title = text,
    `data-toggle` = "tooltip",
    `data-placement` = placement,
    `data-container` = "body",
    style = "color:#3c8dbc; cursor:help; margin-left:4px; font-size:13px;"
  )
}

# =============================================================================
# Well Position Helpers
# =============================================================================

#' Transform Well Position Style
#'
#' @description Removes leading zeros from well position strings.
#'   Converts formats like "A01" to "A1".
#'
#' @param wells Character vector of well positions (e.g., "A01", "B12").
#'
#' @return Character vector with leading zeros removed from numeric portion.
#'
#' @examples
#' transform_well_style(c("A01", "B12", "C03"))
#' # Returns: c("A1", "B12", "C3")
#'
#' @export
transform_well_style <- function(wells) sub("^([A-Za-z]+)0+", "\\1", wells)

# =============================================================================
# Metadata Reading
# =============================================================================

#' Read Plate Metadata from Excel File
#'
#' @description Reads 384-well plate metadata from Excel files with specific
#'   format (sheets containing data in B2:Y17 range). Supports both wide and
#'   long format output.
#'
#' @param file Path to the Excel file.
#' @param format Output format: "wide" (one row per well, columns per sheet)
#'   or "long" (one row per well-sheet combination).
#' @param add_directory Logical. If TRUE, adds a directory column extracted
#'   from the file path.
#' @param colname_prefix String prefix to add to all column names.
#'
#' @return A tibble containing plate metadata. Returns NULL if file doesn't exist.
#'   For wide format: columns are well, sheet names, and optionally directory.
#'   For long format: columns are Metadata_sheet, Metadata_well, and value.
#'
#' @examples
#' \dontrun{
#' read_metadata("plate_info.xlsx", format = "wide")
#' read_metadata("plate_info.xlsx", format = "long", add_directory = TRUE)
#' }
#'
#' @export
read_metadata <- function(file, format = c("wide", "long"), add_directory = FALSE, colname_prefix = "") {
  format <- match.arg(format)
  if (!file.exists(file)) return(NULL)
  print(str_glue("reading: {file}"))
  
  sheets <- readxl::excel_sheets(file) %>% set_names(., nm = .)
  print(str_glue("expected sheets: {str_c(sheets, collapse=', ')}"))
  
  read_sheet <- function(sheet) {
    readxl::read_xlsx(file, sheet = sheet, range = "B2:Y17",
                      col_names = FALSE, col_types = "text",
                      .name_repair = "minimal") %>%
      as.matrix() %>% as.vector()
  }
  
  wells <- map2_chr(rep(LETTERS[1:16], 24), rep(1:24, each = 16), ~ str_c(.x, .y))
  
  if (format == "wide") {
    plate_info <- sheets %>%
      map(~ enframe(read_sheet(.x), name = NULL, value = .x)) %>%
      bind_cols() %>%
      bind_cols(well = wells, .) %>%
      filter(!if_all(!well, is.na))
    
    if (add_directory) {
      plate_info <- plate_info %>%
        mutate(directory = file %>% str_split("[/]", simplify = TRUE) %>%
                 str_subset(".*-Measurement [0-9]") %>%
                 str_split("(__)", simplify = TRUE) %>% .[1], .before = 1)
    }
    
    if (nchar(colname_prefix) > 0) {
      plate_info <- plate_info %>%
        rename_with(~ str_c(colname_prefix, .x), everything())
    }
  } else {
    plate_info <- sheets %>%
      map(~ enframe(read_sheet(.x), name = NULL, value = "value") %>%
            bind_cols(Metadata_well = wells, .) %>%
            filter(!is.na(value))) %>%
      bind_rows(.id = "Metadata_sheet")
  }
  
  return(plate_info)
}

# =============================================================================
# Pattern Matching
# =============================================================================

#' Find Best Match from Pattern List
#'
#' @description Searches through a vector of choices using a list of patterns,
#'   returning the first match found. Useful for auto-selecting reference values
#'   from user data.
#'
#' @param choices Character vector of available choices.
#' @param patterns Character vector of patterns to search for (in order of priority).
#'
#' @return The first choice matching any pattern, or the first choice if no
#'   matches are found. Returns empty string if choices is empty.
#'
#' @examples
#' find_best_match(c("DMSO", "PBS", "Water"), c("blank", "dmso"))
#' # Returns: "DMSO"
#'
#' find_best_match(c("Sample1", "Sample2"), c("blank", "control"))
#' # Returns: "Sample1" (no match, returns first choice)
#'
#' @export
find_best_match <- function(choices, patterns) {
  if (length(choices) == 0) return("")
  
  choices_lower <- tolower(choices)
  for (pattern in patterns) {
    matches <- choices[grepl(pattern, choices_lower, ignore.case = TRUE)]
    if (length(matches) > 0) {
      return(matches[1])
    }
  }
  return(choices[1])
}

# =============================================================================
# Gene ID Conversion
# =============================================================================

#' Convert Gene IDs Between Different Database Identifiers
#'
#' @description Converts gene identifiers between SYMBOL, ENTREZID, and ENSEMBL
#'   formats using Bioconductor annotation packages (org.Hs.eg.db for human,
#'   org.Mm.eg.db for mouse).
#'
#' @param genelist Character vector of gene identifiers to convert.
#' @param fromType Character. The type of input identifiers: "SYMBOL",
#'   "ENSEMBL", or "ENTREZID".
#' @param species Character. Species code: "hs" for human, "mm" for mouse.
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item SYMBOL - Gene symbol
#'     \item ENTREZID - Entrez Gene ID
#'     \item ENSEMBL - Ensembl Gene ID
#'     \item query - Original input identifier
#'   }
#'
#' @examples
#' \dontrun{
#' id_convert(c("TP53", "BRCA1"), fromType = "SYMBOL", species = "hs")
#' }
#'
#' @export
id_convert <- function(genelist, fromType = "SYMBOL", species = "hs") {
  if (!is.character(genelist) || length(genelist) == 0) {
    stop("genelist must be a non-empty character vector")
  }
  
  orgDb <- if (species == "hs") org.Hs.eg.db::org.Hs.eg.db
           else org.Mm.eg.db::org.Mm.eg.db
  
  if (fromType == "SYMBOL") {
    genelist <- if (species == "hs") toupper(genelist)
                else tools::toTitleCase(tolower(genelist))
  }
  
  result <- tryCatch({
    AnnotationDbi::select(orgDb, keys = genelist, keytype = fromType,
                          columns = c("SYMBOL", "ENTREZID", "ENSEMBL"))
  }, error = function(e) {
    data.frame(SYMBOL = character(0), ENTREZID = character(0),
               ENSEMBL = character(0))
  })
  
  result$query <- result[[fromType]]
  result
}
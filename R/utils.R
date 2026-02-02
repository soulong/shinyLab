# shinyBioTools - Utility Functions

# Gene ID conversion using org.Hs.eg.db / org.Mm.eg.db
id_convert <- function(genelist, fromType = "SYMBOL", species = "hs") {
  if (!is.character(genelist) || length(genelist) == 0) {
    stop("genelist must be a non-empty character vector")
  }
  
  orgDb <- if (species == "hs") org.Hs.eg.db::org.Hs.eg.db
           else org.Mm.eg.db::org.Mm.eg.db
  
  # Normalize gene symbols
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
  
  # Add query column for tracking
  result$query <- result[[fromType]]
  result
}

# shRNA primer design for miR-30 based vectors
shRNA_primer <- function(splashRNA_result) {
  loop <- DNAString("TAGTGAAGCCACAGATGTA")
  ovlp5 <- DNAString("TGCTGTTGACAGTGAGCG")
  ovlp3 <- DNAString("TGCCTACTGCCTCGGACT")
  
  primers <- lapply(seq_len(nrow(splashRNA_result)), function(i) {
    anti <- DNAString(splashRNA_result$Antisense[i])
    mm <- DNAString(switch(as.character(anti[22]),
                           "A" = "C", "G" = "A", "C" = "A", "T" = "C"))
    sense <- c(mm, reverseComplement(anti[-22]))
    
    primer_F <- as.character(paste0(loop, anti, ovlp3))
    primer_R <- as.character(reverseComplement(
      DNAString(paste0(ovlp5, sense, loop))))
    
    id <- splashRNA_result$ID[i]
    data.frame(ID = c(paste0(id, "-F"), paste0(id, "-R")),
               Primer = c(primer_F, primer_R), stringsAsFactors = FALSE)
  })
  do.call(rbind, primers)
}

# Query splashRNA database
splashRNA <- function(id, n = 3, anno = "Gene") {
  n <- max(1, min(6, n))
  
  resp <- POST(
    url = "http://splashrna.mskcc.org/show_results",
    body = list(fasta = paste0("> entrezID\n", id),
                n_predictions = as.character(n), removeRE = "on",
                apa = "on", academic = "on", email = "user@example.com"),
    add_headers(Origin = "http://splashrna.mskcc.org",
                Referer = "http://splashrna.mskcc.org/"),
    encode = "form"
  )
  
  html <- read_html(resp)
  base_xpath <- "/html/body/div[2]/div[2]/div/table/tbody/tr["
  
  results <- lapply(seq_len(n), function(i) {
    data.frame(
      ID = paste0(anno, "#", i),
      Antisense = html %>%
        html_node(xpath = paste0(base_xpath, i, "]/td[2]")) %>%
        html_text(trim = TRUE),
      Score = html %>%
        html_node(xpath = paste0(base_xpath, i, "]/td[3]")) %>%
        html_text(trim = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}

# Batch query splashRNA
splashRNA_batch <- function(IDs, N = 3, Anno = NULL) {
  if (is.null(Anno)) Anno <- paste0("Gene_", seq_along(IDs))
  
  results <- lapply(seq_along(IDs), function(i) {
    message("Processing gene ", i, "/", length(IDs), ": ", IDs[i])
    Sys.sleep(2)
    tryCatch(
      splashRNA(IDs[i], N, Anno[i]),
      error = function(e) {
        warning("Failed: ", IDs[i], " - ", e$message)
        data.frame(ID = character(0), Antisense = character(0),
                   Score = character(0), stringsAsFactors = FALSE)
      }
    )
  })
  do.call(rbind, results)
}

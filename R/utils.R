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

# Real-time PCR analysis (delta-delta-CT method)
rtPCR <- function(data, choose_samples, choose_targets, na.do = TRUE,
                  real_file, skip_row = 42) {
  if (is.null(data) || nrow(data) == 0) return(NULL)
  
  data_filtered <- data %>%
    filter(Sample %in% choose_samples, Target %in% choose_targets) %>%
    mutate(CT = as.numeric(CT))
  
  na_samples <- na_targets <- character(0)
  
  if (isTRUE(na.do)) {
    res <- fill_na_values(data_filtered)
    data_filtered <- res$data
    na_samples <- res$na_samples
    na_targets <- res$na_targets
  } else {
    data_filtered <- filter(data_filtered, !is.na(CT))
  }
  
  ref_target <- choose_targets[1]
  plot_data <- calculate_ddct(data_filtered, choose_samples,
                               choose_targets[-1], ref_target)
  melt_data <- read_melt_curve(real_file, skip_row, data_filtered,
                                choose_samples, choose_targets)
  
  list(plot_data, na_samples, na_targets, melt_data)
}

# Fill NA values based on replicate data
fill_na_values <- function(data) {
  na_samples <- na_targets <- character(0)
  
  stats <- data %>%
    reframe(n_total = n(), n_valid = sum(!is.na(ct)),
            mean_ct = mean(ct, na.rm = TRUE), 
            .by=c(sample, target))
  fillable <- filter(stats, n_valid >= 1, n_valid < n_total)
  if (nrow(fillable) > 0) {
    na_samples <- unique(fillable$sample)
    na_targets <- unique(fillable$target)
    for (i in seq_len(nrow(fillable))) {
      r <- fillable[i, ]
      mask <- data$sample == r$sample & data$target == r$target & is.na(data$ct)
      data$ct[mask] <- r$mean_ct
    }
  }
  list(data = data, na_samples = na_samples, na_targets = na_targets)
}

# Calculate delta-delta CT values
calculate_ddct <- function(data, samples, targets, ref_target) {
  results <- lapply(targets, function(tgt) {
    delta_ct <- sapply(samples, function(smp) {
      ref_ct <- mean(data$CT[data$Target == ref_target & data$Sample == smp],
                     na.rm = TRUE)
      tgt_ct <- mean(data$CT[data$Target == tgt & data$Sample == smp],
                     na.rm = TRUE)
      tgt_ct - ref_ct
    })
    ref_mean <- mean(delta_ct, na.rm = TRUE)
    values <- 2^(-(delta_ct - ref_mean))
    data.frame(Samples = samples, Targets = tgt, Values = unname(values),
               SDs = sd(values, na.rm = TRUE), stringsAsFactors = FALSE)
  })
  
  plot_data <- do.call(rbind, results)
  plot_data$Samples <- factor(plot_data$Samples, levels = unique(samples))
  plot_data
}

# Read melt curve data from Excel file
read_melt_curve <- function(real_file, skip_row, data_filtered,
                            samples, targets) {
  tryCatch({
    melt <- read_excel(real_file, sheet = "Melt Curve Raw Data",
                       col_names = TRUE, skip = skip_row) %>%
      mutate(Well = as.character(Well), Target = NA_character_,
             Sample = NA_character_)
    
    for (s in samples) {
      wells <- data_filtered$Well[data_filtered$Sample == s]
      melt$Sample[melt$Well %in% wells] <- s
    }
    for (t in targets) {
      wells <- data_filtered$Well[data_filtered$Target == t]
      melt$Target[melt$Well %in% wells] <- t
    }
    filter(melt, !is.na(Target), !is.na(Sample))
  }, error = function(e) NULL)
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

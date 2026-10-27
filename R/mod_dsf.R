# shinyLab - DSF Module
# Differential Scanning Fluorimetry analysis module
#
# Layout: left column holds upload / reference-overlay / download controls and
# the plot settings; the right side has two tabs - "Data" (result tables) and
# "Plots". The Plots tab starts with well filters (pick any annotation columns
# and the values to keep; empty value selection = no restriction, columns
# combine by AND) and below shows ALL matching groups at once - either merged
# into one plot or as one row per group (fluorescence left, dF/dT right), as
# chosen by the left-sidebar Layout control. Curves can be colored by any
# annotation column; the grey reference overlay optionally matches the
# reference wells by target, so each target's row shows its own reference
# ligand (targets have different Tm baselines).
#
# Two upload flavors are accepted:
#   - the exported .txt tables ('*AnalysisResults.txt' + '*RawData_<plate>.eds.txt')
#   - raw .eds experiment files (ZIP containers with JSON payloads); the
#     per-well analysis results and melt curves are read straight from
#     primary/analysis_result.json inside the container.

library(patchwork)
library(mgcv)
library(gratia)
library(DT)
library(zip)
library(janitor)
# jsonlite is used via explicit namespace (jsonlite::fromJSON) only - attaching
# it would mask shiny::validate and break validate(need(...)) in other modules.

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

# Parse a raw QuantStudio .eds experiment (ZIP container with JSON payloads)
# into the same shapes the .txt exports produce:
#   res: one row per well x target  (well, plate, sample_name, target, ligand, tm)
#   mc : per-well melt curves       (plate, well, reading, temperature,
#                                    fluorescence, derivative)
# QuantStudio stores the melt derivative as -dF/dT; it is negated here so every
# curve in the app is plain dF/dT with the transition as a positive peak. The
# instrument's own peak calls (meltResult$peaks, ordered by confidence) often
# sit on minor tail wiggles in protein-melt runs, so Tm is called at the most
# prominent peak/trough of the derivative instead (ungated here - the Tm gate
# applies when the module re-calls it per the UI control), falling back to the
# instrument value (peak 1, percentage == 100) only when no derivative is
# available.
eds_experiment_from_dir <- function(td, eds_name) {
  get_json <- function(rel) jsonlite::fromJSON(file.path(td, rel), simplifyVector = TRUE)

  smry <- tryCatch(get_json("summary.json"), error = function(e) NULL)
  plate_name <- if (!is.null(smry) && !is.null(smry$name) && nzchar(smry$name)) {
    smry$name
  } else {
    tools::file_path_sans_ext(basename(eds_name))
  }

  ar <- tryCatch(get_json("primary/analysis_result.json"), error = function(e) NULL)
  if (is.null(ar) || is.null(ar$wellResults) || nrow(ar$wellResults) == 0) {
    stop("no well analysis results found (expected 'primary/analysis_result.json' inside the container)")
  }
  wr <- ar$wellResults

  block_type <- if (!is.null(smry$blockType)) smry$blockType else "BLOCK_384W"
  n_total <- suppressWarnings(as.integer(sub("BLOCK_(\\d+)W.*", "\\1", block_type)))
  if (is.na(n_total) || n_total <= 0) n_total <- 384L
  # rows per block: 96-well plates use 8 rows (A-H), 384-well use 16 (A-P)
  n_rows <- if (n_total <= 96L) 8L else 16L
  n_cols <- max(1L, n_total %/% n_rows)
  well_pos <- paste0(LETTERS[(wr$wellIndex %/% n_cols) + 1L],
                     (wr$wellIndex %% n_cols) + 1L)

  res_rows <- vector("list", nrow(wr))
  mc_rows <- vector("list", nrow(wr))
  for (i in seq_len(nrow(wr))) {
    w <- well_pos[[i]]
    rr <- wr$reactionResults[[i]]
    if (is.null(rr) || nrow(rr) == 0) {
      res_rows[[i]] <- data.frame(well = w, sample_name = wr$sampleName[[i]],
                                  target = w, tm = NA_real_)
      next
    }
    res_rows[[i]] <- data.frame(
      well = w,
      sample_name = wr$sampleName[[i]],
      target = ifelse(is.na(rr$targetName), w, rr$targetName),
      tm = vapply(seq_len(nrow(rr)), function(r) {
        mt <- rr$meltResult$meltTemperatures[r][[1]]
        if (is.null(mt) || !length(mt)) NA_real_ else suppressWarnings(as.numeric(mt[1]))
      }, numeric(1))
    )
    mr <- rr$meltResult
    if (is.null(mr) || nrow(mr) == 0) next
    rnT <- mr$rnTemperatures[1][[1]]
    rn <- mr$rn[1][[1]]
    if (is.null(rnT) || length(rnT) < 3) next
    dT <- mr$derivativeRnTemperatures[1][[1]]
    dRn <- mr$derivativeRn[1][[1]]
    # negate: the instrument stores -dF/dT, the app works with plain dF/dT
    deriv <- if (!is.null(dT) && length(dT) > 1 && !is.null(dRn) && length(dRn) == length(dT)) {
      -approx(dT, dRn, xout = rnT, rule = 2)$y
    } else {
      rep(NA_real_, length(rnT))
    }
    # call Tm at the most prominent peak/trough of the derivative (the
    # visible main transition; range edges are trimmed inside
    # refine_extremum). Ungated here - the Tm gate control applies when the
    # module re-calls it for display; the instrument's meltT often misses
    # the main peak in protein-melt runs
    if (any(is.finite(deriv))) {
      res_rows[[i]]$tm <- refine_extremum(rnT, deriv, gate = 0)
    }
    mc_rows[[i]] <- data.frame(
      plate = plate_name, well = w,
      reading = seq_along(rnT),
      temperature = suppressWarnings(as.numeric(rnT)),
      fluorescence = suppressWarnings(as.numeric(rn)),
      derivative = deriv
    )
  }

  res <- list_rbind(res_rows)
  res$plate <- plate_name
  # the sample name doubles as the ligand annotation for typical DSF layouts
  res$ligand <- res$sample_name
  res <- res %>% relocate(tm, .after = last_col())
  mc <- list_rbind(mc_rows)
  list(res = res, mc = mc, plate = plate_name)
}

# Open one uploaded .eds file. Returns:
#   kind = "json" : parsed experiment ($res / $mc / $plate)
#   kind = "txt"  : the container holds the exported .txt tables
#                   ($paths / $names, to be fed into the .txt pipeline)
#   kind = "none" : nothing recognisable inside
expand_eds_file <- function(path, name) {
  td <- tempfile(pattern = "eds_json_")
  dir.create(td)
  unzipped <- tryCatch({
    zip::unzip(path, exdir = td)
    list.files(td, recursive = TRUE, full.names = TRUE)
  }, error = function(e) character())
  if (!length(unzipped)) {
    return(list(kind = "none", name = name,
                note = "could not be opened as a ZIP-based EDS container"))
  }
  # JSON payload (raw instrument experiment)
  if (any(basename(unzipped) == "analysis_result.json")) {
    exp <- tryCatch(eds_experiment_from_dir(td, name), error = function(e) {
      showNotification(paste0("EDS '", name, "' could not be parsed: ", conditionMessage(e)),
                       type = "error", duration = 10)
      NULL
    })
    if (!is.null(exp)) return(list(kind = "json", exp = exp))
  }
  # exported .txt tables packed into a container
  res_f <- unzipped[grepl("AnalysisResults", basename(unzipped))]
  mc_f <- unzipped[grepl("RawData", basename(unzipped))]
  if (length(res_f) || length(mc_f)) {
    return(list(kind = "txt",
                paths = c(res_f, mc_f),
                names = basename(c(res_f, mc_f))))
  }
  list(kind = "none", name = name,
       note = "no analysis_result.json and no exported tables inside the container")
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

# Per-well GAM derivative of the fluorescence trace; one smooth per well.
# `k` caps the smooth's basis dimension (higher = more detail, default 10);
# the derivative is evaluated on a grid fine enough that the parabolic Tm
# refinement below is never coarser than `grid_step` degrees.
compute_derivatives <- function(mc_tidy, k = 10, grid_step = 0.1) {
  if (nrow(mc_tidy) == 0) {
    return(tibble(plate = character(), well = character(),
                  temperature = numeric(), derivative = numeric()))
  }
  k <- max(6L, suppressWarnings(as.integer(k)))
  if (is.na(k)) k <- 10L
  n_grid <- function(tmin, tmax) {
    max(100L, as.integer(ceiling((tmax - tmin) / grid_step)) + 1L)
  }
  mc_tidy %>%
    filter(!is.na(temperature), !is.na(fluorescence)) %>%
    nest(.by = c(plate, well)) %>%
    mutate(pred = purrr::map(data, \(d) {
      tryCatch({
        mod <- mgcv::gam(fluorescence ~ s(temperature, k = k), data = d)
        n <- n_grid(min(d$temperature), max(d$temperature))
        res <- tryCatch({
          r <- gratia::derivatives(mod, n = n, order = 1)
          dcol <- if (".derivative" %in% names(r)) ".derivative" else
                  if ("derivative" %in% names(r)) "derivative" else
                    stop("unexpected output columns: ",
                         paste(names(r), collapse = ", "), call. = FALSE)
          tibble(temperature = as.numeric(r$temperature),
                 derivative = as.numeric(r[[dcol]]))
        }, error = function(e) {
          # gratia missing or incompatible with this mgcv version: fall back
          # to central differences of the GAM predictions on the same grid
          g <- seq(min(d$temperature), max(d$temperature), length.out = n)
          pr <- as.numeric(stats::predict(mod, newdata = data.frame(temperature = g)))
          tibble(temperature = (g[-1] + g[-length(g)]) / 2,
                 derivative = diff(pr) / diff(g))
        })
        filter(res, !is.na(derivative))
      }, error = function(e) NULL)
    })) %>%
    unnest(pred, keep_empty = TRUE) %>%
    dplyr::select(!any_of(c("data", "pred"))) %>%
    { if (!"derivative" %in% colnames(.)) mutate(., derivative = NA_real_) else . } %>%
    filter(!is.na(derivative))
}

# Sub-grid Tm from the derivative: the main transition is the LOCAL EXTREMUM
# OF LARGEST PROMINENCE of dF/dT - the peak or trough that rises furthest
# above the surrounding curve. Prominence (not magnitude) is what makes this
# robust: on a rising melt the post-melt decay can be steeper than the
# transition itself, but its wiggles hug the falling baseline, so the true
# peak still wins. Quenching-style falling curves are covered automatically,
# since their transition is the prominent trough. A winner must also stand
# out from the ambient |dF/dT| level (at least `gate` times the median; a
# gate of 0 or NA disables the check): a steady rise or decay with no
# cooperative transition has no such feature and reports no Tm.
# The first/last `trim` degrees are excluded, since both the instrument
# derivative and the GAM fits carry large unphysical spikes at range edges.
refine_extremum <- function(t, d, trim = 5, gate = 1.2) {
  t <- as.numeric(t); d <- as.numeric(d)
  if (!length(d) || all(is.na(d))) return(NA_real_)
  keep <- which(!is.na(t) & !is.na(d))
  if (!length(keep)) return(NA_real_)
  t <- t[keep]; d <- d[keep]
  if (is.finite(trim) && trim > 0) {
    in_rng <- t > min(t) + trim & t < max(t) - trim
    if (sum(in_rng) >= 3L) { t <- t[in_rng]; d <- d[in_rng] }
  }
  n <- length(d)
  if (n < 3L) return(NA_real_)
  # interior local extrema; strict on the left so a tied run counts once
  peaks <- which(d[2:(n - 1L)] > d[1:(n - 2L)] & d[2:(n - 1L)] >= d[3:n]) + 1L
  troughs <- which(d[2:(n - 1L)] < d[1:(n - 2L)] & d[2:(n - 1L)] <= d[3:n]) + 1L
  if (!length(peaks) && !length(troughs)) return(NA_real_)
  # prominence of one extremum: its height above the higher of the two
  # saddles crossed when walking away from it (up to higher terrain or the
  # range edge)
  prom <- function(x, i) {
    h <- x[i]; lmin <- h
    if (i > 1L) for (j in (i - 1L):1L) { if (x[j] > h) break; lmin <- min(lmin, x[j]) }
    rmin <- h
    if (i < n)  for (j in (i + 1L):n)  { if (x[j] > h) break; rmin <- min(rmin, x[j]) }
    h - max(lmin, rmin)
  }
  cand <- sort(union(peaks, troughs))
  p <- vapply(cand, function(i) if (i %in% peaks) prom(d, i) else prom(-d, i),
              numeric(1))
  w <- cand[which.max(p)]
  # a real transition is steeper than the ambient slope of the trace; on a
  # steady rise/decay the most prominent wiggle merely rides on that slope,
  # and the honest result is no Tm (gate 0 / NA disables the check)
  if (is.finite(gate) && gate > 0 && abs(d[w]) < gate * median(abs(d)))
    return(NA_real_)
  if (w <= 1L || w >= n) return(t[w])
  d0 <- d[w - 1L]; d1 <- d[w]; d2 <- d[w + 1L]
  den <- d0 - 2 * d1 + d2
  delta <- if (!is.finite(den) || den == 0) 0 else 0.5 * (d0 - d2) / den
  delta <- min(1, max(-1, delta))
  t[w] + delta * (t[w + 1L] - t[w])
}

# Smooth a fluorescence curve table for DISPLAY: a loess fit per well,
# evaluated on a thinner temperature grid (~n_out points). Only the
# temperature/fluorescence columns are touched - dF/dT curves and Tm calls
# keep using the raw data.
smooth_display <- function(tbl, span = 0.2, n_out = 120) {
  keep_cols <- setdiff(colnames(tbl),
                       c("plate", "well", "reading", "temperature",
                         "fluorescence", "derivative"))
  annot <- tbl %>%
    group_by(plate, well) %>%
    summarise(across(all_of(keep_cols), ~ dplyr::first(.x)), .groups = "drop")
  tbl %>%
    dplyr::select(plate, well, temperature, fluorescence) %>%
    filter(!is.na(temperature), !is.na(fluorescence)) %>%
    nest(.by = c(plate, well)) %>%
    mutate(sm = purrr::map(data, function(d) {
      tryCatch({
        g <- seq(min(d$temperature), max(d$temperature),
                 length.out = max(20L, min(n_out, nrow(d))))
        fit <- stats::loess(fluorescence ~ temperature, data = d, span = span,
                            degree = 2,
                            control = stats::loess.control(surface = "direct"))
        tibble(temperature = g,
               fluorescence = as.numeric(predict(fit, data.frame(temperature = g))))
      }, error = function(e) d)
    })) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(cols = sm) %>%
    inner_join(annot, by = c("plate", "well"))
}

# Thin a curve table for on-screen display: per well, split the temperature
# range into n_bins buckets and keep the rows holding each bucket's min and
# max display value. Peaks and transitions survive while the point count
# drops ~3x, which makes the plots much cheaper to render.
thin_display <- function(tbl, value = "fluorescence", n_bins = 200) {
  if (!all(c("plate", "well", "temperature") %in% colnames(tbl)) ||
      !value %in% colnames(tbl)) return(tbl)
  tbl %>%
    filter(!is.na(temperature), !is.na(.data[[value]])) %>%
    group_by(plate, well) %>%
    mutate(.bin = ceiling((temperature - min(temperature)) /
                            pmax(max(temperature) - min(temperature), 1e-9) * n_bins)) %>%
    group_by(plate, well, .bin) %>%
    filter(.data[[value]] %in% range(.data[[value]])) %>%
    ungroup() %>%
    dplyr::select(-.bin)
}

# Combined key over one or more columns ("PK1 | DMSO"). With collapse_na the
# key becomes NA when any component is NA (used for grouping, so NA wells are
# not plotted); otherwise NA renders as "NA" (used for coloring).
paste_key <- function(df, cols, collapse_na = TRUE) {
  if (length(cols) == 0) return(rep(NA_character_, nrow(df)))
  cols <- cols[cols %in% colnames(df)]
  if (length(cols) == 0) return(rep(NA_character_, nrow(df)))
  vals <- as.data.frame(df[, cols, drop = FALSE])
  vals[] <- lapply(vals, as.character)
  out <- do.call(paste, c(vals, sep = " | "))
  if (collapse_na) out[rowSums(is.na(vals)) > 0] <- NA_character_
  out
}

# All rows of a curve table whose combined key over `cols` equals `g`. `key`
# may carry a precomputed paste_key(tbl, cols) so callers drawing many groups
# reuse one full-table pass instead of re-scanning the table per plot.
group_curve_df <- function(tbl, cols, g, key = NULL) {
  if (length(cols) == 0) return(tbl[0, , drop = FALSE])
  k <- if (is.null(key)) paste_key(tbl, cols) else key
  tbl[!is.na(k) & k %in% g, , drop = FALSE]
}

# Reference-group overlay rows, split in two: ref_rows() does the (expensive,
# group-independent) row filter and runs once per data/reference change, while
# the cheap per-group self check lives in ref_is_self(). An overlay is drawn
# only when ref_rows() returns data and ref_is_self() is FALSE.
ref_rows <- function(tbl, ref) {
  if (is.null(ref) || !length(ref$cols)) return(NULL)
  idx <- rep(TRUE, nrow(tbl))
  for (j in seq_along(ref$cols)) {
    if (!ref$cols[[j]] %in% colnames(tbl)) return(NULL)
    idx <- idx & !is.na(tbl[[ref$cols[[j]]]]) & tbl[[ref$cols[[j]]]] == ref$vals[[j]]
  }
  out <- tbl[idx, , drop = FALSE]
  if (nrow(out) == 0) NULL else out
}

ref_is_self <- function(ref, group_cols, g) {
  if (is.null(ref) || !length(ref$cols)) return(FALSE)
  m <- match(ref$cols, group_cols)
  if (anyNA(m)) return(FALSE)
  parts <- strsplit(g, " | ", fixed = TRUE)[[1]]
  length(parts) == length(group_cols) && all(parts[m] == ref$vals)
}

# Prepare the color mapping: a single numeric column colors continuously,
# anything else (single discrete column or a combination) colors discretely
# via a combined key.
color_prep <- function(df, cols) {
  cols <- cols[cols %in% colnames(df)]
  if (length(cols) == 0) return(list(df = df, col = NULL, name = ""))
  if (length(cols) == 1 && is.numeric(df[[cols]])) {
    list(df = df, col = cols, name = cols)
  } else {
    df$.color_key <- paste_key(df, cols, collapse_na = FALSE)
    list(df = df, col = ".color_key", name = paste(cols, collapse = " \u00d7 "))
  }
}

# Build the melt-curve ggplot used for both the on-screen view and the
# vector-PDF export. `data` must contain temperature, well and the y_var
# column ("fluorescence" or "derivative").
build_curve_ggplot <- function(data, y_var = c("fluorescence", "derivative"),
                               color_by = "", color_name = "", temp_range = NULL,
                               ref_data = NULL, tm_lines = NULL, ref_label = NULL) {
  y_var <- match.arg(y_var)
  y_lab <- if (y_var == "derivative") "dF/dT" else "Fluorescence"

  if ((is.null(data) || nrow(data) == 0) && (is.null(ref_data) || nrow(ref_data) == 0)) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, label = "No data for this group", size = 4) +
             theme_void())
  }

  has_color <- !is.null(color_by) && !identical(color_by, "") &&
               !is.null(data) && nrow(data) > 0 && color_by %in% colnames(data)
  color_discrete <- has_color && !is.numeric(data[[color_by]])
  has_ref <- !is.null(ref_data) && nrow(ref_data) > 0
  if (has_ref && is.null(ref_label)) ref_label <- "ref"

  # ONE explicit colour scale shared by curves, the grey reference, the Tm
  # marker lines and their labels: ligand groups take the magma palette, the
  # reference takes grey60 under its own "<value> (ref)" legend entry.
  scale_vals <- NULL
  if (color_discrete) {
    lv <- as.character(unique(data[[color_by]]))
    scale_vals <- setNames(viridisLite::magma(length(lv)), lv)
    data[[color_by]] <- as.character(data[[color_by]])
  }
  if (has_ref && !is.null(ref_label)) {
    if (color_discrete) {
      ref_data[[color_by]] <- ref_label
      scale_vals <- c(scale_vals, setNames("grey60", ref_label))
    } else if (!has_color) {
      ref_data$.ref_key <- ref_label
      scale_vals <- c(scale_vals, setNames("grey60", ref_label))
    }
  }

  p <- ggplot()

  if (has_ref) {
    ref_aes <- if (color_discrete) {
      aes(temperature, .data[[y_var]], group = well, colour = .data[[color_by]])
    } else if (!has_color && !is.null(ref_label)) {
      aes(temperature, .data[[y_var]], group = well, colour = .data[[".ref_key"]])
    } else {
      NULL
    }
    p <- p + if (!is.null(ref_aes)) {
      geom_path(data = ref_data, mapping = ref_aes, linewidth = 0.4, na.rm = TRUE)
    } else {
      geom_path(data = ref_data, mapping = aes(temperature, .data[[y_var]], group = well),
                colour = "grey60", linewidth = 0.4, inherit.aes = FALSE)
    }
  }

  if (!is.null(data) && nrow(data) > 0) {
    main_aes <- if (has_color) {
      aes(temperature, .data[[y_var]], group = well, colour = .data[[color_by]])
    } else {
      aes(temperature, .data[[y_var]], group = well)
    }

    p <- p + if (has_color) {
      geom_path(data = data, mapping = main_aes, linewidth = 0.6, na.rm = TRUE)
    } else {
      geom_path(data = data, mapping = main_aes, colour = "darkred",
                linewidth = 0.6, na.rm = TRUE)
    }

    if (color_discrete) {
      # legend keys: all subgroups, or (when there are too many to read, e.g.
      # colouring by well) just the grey reference entry
      brks <- if (length(lv) > 20L) ref_label else c(lv, ref_label)
      p <- p + scale_colour_manual(values = scale_vals, breaks = brks,
                                   limits = names(scale_vals),
                                   name = color_name, na.value = "darkred")
    } else if (has_color) {
      p <- p + scale_colour_viridis_c(option = "C", na.value = "darkred", name = color_name)
    } else if (!is.null(ref_label)) {
      p <- p + scale_colour_manual(values = scale_vals, breaks = ref_label,
                                   name = "", na.value = "grey30")
    }
  }

  # Median-Tm markers for the dF/dT view: one dashed line per colour subgroup
  # (and the reference), each carrying its subgroup colour, with a horizontal
  # value label placed just to the right of the line. Labels alternate between
  # two heights so neighbours do not fully overlap.
  if (!is.null(tm_lines) && nrow(tm_lines) > 0) {
    tm_colored <- color_discrete
    tm_lines$.vj <- ifelse(seq_len(nrow(tm_lines)) %% 2 == 1, 1.2, 2.8)
    tm_args <- list(linetype = "dashed", linewidth = 0.4, show.legend = FALSE)
    txt_args <- list(y = Inf, label = sprintf("%.1f", tm_lines$tm),
                     hjust = -0.08, size = 3, show.legend = FALSE)
    if (tm_colored) {
      p <- p + do.call(geom_vline, c(list(data = tm_lines,
                                          mapping = aes(xintercept = tm,
                                                        colour = .data[["value"]])), tm_args))
      p <- p + do.call(geom_text, c(list(data = tm_lines,
                                         mapping = aes(x = tm, colour = .data[["value"]],
                                                       vjust = .data[[".vj"]])),
                                     txt_args))
    } else {
      p <- p + do.call(geom_vline, c(list(data = tm_lines,
                                          mapping = aes(xintercept = tm),
                                          colour = "grey30"), tm_args))
      p <- p + do.call(geom_text, c(list(data = tm_lines,
                                         mapping = aes(x = tm, vjust = .data[[".vj"]]),
                                         colour = "grey30"),
                                     txt_args))
    }
  }

  p + theme_bw(base_size = 11) +
    theme(legend.position = "top") +
    labs(x = "Temperature (\u00b0C)", y = y_lab) +
    { if (!is.null(temp_range) && length(temp_range) == 2 && !anyNA(temp_range))
        coord_cartesian(xlim = as.numeric(temp_range)) else NULL }
}

# Legend label for the grey reference overlay: "<value> (ref)" when the
# reference wells carry the colour column, otherwise plain "ref".
ref_legend_label <- function(ref_df, color_col) {
  if (is.null(ref_df) || nrow(ref_df) == 0 || is.null(color_col) ||
      length(color_col) != 1 || !color_col %in% colnames(ref_df)) return("ref")
  v <- unique(ref_df[[color_col]]); v <- v[!is.na(v)]
  if (length(v) == 1) paste0(v, " (ref)") else "ref"
}

# "PK1 | DMSO  |  6 wells  |  median Tm = 52.3 C" for a row header.
group_header_label <- function(df, cols, g) {
  n_wells <- 0L
  tm_med <- NA_real_
  if (length(cols)) {
    k <- paste_key(df, cols)
    idx <- !is.na(k) & k %in% g
    n_wells <- sum(idx)
    if ("tm" %in% colnames(df)) {
      tm_med <- suppressWarnings(median(df$tm[idx], na.rm = TRUE))
    }
  }
  tm_part <- if (is.na(tm_med)) "" else sprintf("  |  median Tm = %.2f \u00b0C", tm_med)
  sprintf("%s  |  %d wells%s", g, n_wells, tm_part)
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
        fileInput(ns("file_analysis"), "QuantStudio DSF exports (.txt / .eds)",
                  accept = c(".txt", ".tsv", ".eds"), multiple = TRUE),
        helpText("Upload '*AnalysisResults.txt' + '*RawData_<plate>.eds.txt' - or the .eds files.",
                 hint_tip("Plain .txt exports are read directly. Raw .eds experiment files are parsed from the JSON payload inside the container (per-well analysis results and melt curves); EDS runs carry no Protein/Ligand properties, so the sample name is used as the 'ligand' annotation unless a metadata xlsx supplies real annotations. Plate names come from the 'Experiment File Name' column (txt) or the experiment name in summary.json (EDS); each plate's curves are read from the RawData file whose name contains that plate name.")),
        fileInput(ns("file_metadata"), "Well metadata (.xlsx, optional)", accept = ".xlsx"),
        helpText("Optional target/ligand/etc. annotations per well.",
                 hint_tip("Name the file after its raw-data file (e.g. '...RawData_<plate>.eds.xlsx' for '...RawData_<plate>.eds.txt', or '<experiment>.xlsx' for '<experiment>.eds') so annotations are assigned to the correct plate - preferred for multi-plate runs. A name matching the AnalysisResults file assigns it to the single plate; if nothing matches, annotations apply to all plates by well. Each sheet name becomes a well-annotation column ('target' and 'ligand' required; optional 'conc', 'smiles'). Fill every sheet with 384-well values in range B2:Y17. Without a metadata file, the export's own annotations are used: Protein/Ligand properties for .txt uploads, or the sample name (as ligand) and target for .eds uploads. Wells left empty in a sheet keep NA for that annotation and stay out of groups keyed on it (e.g. protein-only controls with an empty 'ligand' cell).")),
        radioButtons(ns("tm_type"), tagList("Tm Calculation",
                     hint_tip("Derivative / Boltzmann: Tm values from the uploaded files (txt exports carry both; for .eds files Tm is picked at the most prominent peak/trough of the instrument derivative - the feature that rises furthest above the surrounding curve - since the instrument's own peak calls can miss the main transition of protein-melt runs). Both derivative-based sources ignore the first/last 5 \u00b0C of the temperature range, where the traces carry edge spikes, and only call a Tm when the feature stands out from the ambient dF/dT level (threshold: the 'Tm gate' control below; traces with a steady rise or decay but no cooperative transition report no Tm). App dF/dT: Tm computed inside the app - the loaded curves are GAM-smoothed (see Curve detail below) and Tm is the most prominent peak/trough of that smoothed dF/dT; bundled Tm columns are ignored.")),
                     choices = c("Derivative" = "tm_d", "Boltzmann" = "tm_b",
                                 "App dF/dT" = "app_dfdt"),
                     selected = "app_dfdt", inline = TRUE),
        numericInput(ns("tm_gate"), tagList("Tm gate (\u00d7 ambient dF/dT)",
                     hint_tip("A Tm is only called when the prominent peak/trough exceeds the ambient |dF/dT| level (per-trace median) by this factor: higher values call only clear transitions, lower values also call weak ones, 0 disables the check (every well with any prominent feature gets a Tm). Applies to both Derivative and App dF/dT sources; changes apply immediately, no re-Analyze needed. Default 1.2.")),
                     value = 1.2, min = 0, max = 100, step = 0.1),
        conditionalPanel(
          condition = "input.tm_type === 'app_dfdt'",
          ns = ns,
          numericInput(ns("dfdt_k"), tagList("Curve detail (k)",
                       hint_tip("Upper limit on the GAM basis dimension for the app-computed dF/dT: higher k follows sharp or double transitions more closely (wigglier dF/dT), lower k over-smooths. Default 10; try 20-40 for noisy multi-peak wells. Changing it recomputes the Tm values and the dF/dT plots.")),
                       value = 10, min = 6, max = 100, step = 1)
        ),
        div(style = "display: flex; justify-content: space-between; align-items: center; gap: 4px; flex-wrap: nowrap;",
          actionButton(ns("submit"), "Analyze", icon = icon("play"), class = "btn-default btn-sm",
                       title = "Parse the uploads, annotate wells and compute curve derivatives."),
          div(style = "display: flex; gap: 4px; flex-wrap: nowrap;",
            downloadButton(ns("download_data"), "Data", icon = icon("table"), class = "btn-default btn-sm",
                           title = "Per-well results, Tm summary, metadata, raw fluorescence and derivatives as one .xlsx, one sheet each."),
            downloadButton(ns("download_plots"), "Plots", icon = icon("file-image"), class = "btn-default btn-sm",
                           title = "Per-target vector PDFs (one fixed-size page per ligand: RFU left, dF/dT right; ligand wells dark blue, reference grey; dashed lines mark the ligand Tm (labelled) and the reference Tm; the title reports Tm and delta Tm) plus a Delta Tm overview PDF, zipped. Falls back to per-group pages when target/ligand annotations are absent.")))
      ),
      box(title = "Plot Settings", width = NULL, status = "primary", solidHeader = TRUE,
          collapsible = TRUE,
        selectInput(ns("row_group"), tagList("Group curves by",
                    hint_tip("One or more columns from the per-well results. In the 'one row per group' layout every unique combination gets its own row (fluorescence left, dF/dT right); in the merged layout it defines the curve groups (used for the default colours and the PDF export). Available after Analyze.")),
                    choices = character(0), multiple = TRUE),
        selectInput(ns("color_var"), tagList("Color curves by",
                    hint_tip("One column, or a combination, to color the curves. A single numeric column (e.g. conc or Tm) uses a continuous scale.")),
                    choices = character(0), multiple = TRUE),
        checkboxInput(ns("smooth_curves"), tagList("Smooth fluorescence curves",
                      hint_tip("Loess-smooths the displayed fluorescence traces and thins them to ~120 points per well (less jagged). Display only: the dF/dT curves always use the raw data, and App dF/dT keeps smoothing with its own k.")),
                      value = TRUE),
        radioButtons(ns("facet_mode"), tagList("Layout",
                     hint_tip("Merged: all filtered curves in one plot. One row per group: every unique combination of 'Group curves by' gets its own row - fluorescence on the left, dF/dT on the right. All rows are shown at once.")),
                     choices = c("Merged (one plot)" = "merged", "One row per group" = "facet"),
                     selected = "facet", inline = TRUE),
        sliderInput(ns("temp_range"), tagList("Temperature range (\u00b0C)",
                    hint_tip("Drag the handles, then release - the plots update once on release, not continuously while dragging.")),
                    min = 20, max = 100,
                    value = c(20, 100), step = 1, dragRange = FALSE),
        strong("Reference overlay"),
        helpText("Grey curves drawn into the other groups.",
                 hint_tip("Pick one or more columns and a value for each; wells matching all values form the reference group, whose curves are overlaid in grey.")),
        selectInput(ns("ref_cols"), "Reference columns", choices = character(0), multiple = TRUE),
        uiOutput(ns("ref_values_ui")),
        checkboxInput(ns("overlay_ref"), "Overlay reference group (grey)", value = TRUE),
        checkboxInput(ns("overlay_by_target"), tagList("Match reference by target",
                      hint_tip("Each target's panel overlays its OWN reference wells (same target, e.g. its own DMSO control) instead of one fixed reference. Targets melt at different temperatures, so comparing every target against a single target's reference would be misleading.")),
                      value = TRUE)
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
            # ---- split uploads into .txt exports and raw .eds containers -----
            is_eds <- tolower(tools::file_ext(input$file_analysis$name)) == "eds"
            txt_paths <- input$file_analysis$datapath[!is_eds]
            txt_names <- input$file_analysis$name[!is_eds]
            eds_paths <- input$file_analysis$datapath[is_eds]
            eds_names <- input$file_analysis$name[is_eds]

            # ---- parse raw .eds containers -----------------------------------
            # kind "json" = instrument experiment (parsed here); kind "txt" =
            # the container holds the exported tables (fed to the .txt path).
            eds_res <- list(); eds_mc <- list()
            eds_name_to_plate <- character(0)
            if (length(eds_paths)) {
              incProgress(0.1, message = "Parsing EDS containers...")
              for (k in seq_along(eds_paths)) {
                res_k <- tryCatch(
                  expand_eds_file(eds_paths[[k]], eds_names[[k]]),
                  error = function(e) {
                    showNotification(paste0("EDS '", eds_names[[k]], "' could not be read: ",
                                            conditionMessage(e)),
                                     type = "error", duration = 10)
                    NULL
                  })
                if (is.null(res_k)) next
                if (identical(res_k$kind, "json")) {
                  eds_res[[length(eds_res) + 1L]] <- res_k$exp$res
                  eds_mc[[length(eds_mc) + 1L]] <- res_k$exp$mc
                  eds_name_to_plate[eds_names[[k]]] <- res_k$exp$plate
                } else if (identical(res_k$kind, "txt")) {
                  txt_paths <- c(txt_paths, res_k$paths)
                  txt_names <- c(txt_names, res_k$names)
                } else {
                  showNotification(paste0("EDS '", res_k$name, "' ", res_k$note,
                                          " - please upload the QuantStudio-exported .txt tables instead."),
                                   type = "warning", duration = 10)
                }
              }
              if (length(eds_res)) {
                showNotification("EDS: sample names are copied to the 'ligand' column; group or overlay by sample_name instead if you prefer.",
                                 type = "message", duration = 8)
              }
            }
            eds_res_all <- list_rbind(eds_res)
            eds_mc_all <- list_rbind(eds_mc)

            all_files <- txt_paths %>% set_names(., nm = txt_names)
            res_idx <- which(str_detect(names(all_files), "AnalysisResults"))
            if (length(res_idx) == 0 &&
                (is.null(eds_res_all) || nrow(eds_res_all) == 0)) {
              stop("No '*AnalysisResults*' file and no parseable .eds file found among the uploads (found: ",
                   paste(input$file_analysis$name, collapse = ", "), ").")
            }

            incProgress(0.05, message = "Reading analysis results...")
            txt_merged <- NULL
            mc_txt <- NULL
            meta <- NULL

            if (length(res_idx) > 0) {
              # ---- .txt exports --------------------------------------------
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

              # App dF/dT recomputes Tm from the curves, so a missing or
              # unparseable bundled Tm column is not fatal in that mode.
              app_tm <- identical(input$tm_type, "app_dfdt")
              tm_col <- tryCatch(
                resolve_tm_column(res, if (app_tm) "tm_d" else input$tm_type),
                error = function(e) {
                  if (!app_tm) stop(e)
                  showNotification("App dF/dT is selected - Tm is computed from the curves; no bundled Tm column is needed.",
                                   type = "message", duration = 8)
                  structure(NA_character_, matched = "none")
                })
              if (identical(attr(tm_col, "matched"), "fallback")) {
                showNotification(paste0("Requested Tm method column not found - falling back to column '", tm_col, "'."),
                                 type = "warning", duration = 8)
              }

              # Plate name: QuantStudio exports 'Experiment File Name'; fall back
              # to the results file name (e.g. '<plate>.eds_AnalysisResults.txt').
              # Some exports also carry a numeric 'Plate' column - drop it first
              # so the rename cannot collide into duplicated names.
              if ("experiment_file_name" %in% colnames(res)) {
                res <- res %>% dplyr::select(-any_of("plate")) %>%
                  dplyr::rename(plate = experiment_file_name)
              } else {
                fallback_plate <- res_name %>%
                  str_remove("_?AnalysisResults\\.txt$") %>% str_remove("\\.eds$")
                showNotification(paste0("AnalysisResults has no 'Experiment File Name' column - using '", fallback_plate, "' as the plate name."),
                                 type = "warning", duration = 8)
                res <- res %>% mutate(plate = fallback_plate)
              }

              res <- res %>%
                mutate(tm = if (!is.na(tm_col))
                         suppressWarnings(as.numeric(.data[[tm_col]]))
                       else NA_real_,
                       well = transform_well_style(well),
                       plate = str_remove(plate, "\\.eds$"))

              if (all(is.na(res$tm)) && !app_tm) {
                stop(paste0("Tm column '", tm_col, "' could not be parsed as numeric."))
              }

              incProgress(0.15, message = "Reading metadata...")
              # Custom metadata (if provided) is joined after the .txt and .eds
              # results are combined, so it annotates both upload flavors.
              # Without it, use the export's own Protein/Ligand annotations,
              # falling back to the well label where they are missing/empty.
              txt_merged <- res
              if (is.null(input$file_metadata)) {
                if ("protein" %in% colnames(txt_merged)) {
                  txt_merged <- txt_merged %>% mutate(target = if_else(is.na(protein) | trimws(protein) == "", well, protein))
                } else {
                  txt_merged <- txt_merged %>% mutate(target = well)
                }
                if ("ligand" %in% colnames(txt_merged)) {
                  txt_merged <- txt_merged %>% mutate(ligand = if_else(is.na(ligand) | trimws(ligand) == "", well, ligand))
                } else {
                  txt_merged <- txt_merged %>% mutate(ligand = well)
                }
                meta <- txt_merged %>%
                  dplyr::select(any_of(c("plate", "well", "target", "ligand", "conc", "smiles"))) %>%
                  distinct()
              }
            } else {
              showNotification("No '.txt' exports uploaded - reading everything from the .eds file(s).",
                               type = "message", duration = 6)
            }

            # ---- combined per-well annotations --------------------------------
            eds_res_ok <- if (!is.null(eds_res_all) && nrow(eds_res_all) > 0) eds_res_all
            res_all <- bind_rows(purrr::compact(list(txt_merged, eds_res_ok)))
            if (!nrow(res_all)) stop("No per-well results could be read.")

            if (!is.null(input$file_metadata)) {
              meta_files <- input$file_metadata$datapath %>%
                set_names(., input$file_metadata$name)
              meta <- meta_files %>%
                map(read_metadata) %>%
                list_rbind(names_to = "plate")
              if (!all(c("target", "ligand") %in% colnames(meta))) {
                stop("Custom metadata must contain 'target' and 'ligand' columns")
              }

              # Align the metadata's plate name with the data plates by matching the
              # metadata file name to the uploaded data files (.txt raw exports AND
              # raw .eds containers):
              #   - PRIMARY:   match a raw-data file (e.g. '...RawData_<plate>.eds.xlsx'
              #               for '...RawData_<plate>.eds.txt', or '<experiment>.eds')
              #               -> assign that file's plate.
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
              raw_names <- c(names(all_files)[str_detect(names(all_files), "\\.eds\\.txt$")],
                             eds_names)
              analysis_name <- if (length(res_idx) > 0) names(all_files)[res_idx] else NA_character_
              plates_all <- unique(res_all$plate)
              raw_file_to_plate <- function(raw_name) {
                if (raw_name %in% eds_names) return(eds_name_to_plate[[raw_name]])
                for (p in plates_all) {
                  if (!is.na(p) && str_detect(raw_name, fixed(p))) return(p)
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
                } else if (!is.na(analysis_name) && norm_stem(analysis_name) == ms) {
                  meta_target_plate[m] <- if (length(plates_all) == 1) plates_all else NA_character_
                } else {
                  meta_target_plate[m] <- NA_character_
                }
              }
              stem_to_plate <- setNames(meta_target_plate, norm_stem(names(meta_target_plate)))
              meta <- meta %>% mutate(
                plate = dplyr::coalesce(unname(stem_to_plate[ norm_stem(plate) ]), plate)
              )

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

              # Custom metadata provides target/ligand(/conc/smiles); drop the
              # export's own annotation columns first so the join can't collide
              # into suffixed columns (ligand.x / ligand.y).
              merged <- res_all %>%
                dplyr::select(-any_of(c("protein", "target", "ligand", "conc", "smiles"))) %>%
                inner_join(meta, by = c("plate", "well"))
              if (nrow(merged) == 0) {
                showNotification("Metadata plate name did not match the experiment plate; matched annotations by well only.",
                                 type = "warning", duration = 8)
                merged <- res_all %>%
                  dplyr::select(-any_of(c("protein", "target", "ligand", "conc", "smiles"))) %>%
                  inner_join(meta %>% dplyr::select(-any_of("plate")), by = "well")
              } else if (!setequal(unique(res_all$plate), unique(meta$plate))) {
                # Plates the metadata does not cover keep their export
                # annotations instead of being silently dropped by the join.
                missing_plates <- setdiff(unique(res_all$plate), unique(meta$plate))
                if (length(missing_plates)) {
                  addback <- res_all %>% filter(plate %in% missing_plates)
                  for (cn in c("target", "ligand")) {
                    if (cn %in% colnames(addback)) {
                      addback <- addback %>%
                        mutate(!!cn := if_else(is.na(.data[[cn]]) | trimws(.data[[cn]]) == "",
                                               well, .data[[cn]]))
                    } else {
                      addback[[cn]] <- addback$well
                    }
                  }
                  merged <- bind_rows(merged, addback)
                  showNotification(paste0("No metadata for plate(s): ",
                                          paste(missing_plates, collapse = ", "),
                                          " - export annotations kept for those wells."),
                                   type = "warning", duration = 8)
                }
              }
            } else {
              merged <- res_all
            }

            merged <- merged %>%
              { if ("conc" %in% colnames(.)) mutate(., conc = suppressWarnings(as.numeric(conc))) else . } %>%
              # every Tm column (the selected one plus any extra export columns
              # like tm_boltzmann / tm1 / tm2) is kept at 2 decimals on screen
              # and in the saved workbook
              mutate(across(matches("^tm"), ~ if (is.numeric(.x)) round(.x, 2) else .x)) %>%
              relocate(tm, .after = last_col())

            incProgress(0.2, message = "Processing melting curve...")
            match_mc_files <- c()
            for (x in unique(txt_merged$plate)) {
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

            missing_plates <- setdiff(unique(txt_merged$plate), names(match_mc_files))
            if (length(missing_plates) > 0) {
              showNotification(paste0("No raw data file found for plate(s): ",
                                      paste(missing_plates, collapse = ", "),
                                      " - those plates are skipped in the curve plots."),
                               type = "warning", duration = 10)
            }
            if (length(match_mc_files)) {
              mc_txt <- match_mc_files %>%
                map(read_fluorescence_table) %>%
                list_rbind(names_to = 'plate') %>%
                janitor::clean_names() %>%
                dplyr::select(!c(well)) %>%
                dplyr::rename(well = well_position) %>%
                mutate(well = transform_well_style(well)) %>%
                mutate(temperature = as.numeric(temperature),
                       fluorescence = as.numeric(fluorescence))
            }

            # ---- combined curve table ------------------------------------------
            eds_mc_ok <- if (!is.null(eds_mc_all) && nrow(eds_mc_all) > 0) eds_mc_all
            mc_all <- bind_rows(purrr::compact(list(mc_txt, eds_mc_ok)))
            if (!nrow(mc_all)) stop("No melt-curve data could be read from the uploads.")
            # .eds curves carry the instrument derivative; .txt curves get a
            # GAM-smoothed one below (the marker column drives the split).
            if (!"derivative" %in% colnames(mc_all)) mc_all$derivative <- NA_real_

            incProgress(0.25, detail = "Merging with raw data...")
            annot <- merged %>%
              dplyr::select(-any_of(c("reading", "temperature", "fluorescence", "derivative"))) %>%
              distinct()
            mc_tidy <- mc_all %>%
              inner_join(annot, by = c("plate", "well"))

            incProgress(0.3, detail = "Computing derivatives...")
            # .eds curves already carry the instrument derivative; .txt curves
            # get a GAM-smoothed derivative.
            derivatives <- bind_rows(
              compute_derivatives(mc_tidy %>% filter(is.na(derivative)) %>%
                                    dplyr::select(-derivative)),
              mc_tidy %>% filter(!is.na(derivative)) %>%
                dplyr::select(plate, well, temperature, derivative)
            ) %>%
              left_join(annot, by = c("plate", "well"))

            # Columns available for grouping / coloring / reference: every
            # per-well result column carried onto the curves (including plate,
            # well and Tm) except the pure curve-value columns.
            annot_vars <- setdiff(colnames(mc_tidy),
                                  c("reading", "temperature",
                                    "fluorescence", "derivative"))

            incProgress(0.35, detail = "Done!")

            showNotification(
              sprintf("Processed %d wells on %d plate(s); %d wells have curve data.",
                      nrow(merged), dplyr::n_distinct(merged$plate),
                      nrow(dplyr::distinct(mc_tidy, plate, well))),
              type = "message", duration = 8)

            list(merged = merged,
                 mc_tidy = mc_tidy,
                 meta = meta,
                 derivatives = derivatives,
                 annot_vars = annot_vars)
          })
        })
      }, error = function(e) {
        showNotification(paste("Error processing data:", e$message),
                         type = "error", duration = 10)
        return(NULL)
      })
    })

    # =========================================================================
    # Tm source switching: bundled Tm columns (txt exports / EDS instrument
    # peaks) or Tm computed inside the app from a smooth dF/dT of the loaded
    # curves. The per-well GAM fits are cached, so switching k or changing the
    # well filter only fits wells that have not been computed yet.
    # =========================================================================
    gam_cache <- reactiveValues(store = list())
    observeEvent(data(), { gam_cache$store <- list() })

    smooth_derivatives <- reactive({
      d <- data()
      req(d)
      k <- if (is.null(input$dfdt_k) || is.na(input$dfdt_k)) 10L else as.integer(input$dfdt_k)
      first_err <- NULL
      withProgress(message = "Computing smooth dF/dT...", value = 0, {
        # one self-contained chunk per well (keys kept in - compute_derivatives
        # groups by plate/well internally)
        chunks <- d$mc_tidy %>%
          filter(!is.na(temperature), !is.na(fluorescence)) %>%
          dplyr::select(plate, well, temperature, fluorescence) %>%
          group_by(plate, well) %>%
          group_split(.keep = TRUE)
        parts <- vector("list", length(chunks))
        for (i in seq_along(chunks)) {
          ck <- paste(chunks[[i]]$plate[[1]], chunks[[i]]$well[[1]], k, sep = "\r")
          if (is.null(gam_cache$store[[ck]])) {
            res <- tryCatch(compute_derivatives(chunks[[i]], k = k),
                            error = function(e) e)
            if (inherits(res, "error")) {
              if (is.null(first_err)) first_err <- conditionMessage(res)
              gam_cache$store[[ck]] <- NULL
            } else {
              gam_cache$store[[ck]] <- res
            }
            incProgress(1 / max(length(chunks), 1),
                        detail = sprintf("%s (%d/%d)", chunks[[i]]$well[[1]],
                                         i, length(chunks)))
          }
          parts[[i]] <- gam_cache$store[[ck]]
        }
        out <- list_rbind(parts)
        if (nrow(out) == 0) {
          # every fit failed (or no wells): report why instead of crashing in
          # the downstream joins with an empty, column-less table
          msg <- paste0("App dF/dT failed on every well",
                        if (is.null(first_err)) "." else paste0(": ", first_err))
          showNotification(msg, type = "error", duration = 12)
          out <- tibble(plate = character(), well = character(),
                        temperature = numeric(), derivative = numeric())
        }
        annot <- d$merged %>%
          dplyr::select(-any_of(c("reading", "temperature", "fluorescence", "derivative"))) %>%
          distinct()
        out %>% left_join(annot, by = c("plate", "well"))
      })
    })

    # dF/dT curve source for plots: the app-computed smooth derivative when
    # App dF/dT is selected; otherwise the uploaded derivative (instrument
    # derivative for .eds, GAM-smoothed for .txt uploads which carry none).
    dfdt_source <- reactive({
      d <- data()
      req(d)
      if (identical(input$tm_type, "app_dfdt")) {
        sd <- smooth_derivatives()
        req(sd, nrow(sd) > 0)
        sd
      } else {
        d$derivatives
      }
    })

    # per-well results with the Tm column swapped in per the Tm Calculation choice
    merged_tm <- reactive({
      d <- data()
      req(d)
      gate <- suppressWarnings(as.numeric(input$tm_gate))
      if (is.na(gate)) gate <- 1.2
      if (identical(input$tm_type, "app_dfdt")) {
        withProgress(message = "Calling Tm from app dF/dT...", value = 0.5, {
          sd <- smooth_derivatives()
          tm_by_well <- sd %>%
            group_by(plate, well) %>%
            summarise(tm = round(refine_extremum(temperature, derivative, gate = gate), 2),
                      .groups = "drop")
          d$merged %>%
            dplyr::select(-tm) %>%
            left_join(tm_by_well, by = c("plate", "well")) %>%
            relocate(tm, .after = last_col())
        })
      } else {
        # .eds wells carry the instrument derivative, so their Tm is re-called
        # live under the current gate (no re-Analyze needed); wells without
        # derivative data (.txt uploads) keep their uploaded Tm columns
        tm_by_well <- d$mc_tidy %>%
          filter(!is.na(derivative)) %>%
          group_by(plate, well) %>%
          summarise(tm = round(refine_extremum(temperature, derivative, gate = gate), 2),
                    .groups = "drop")
        d$merged %>%
          left_join(tm_by_well, by = c("plate", "well"), suffix = c("", ".re")) %>%
          mutate(tm = dplyr::coalesce(tm.re, tm)) %>%
          dplyr::select(-tm.re)
      }
    })

    # curve table with the matching tm column (used when grouping/coloring by tm)
    mc_tidy_tm <- reactive({
      d <- data()
      req(d)
      if (!identical(input$tm_type, "app_dfdt")) return(d$mc_tidy)
      tm_map <- merged_tm() %>% dplyr::select(plate, well, tm)
      d$mc_tidy %>%
        dplyr::select(-tm) %>%
        left_join(tm_map, by = c("plate", "well"))
    })

    # =========================================================================
    # Reference group (any columns / any combination) for the grey overlay
    # =========================================================================
    ref_def <- reactive({
      cols <- input$ref_cols
      if (length(cols) == 0) return(NULL)
      vals <- purrr::map_chr(seq_along(cols), function(i) {
        v <- input[[paste0("ref_val_", i)]]
        if (is.null(v)) "" else v
      })
      if (any(!nzchar(vals))) return(NULL)
      list(cols = cols, vals = vals)
    })

    output$ref_values_ui <- renderUI({
      cols <- input$ref_cols
      req(data(), length(cols) > 0)
      mc <- mc_tidy_tm()
      tagList(lapply(seq_along(cols), function(i) {
        vals <- sort(unique(stats::na.omit(mc[[cols[i]]])))
        if (!length(vals)) return(NULL)
        sel <- vals[1]
        if (identical(cols[[i]], "ligand") || identical(cols[[i]], "sample_name")) {
          cand <- find_best_match(vals, c("blank", "dmso", "pbs"))
          if (length(cand) == 1 && grepl("blank|dmso|pbs", tolower(cand))) sel <- cand
        }
        selectInput(ns(paste0("ref_val_", i)),
                    label = tagList("Reference:", tags$code(cols[[i]])),
                    choices = vals, selected = sel)
      }))
    })

    # =========================================================================
    # Grouping columns and reference rows, cached per change
    # =========================================================================
    group_cols_active <- reactive({
      tbl <- mc_tidy_tm()
      input$row_group[input$row_group %in% colnames(tbl)]
    })

    mc_ref_rows <- reactive({
      req(data())
      if (!isTRUE(input$overlay_ref)) return(NULL)
      ref_rows(mc_tidy_tm(), ref_def())
    })

    deri_ref_rows <- reactive({
      req(data())
      if (!isTRUE(input$overlay_ref)) return(NULL)
      ref_rows(dfdt_source(), ref_def())
    })

    # Fluorescence reference rows for display: smoothed along with the main
    # curves when the display-smoothing checkbox is on. The dF/dT reference
    # (deri_ref_rows) always stays raw.
    mc_ref_rows_display <- reactive({
      rd <- mc_ref_rows()
      if (!is.null(rd) && isTRUE(input$smooth_curves)) rd <- smooth_display(rd)
      rd
    })

    # =========================================================================
    # Well filter (top of the Plots tab): one or more annotation columns, each
    # with a value multi-select. Empty value selection = keep all values of
    # that column; the chosen columns combine by AND.
    # =========================================================================
    filter_cols_active <- reactive({
      tbl <- mc_tidy_tm()
      input$filter_cols[input$filter_cols %in% colnames(tbl)]
    })

    filtered_wells <- reactive({
      mc <- mc_tidy_tm()
      keep <- rep(TRUE, nrow(mc))
      for (cn in filter_cols_active()) {
        sel <- input[[paste0("filter_val_", cn)]]
        if (length(sel)) keep <- keep & !is.na(mc[[cn]]) & mc[[cn]] %in% sel
      }
      dplyr::distinct(mc[keep, c("plate", "well")])
    })

    output$filter_ui <- renderUI({
      req(data())
      div(style = "border-bottom: 1px solid #eee; padding-bottom: 6px; margin-bottom: 4px;",
        fluidRow(
          column(width = 4, selectInput(ns("filter_cols"), tagList("Filter by",
                      hint_tip("Pick one or more columns, then the values to keep for each. Leave a column's values empty to keep all of them; columns combine by AND. The plots below always show every group that matches the filter.")),
                      choices = data()$annot_vars, multiple = TRUE))
        ),
        uiOutput(ns("filter_vals_ui"))
      )
    })

    output$filter_vals_ui <- renderUI({
      cols <- filter_cols_active()
      req(length(cols) > 0)
      mc <- mc_tidy_tm()
      fluidRow(lapply(cols, function(cn) {
        vals <- sort(unique(stats::na.omit(mc[[cn]])))
        prev <- input[[paste0("filter_val_", cn)]]
        column(width = 3, selectizeInput(ns(paste0("filter_val_", cn)), label = cn,
                 choices = vals, multiple = TRUE,
                 selected = if (is.null(prev)) character(0) else intersect(prev, vals),
                 options = list(placeholder = "All")))
      }))
    })

    # Group keys are annotation-only, so they are computed once per WELL
    # (hundreds of rows) instead of once per reading row (hundreds of
    # thousands), then joined onto the curve tables.
    well_keys <- reactive({
      mf <- merged_tm() %>% inner_join(filtered_wells(), by = c("plate", "well"))
      mf$.group_key <- paste_key(mf, group_cols_active())
      dplyr::distinct(mf, plate, well, .group_key)
    })

    # Filtered, display-thinned curve tables carrying the group key. The
    # optional loess smoothing here is DISPLAY-ONLY: deri_plot() and the
    # App dF/dT smooth keep working from the raw curves.
    mc_plot <- reactive({
      tbl <- mc_tidy_tm() %>%
        inner_join(filtered_wells(), by = c("plate", "well")) %>%
        inner_join(well_keys(), by = c("plate", "well"))
      if (isTRUE(input$smooth_curves)) tbl <- smooth_display(tbl)
      thin_display(tbl, "fluorescence", n_bins = 120)
    })

    deri_plot <- reactive({
      tbl <- dfdt_source() %>%
        inner_join(filtered_wells(), by = c("plate", "well")) %>%
        inner_join(well_keys(), by = c("plate", "well"))
      thin_display(tbl, "derivative", n_bins = 120)
    })

    # Filtered per-well results (well counts / Tm stats in the row headers)
    merged_filt <- reactive({
      merged_tm() %>% inner_join(filtered_wells(), by = c("plate", "well"))
    })

    groups <- reactive({
      req(data())
      sort(unique(stats::na.omit(mc_plot()$.group_key)))
    })

    # =========================================================================
    # Per-group reference rows. With "Match reference by target" on, only the
    # group's own target's reference wells are used: each target has its own
    # baseline, so comparing every target against one fixed target's reference
    # would be misleading. NULL when the group is the reference itself.
    # =========================================================================
    group_ref_rows <- function(rd_full, tbl, cols, g, key) {
      if (is.null(rd_full) || !length(cols)) return(NULL)
      if (ref_is_self(ref_def(), cols, g)) return(NULL)
      rows <- rd_full
      if (isTRUE(input$overlay_by_target) && "target" %in% colnames(rows)) {
        idx <- !is.na(key) & key %in% g
        tg <- unique(tbl$target[idx]); tg <- tg[!is.na(tg)]
        if (length(tg) == 1) {
          rows <- rows[!is.na(rows$target) & rows$target == tg, , drop = FALSE]
        }
      }
      if (nrow(rows) == 0) NULL else rows
    }

    # Median-Tm marker lines for one group's dF/dT plot: one entry per colour
    # subgroup (e.g. each ligand) plus one for the row's reference group,
    # each the median of the method-selected Tm over that subgroup's wells.
    # `df` / `ref_df` are the colour-prepared curve rows of the group;
    # `ref_label` matches the grey reference's legend entry.
    tm_lines_for <- function(df, ref_df, ref_label = NULL) {
      wells <- dplyr::distinct(df, plate, well)
      mf <- merged_filt() %>% inner_join(wells, by = c("plate", "well"))
      sub <- if (".color_key" %in% colnames(df)) {
        wk <- dplyr::distinct(df, plate, well, .color_key)
        mf %>% inner_join(wk, by = c("plate", "well")) %>%
          filter(!is.na(.color_key)) %>%
          group_by(.color_key) %>%
          summarise(tm = suppressWarnings(median(tm, na.rm = TRUE)), .groups = "drop") %>%
          filter(!is.na(tm)) %>%
          rename(value = .color_key)
      } else {
        mf %>%
          summarise(tm = suppressWarnings(median(tm, na.rm = TRUE))) %>%
          filter(!is.na(tm)) %>%
          mutate(value = "median")
      }
      if (!is.null(ref_df) && nrow(ref_df) > 0) {
        rw <- dplyr::distinct(ref_df, plate, well)
        ref_tm <- suppressWarnings(median(
          (merged_filt() %>% inner_join(rw, by = c("plate", "well")))$tm, na.rm = TRUE))
        if (is.finite(ref_tm)) {
          rv <- if (!is.null(ref_label)) ref_label else "ref"
          # when the reference is itself one of the row's subgroups (e.g. the
          # DMSO wells of this target), its own line already marks the spot
          if (!(sub(" \\(ref\\)$", "", rv) %in% sub$value)) {
            sub <- bind_rows(sub, tibble(tm = ref_tm, value = rv))
          }
        }
      }
      sub
    }

    # =========================================================================
    # Per-well results / Tm summary
    # =========================================================================
    tm_summary <- reactive({
      req(data())
      df <- merged_tm()
      grp <- intersect(c("plate", "target", "ligand"), colnames(df))
      df %>%
        group_by(across(all_of(grp))) %>%
        summarise(
          n_wells = dplyr::n(),
          tm_median = round(suppressWarnings(median(tm, na.rm = TRUE)), 2),
          .groups = "drop"
        ) %>%
        mutate(tm_median = dplyr::if_else(is.nan(tm_median), NA_real_, tm_median)) %>%
        arrange(across(all_of(grp)))
    })

    # =========================================================================
    # Update controls after data is loaded
    # =========================================================================
    observeEvent(data(), {
      d <- data()

      tr <- suppressWarnings(range(d$mc_tidy$temperature, na.rm = TRUE))
      if (all(is.finite(tr))) {
        updateSliderInput(session, "temp_range",
                          min = floor(tr[1]), max = ceiling(tr[2]),
                          value = c(floor(tr[1]), ceiling(tr[2])))
      }

      av <- d$annot_vars
      default_group <- intersect(c("target", "ligand", "well"), av)
      updateSelectInput(session, "row_group",
                        choices = av,
                        selected = if (length(default_group)) default_group[1] else av[1])
      default_color <- intersect(c("ligand", "target"), av)
      updateSelectInput(session, "color_var", choices = av,
                        selected = if (length(default_color)) default_color[1] else character(0))
      updateSelectInput(session, "ref_cols",
                        choices = av,
                        selected = intersect(c("ligand", "sample_name"), av)[1])
    })

    # =========================================================================
    # Right side: Data tab + Plots tab
    # Re-renders only when new data arrives, so changing plot settings keeps
    # the current tab (and scroll position) instead of jumping back to Data.
    # =========================================================================
    output$results_ui <- renderUI({
      req(data())

      tabsetPanel(type = "tabs",
        tabPanel("Data", br(),
          tabsetPanel(type = "tabs",
            tabPanel("Per-well Results", br(), DT::dataTableOutput(ns("data_preview"))),
            tabPanel("Tm Summary", br(), DT::dataTableOutput(ns("tm_summary_tbl")))
          )
        ),
        tabPanel("Plots", br(),
          uiOutput(ns("filter_ui")),
          uiOutput(ns("plots_ui"))
        )
      )
    })

    # Row layout is paginated (5 group rows per page) so switching to the
    # Plots tab only waits for ~10 plots instead of every group at once.
    # Every page still obeys the current well filter.
    page_size <- 5L
    rv <- reactiveValues(page = 1L)
    n_pages <- reactive({
      req(data())
      max(1L, ceiling(length(groups()) / page_size))
    })
    visible_groups <- reactive({
      gs <- groups()
      req(length(gs) > 0)
      pg <- min(max(as.integer(rv$page), 1L), n_pages())
      gs[((pg - 1L) * page_size + 1L):pmin(pg * page_size, length(gs))]
    })
    observeEvent(data(), { rv$page <- 1L })
    observeEvent(input$row_group, { rv$page <- 1L })
    observe({
      req(data())
      pg <- min(max(as.integer(rv$page), 1L), n_pages())
      if (!identical(as.integer(rv$page), pg)) rv$page <- pg
    })
    observeEvent(input$prev_page, {
      pg <- max(1L, as.integer(rv$page) - 1L)
      if (!identical(as.integer(rv$page), pg)) rv$page <- pg
    })
    observeEvent(input$next_page, {
      pg <- min(n_pages(), as.integer(rv$page) + 1L)
      if (!identical(as.integer(rv$page), pg)) rv$page <- pg
    })

    output$plots_ui <- renderUI({
      req(data())
      n <- length(groups())
      if (n == 0) {
        return(helpText("No curves match the current filter - relax the filter values or clear some filter columns."))
      }
      if (identical(input$facet_mode, "facet")) {
        pg <- min(max(as.integer(rv$page), 1L), n_pages())
        rows <- lapply(seq_along(visible_groups()), function(i) {
          div(style = "border-top: 2px solid #eee; margin-top: 14px; padding-top: 6px;",
              h5(textOutput(ns(paste0("row_header_", i))), style = "margin: 0 0 2px 0;"),
              fluidRow(
                column(width = 6, plotOutput(ns(paste0("raw_", i)), height = "330px")),
                column(width = 6, plotOutput(ns(paste0("deri_", i)), height = "330px"))
              ))
        })
        tagList(
          fluidRow(
            column(width = 1, actionButton(ns("prev_page"), "", icon = icon("chevron-left"))),
            column(width = 3, div(style = "text-align: center; padding-top: 8px;",
                                  textOutput(ns("page_indicator")))),
            column(width = 1, actionButton(ns("next_page"), "", icon = icon("chevron-right")))
          ),
          div(style = "max-height: 900px; overflow-y: auto;", rows)
        )
      } else {
        fluidRow(
          column(width = 6, plotOutput(ns("plot_merged_raw"), height = "430px")),
          column(width = 6, plotOutput(ns("plot_merged_deri"), height = "430px"))
        )
      }
    })

    output$page_indicator <- renderText({
      req(data())
      sprintf("Page %d / %d  (%d groups)", rv$page, n_pages(), length(groups()))
    })

    # ---- merged layout: one plot with every filtered curve ------------------
    merged_color_cols <- reactive({
      if (length(input$color_var)) input$color_var
      else if (length(group_cols_active())) ".group_key"
      else character(0)
    })

    output$plot_merged_raw <- renderPlot({
      df <- mc_plot()
      req(nrow(df) > 0)
      cp <- color_prep(df, merged_color_cols())
      if (identical(merged_color_cols(), ".group_key")) cp$name <- "Group"
      ref <- mc_ref_rows_display()
      ref_cc <- if (length(merged_color_cols()) == 1 && !is.null(ref) &&
                    merged_color_cols() %in% colnames(ref) &&
                    !is.numeric(df[[merged_color_cols()]])) merged_color_cols() else NULL
      build_curve_ggplot(cp$df, "fluorescence", cp$col, cp$name,
                         input$temp_range, ref,
                         ref_label = ref_legend_label(ref, ref_cc))
    }, res = 96)

    output$plot_merged_deri <- renderPlot({
      df <- deri_plot()
      req(nrow(df) > 0)
      cp <- color_prep(df, merged_color_cols())
      if (identical(merged_color_cols(), ".group_key")) cp$name <- "Group"
      ref <- deri_ref_rows()
      ref_cc <- if (length(merged_color_cols()) == 1 && !is.null(ref) &&
                    merged_color_cols() %in% colnames(ref) &&
                    !is.numeric(df[[merged_color_cols()]])) merged_color_cols() else NULL
      build_curve_ggplot(cp$df, "derivative", cp$col, cp$name,
                         input$temp_range, ref,
                         ref_label = ref_legend_label(ref, ref_cc))
    }, res = 96)

    # ---- one row per group: fluorescence left, dF/dT right ------------------
    # Only the current page's slots (max 5) are registered; the render
    # expressions read visible_groups() reactively, so a page flip swaps the
    # contents in place and every page still obeys the well filter.

    # Row header: group key + well count + median Tm, plus the reference Tm
    # and the group's delta Tm vs it. Shared by the on-screen row headers and
    # the PDF page titles so both views match.
    row_header_text <- function(g) {
      base <- group_header_label(merged_filt(), group_cols_active(), g)
      tbl <- mc_plot()
      ref <- group_ref_rows(mc_ref_rows(), tbl, group_cols_active(), g, tbl$.group_key)
      if (is.null(ref)) return(base)
      rw <- dplyr::distinct(ref, plate, well)
      ref_tm <- suppressWarnings(median(
        (merged_filt() %>% inner_join(rw, by = c("plate", "well")))$tm, na.rm = TRUE))
      mf <- merged_filt()
      k <- paste_key(mf, group_cols_active())
      grp_tm <- suppressWarnings(median(mf$tm[!is.na(k) & k %in% g], na.rm = TRUE))
      if (!is.finite(ref_tm) || !is.finite(grp_tm)) return(base)
      paste0(base, sprintf("  |  ref Tm = %.2f \u00b0C  |  delta Tm = %+.2f \u00b0C",
                           ref_tm, grp_tm - ref_tm))
    }

    observeEvent(list(input$facet_mode, length(groups())), {
      n <- length(groups())
      if (!n || !identical(input$facet_mode, "facet")) return()
      for (ii in seq_len(min(page_size, n))) {
        local({
          i <- ii
          output[[paste0("row_header_", i)]] <- renderText({
            vg <- visible_groups()
            req(length(vg) >= i)
            row_header_text(vg[[i]])
          })
          output[[paste0("raw_", i)]] <- renderPlot({
            vg <- visible_groups()
            req(length(vg) >= i)
            g <- vg[[i]]
            tbl <- mc_plot()
            cols <- group_cols_active()
            df <- group_curve_df(tbl, cols, g, key = tbl$.group_key)
            ref <- group_ref_rows(mc_ref_rows_display(), tbl, cols, g, tbl$.group_key)
            cp <- color_prep(df, input$color_var)
            ref_cc <- if (length(input$color_var) == 1 && !is.null(ref) &&
                          input$color_var %in% colnames(ref) &&
                          !is.numeric(df[[input$color_var]])) input$color_var else NULL
            build_curve_ggplot(cp$df, "fluorescence", cp$col, cp$name,
                               input$temp_range, ref,
                               ref_label = ref_legend_label(ref, ref_cc))
          }, res = 96)
          output[[paste0("deri_", i)]] <- renderPlot({
            vg <- visible_groups()
            req(length(vg) >= i)
            g <- vg[[i]]
            tbl <- deri_plot()
            cols <- group_cols_active()
            df <- group_curve_df(tbl, cols, g, key = tbl$.group_key)
            ref <- group_ref_rows(deri_ref_rows(), tbl, cols, g, tbl$.group_key)
            cp <- color_prep(df, input$color_var)
            ref_cc <- if (length(input$color_var) == 1 && !is.null(ref) &&
                          input$color_var %in% colnames(ref) &&
                          !is.numeric(df[[input$color_var]])) input$color_var else NULL
            ref_lab <- ref_legend_label(ref, ref_cc)
            lines <- tm_lines_for(cp$df, ref, ref_lab)
            build_curve_ggplot(cp$df, "derivative", cp$col, cp$name,
                               input$temp_range, ref, tm_lines = lines,
                               ref_label = ref_lab)
          }, res = 96)
        })
      }
    })

    # =========================================================================
    # Data tab tables
    # =========================================================================
    output$data_preview <- DT::renderDataTable({
      req(data())
      df <- merged_tm() %>%                   # tm columns are already rounded to 2 dp
        mutate(across(any_of("conc"), ~ round(.x, 3)))
      tbl <- DT::datatable(df, rownames = FALSE,
                           options = list(pageLength = 25, scrollX = TRUE))
      tm_cols <- grep("^tm", names(df), value = TRUE)
      if (length(tm_cols)) DT::formatRound(tbl, columns = tm_cols, digits = 2) else tbl
    })

    output$tm_summary_tbl <- DT::renderDataTable({
      req(data())
      DT::datatable(tm_summary(), rownames = FALSE,
                    options = list(pageLength = 25, scrollX = TRUE)) %>%
        DT::formatRound(columns = "tm_median", digits = 2)
    })

    # =========================================================================
    # Download Data (raw + analysis + plot data)
    # =========================================================================
    output$download_data <- downloadHandler(
      filename = function() {
        paste0('DSF_', Sys.Date(), '_data.xlsx')
      },
      content = function(file) {
        req(data())
        per_well <- merged_tm()
        tm_sum <- tm_summary()

        # With a reference configured (overlay on + columns/values chosen),
        # add each row's reference Tm and its delta Tm. Without a reference,
        # no ref information is emitted anywhere.
        ref_rows_tbl <- mc_ref_rows()
        ref_set <- !is.null(ref_rows_tbl) && nrow(ref_rows_tbl) > 0 &&
          "target" %in% colnames(merged_tm())
        if (ref_set) {
          rw <- dplyr::distinct(ref_rows_tbl, plate, well)
          ref_tm_tbl <- merged_tm() %>% inner_join(rw, by = c("plate", "well"))
          if (isTRUE(input$overlay_by_target)) {
            # each target is compared against its own reference wells
            ref_tm_tbl <- ref_tm_tbl %>%
              group_by(target) %>%
              summarise(`ref Tm` = round(suppressWarnings(median(tm, na.rm = TRUE)), 2),
                        .groups = "drop")
            per_well <- per_well %>% left_join(ref_tm_tbl, by = "target")
            tm_sum <- tm_sum %>% left_join(ref_tm_tbl, by = "target")
          } else {
            # a single fixed reference for every well
            one <- tibble(`ref Tm` = round(
              suppressWarnings(median(ref_tm_tbl$tm, na.rm = TRUE)), 2))
            per_well <- per_well %>% crossing(one)
            tm_sum <- tm_sum %>% crossing(one)
          }
          per_well <- per_well %>%
            mutate(`delta Tm` = round(tm - `ref Tm`, 2)) %>%
            relocate(`ref Tm`, `delta Tm`, .after = last_col())
          tm_sum <- tm_sum %>%
            mutate(`delta Tm` = round(tm_median - `ref Tm`, 2)) %>%
            relocate(`ref Tm`, `delta Tm`, .after = last_col())
        }

        write_xlsx_sheets(list(
          "Per-well Results" = per_well,
          "Tm Summary" = tm_sum,
          "Metadata" = data()$meta,
          "Raw Fluorescence" = mc_tidy_tm(),
          "1st Derivative" = data()$derivatives
        ), file)
      }
    )

    # =========================================================================
    # Download Plots. Zip with two output families (when the data carries
    # target/ligand annotations):
    #   DeltaTm_overview.pdf  one page per target: each ligand's delta Tm vs
    #                         the reference (wells in blue, median in black)
    #   Plots_<target>.pdf    one fixed-size page per ligand: RFU left, dF/dT
    #                         right, all wells of the ligand in dark blue, the
    #                         reference in grey, dashed median-Tm marker, the
    #                         combined Tm in the page title
    # Mirrors the on-screen view: current well filter, smoothing, reference
    # overlay (with match-by-target) and temperature range all apply. Without
    # target/ligand annotations the previous per-group export is used instead.
    # =========================================================================
    build_pdf_page <- function(g, cols, key_raw, key_deri, ref_raw, ref_deri,
                               mc_df, deri_df, color_var, temp_range, title) {
      raw_df <- group_curve_df(mc_df, cols, g, key = key_raw)
      deri_df2 <- group_curve_df(deri_df, cols, g, key = key_deri)
      ref_r <- group_ref_rows(ref_raw, mc_df, cols, g, key_raw)
      ref_d <- group_ref_rows(ref_deri, deri_df, cols, g, key_deri)
      ref_cc <- if (length(color_var) == 1 && !is.null(ref_r) &&
                    color_var %in% colnames(ref_r) &&
                    !is.numeric(raw_df[[color_var]])) color_var else NULL
      ref_lab_d <- ref_legend_label(ref_d, ref_cc)
      cp_raw <- color_prep(raw_df, color_var)
      cp_deri <- color_prep(deri_df2, color_var)
      p1 <- build_curve_ggplot(cp_raw$df, "fluorescence", cp_raw$col, cp_raw$name,
                               temp_range, ref_r,
                               ref_label = ref_legend_label(ref_r, ref_cc))
      p2 <- build_curve_ggplot(cp_deri$df, "derivative", cp_deri$col, cp_deri$name,
                               temp_range, ref_d,
                               tm_lines = tm_lines_for(cp_deri$df, ref_d, ref_lab_d),
                               ref_label = ref_lab_d)
      p1 + p2 + patchwork::plot_annotation(
        title = title,
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

        mc_df <- mc_plot()
        deri_df <- deri_plot()
        if (!nrow(mc_df) || !nrow(deri_df)) {
          showNotification("No curves available to plot.", type = "error", duration = 8)
          return(NULL)
        }
        mt <- merged_filt()                      # per-well results, current filter
        has_tl <- all(c("target", "ligand") %in% colnames(mt)) &&
          all(c("target", "ligand") %in% colnames(mc_df))

        # ---- helpers for the per-target / delta-Tm export -----------------
        # Reference curves for one target: match-by-target when enabled,
        # otherwise the full reference set; NULL when this page IS the
        # reference ligand (its own grey overlay would double-draw it).
        ref_rows_for_target <- function(rd_full, tgt, lg) {
          if (is.null(rd_full) || !nrow(rd_full)) return(NULL)
          rd <- ref_def()
          if (!is.null(rd)) {
            i <- which(rd$cols == "ligand")
            if (length(i) == 1 && identical(rd$vals[[i]], lg)) return(NULL)
          }
          if (isTRUE(input$overlay_by_target) && "target" %in% colnames(rd_full)) {
            out <- rd_full[!is.na(rd_full$target) & rd_full$target == tgt, , drop = FALSE]
            if (!nrow(out)) return(NULL)
            out
          } else rd_full
        }
        nat_lig_order <- function(lgs) {
          ls <- unique(lgs)
          ord <- stringr::str_sort(ls, numeric = TRUE, locale = "en")
          if ("DMSO" %in% ls) ord <- c("DMSO", setdiff(ord, "DMSO"))
          ord
        }
        # one page per ligand: RFU left, dF/dT right; every replicate of the
        # ligand group in one dark blue, the reference in grey. The dF/dT
        # panel marks the ligand's median Tm (dashed, labelled) and the
        # reference Tm (dotdash); the title reports Tm and delta Tm. Without
        # a reference for the target, every reference element is skipped.
        plot_ligand_page <- function(tgt, lg, tm_med, delta, ref_tm, tr) {
          rows <- mc_df %>% filter(.data$target == !!tgt, .data$ligand == !!lg)
          dvs  <- deri_df %>% filter(.data$target == !!tgt, .data$ligand == !!lg)
          if (!nrow(rows) || !nrow(dvs)) return(NULL)
          ref_r <- ref_rows_for_target(mc_ref_rows_display(), tgt, lg)
          ref_d <- ref_rows_for_target(deri_ref_rows(), tgt, lg)
          has_ref <- is.finite(ref_tm)
          ttl <- paste0(tgt, " \u2014 ", lg)
          if (is.finite(tm_med)) {
            ttl <- paste0(ttl, "   (Tm ", sprintf("%.1f", tm_med), " \u00b0C")
            if (has_ref)
              ttl <- paste0(ttl, ", \u0394Tm ", sprintf("%+.1f", tm_med - ref_tm), " \u00b0C")
            ttl <- paste0(ttl, ")")
          }
          pth <- function(df, y, col, wd) {
            geom_path(data = df,
                      mapping = aes(x = .data$temperature, y = .data[[y]],
                                    group = interaction(plate, well)),
                      colour = col, linewidth = wd, na.rm = TRUE)
          }

          p_fl <- ggplot() +
            { if (!is.null(ref_r)) pth(ref_r, "fluorescence", "grey60", 0.35) } +
            pth(rows, "fluorescence", "#08519C", 0.45) +
            labs(title = ttl, x = NULL, y = "RFU") +
            { if (length(tr) == 2 && !anyNA(tr)) coord_cartesian(xlim = tr) } +
            theme_bw(base_size = 10) +
            theme(plot.title = element_text(face = "bold", size = 13),
                  axis.text.x = element_blank(), axis.ticks.x = element_blank())
          p_dv <- ggplot() +
            { if (!is.null(ref_d)) pth(ref_d, "derivative", "grey60", 0.35) } +
            pth(dvs, "derivative", "#08519C", 0.45) +
            { if (has_ref) geom_vline(xintercept = ref_tm, linetype = "longdash",
                                      colour = "grey45", linewidth = 0.4) } +
            { if (is.finite(tm_med)) list(
                geom_vline(xintercept = tm_med, linetype = "dashed",
                           colour = "#08519C", linewidth = 0.4),
                geom_text(data = data.frame(x = tm_med),
                          aes(x = x, y = Inf, label = sprintf("%.1f", tm_med)),
                          inherit.aes = FALSE, vjust = 1.4, hjust = -0.08,
                          size = 3.4, colour = "#08519C")) } +
            labs(x = "Temperature (\u00b0C)", y = "dF/dT",
                 caption = paste0("dark blue = ligand (all wells overlaid); ",
                                  "grey = reference; dashed = ligand Tm (labelled)",
                                  if (has_ref) " and reference Tm (dotdash)")) +
            { if (length(tr) == 2 && !anyNA(tr)) coord_cartesian(xlim = tr) } +
            theme_bw(base_size = 10) +
            theme(plot.caption = element_text(colour = "grey40", size = 7))
          p_fl + p_dv + patchwork::plot_layout(widths = c(1, 1))
        }
        # one page per target: delta Tm per well/ligand vs the reference
        delta_overview_page <- function(tgt, dtab, ref_lg) {
          dd <- dtab %>% filter(.data$target == !!tgt, !is.na(ligand),
                                !is.na(.data[["delta Tm"]]))
          if (!is.na(ref_lg)) dd <- filter(dd, ligand != ref_lg)
          ls <- nat_lig_order(unique(dd$ligand))
          if (!length(ls)) return(NULL)
          dd <- dd %>% mutate(ligand = factor(ligand, levels = rev(ls)))
          dm <- dd %>% group_by(ligand) %>%
            summarise(d = suppressWarnings(median(.data[["delta Tm"]], na.rm = TRUE)),
                      .groups = "drop")
          ggplot(dd, aes(.data[["delta Tm"]], ligand)) +
            geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
            geom_point(colour = "#08519C", size = 1.8, alpha = 0.7, na.rm = TRUE) +
            geom_point(data = dm, aes(x = d), colour = "black", size = 2.4,
                       na.rm = TRUE) +
            labs(title = paste0(tgt, " \u2014 \u0394Tm vs reference (\u00b0C; ",
                                "black = ligand median, blue = wells)"),
                 x = "\u0394Tm (\u00b0C)", y = NULL) +
            theme_bw(base_size = 10)
        }

        # all groups matching the current filter, one page per group
        gs <- groups()
        if (length(gs) == 0) {
          showNotification("No groups available to plot.", type = "error", duration = 8)
          return(NULL)
        }

        cols <- group_cols_active()
        key_raw <- mc_df$.group_key
        key_deri <- deri_df$.group_key
        ref_raw <- mc_ref_rows_display()
        ref_deri <- deri_ref_rows()
        color_var <- input$color_var
        tr <- input$temp_range

        plot_dir <- tempfile(pattern = "dsf_plots_")
        dir.create(plot_dir, recursive = TRUE)

        # progress feedback: download handlers run outside the normal render
        # loop, so drive a Progress panel bound to the session explicitly and
        # keep a sticky notification in case the panel cannot be displayed
        msg_id <- paste0("pdf_gen_", session$ns(""))
        showNotification("Generating PDF export - this can take a while...",
                         type = "message", duration = NULL, id = msg_id)
        progress <- shiny::Progress$new(session = session, min = 0, max = 1)
        progress$set(message = "Generating vector PDFs", value = 0)

        tryCatch({
          dev_fn <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

          if (has_tl) {
            # ---- reference Tm (shared by the ligand pages and the overview) --
            ref_rows_tbl <- mc_ref_rows()
            ref_set <- !is.null(ref_rows_tbl) && nrow(ref_rows_tbl) > 0
            ref_lg <- NA_character_
            ref_tm_tbl <- NULL
            if (ref_set) {
              rd <- ref_def()
              i <- which(rd$cols == "ligand")
              if (length(i) == 1) ref_lg <- rd$vals[[i]]
              rw <- dplyr::distinct(ref_rows_tbl, plate, well)
              ref_tbl <- mt %>% inner_join(rw, by = c("plate", "well"))
              if (isTRUE(input$overlay_by_target)) {
                ref_tm_tbl <- ref_tbl %>% group_by(target) %>%
                  summarise(`ref Tm` = round(suppressWarnings(median(tm, na.rm = TRUE)), 2),
                            .groups = "drop")
              } else {
                ref_tm_tbl <- tibble(target = NA_character_,
                                     `ref Tm` = round(suppressWarnings(median(ref_tbl$tm, na.rm = TRUE)), 2))
              }
            }
            ref_tm_for <- function(tgt) {
              if (!ref_set || is.null(ref_tm_tbl)) return(NA_real_)
              hit <- ref_tm_tbl$`ref Tm`[is.na(ref_tm_tbl$target) |
                                           ref_tm_tbl$target == tgt]
              if (length(hit)) hit[1] else NA_real_
            }

            # ---- per-target files: one page per ligand --------------------
            targets <- sort(unique(mc_df$target[!is.na(mc_df$target)]))
            for (tgt in targets) {
              ligs <- nat_lig_order(unique(mc_df$ligand[!is.na(mc_df$ligand) &
                                                          mc_df$target == tgt]))
              ref_tm <- ref_tm_for(tgt)
              dev_fn(file.path(plot_dir, paste0("Plots_", sanitize_filename(tgt), ".pdf")),
                     width = 10, height = 4)
              tryCatch(
                for (lg in ligs) {
                  tm_med <- suppressWarnings(median(
                    mt$tm[!is.na(mt$target) & mt$target == tgt &
                            !is.na(mt$ligand) & mt$ligand == lg], na.rm = TRUE))
                  page_ref_tm <- if (identical(lg, ref_lg)) NA_real_ else ref_tm
                  delta <- if (is.finite(tm_med) && is.finite(page_ref_tm))
                    round(tm_med - page_ref_tm, 2) else NA_real_
                  p <- plot_ligand_page(tgt, lg, tm_med, delta, page_ref_tm, tr)
                  if (!is.null(p)) print(p)
                },
                error = function(e) showNotification(
                  paste0("Some pages failed for ", tgt, ": ", e$message),
                  type = "warning", duration = 8),
                finally = grDevices::dev.off())
              done <- 0.9 * which(targets == tgt) / length(targets)
              progress$set(value = done, detail = tgt)
            }

            # ---- delta Tm overview (needs a configured reference) ---------
            if (ref_set) {
              dtab <- if (isTRUE(input$overlay_by_target)) {
                mt %>% left_join(ref_tm_tbl %>% filter(!is.na(target)), by = "target") %>%
                  mutate(`delta Tm` = round(tm - `ref Tm`, 2))
              } else {
                mt %>% mutate(`delta Tm` = round(tm - ref_tm_tbl$`ref Tm`[1], 2))
              }
              max_ligs <- max(vapply(targets, function(tg) {
                d <- dtab %>% filter(.data$target == tg, !is.na(ligand))
                if (!is.na(ref_lg)) d <- filter(d, ligand != ref_lg)
                nrow(distinct(d, ligand))
              }, numeric(1)), 1)
              progress$set(value = 0.95, detail = "Delta Tm overview...")
              dev_fn(file.path(plot_dir, "DeltaTm_overview.pdf"),
                     width = 8, height = max(6, 0.28 * max_ligs + 1.5))
              tryCatch(
                for (tgt in targets) {
                  p <- delta_overview_page(tgt, dtab, ref_lg)
                  if (!is.null(p)) print(p)
                },
                error = function(e) showNotification(
                  paste("Delta Tm overview failed:", e$message),
                  type = "warning", duration = 8),
                finally = grDevices::dev.off())
            }
            progress$set(value = 0.97, detail = "Writing PDFs...")
            pdf_count <- length(list.files(plot_dir, pattern = "\\.pdf$"))
          } else {
            # ---- no target/ligand annotations: previous per-group export ---
            pages <- vector("list", length(gs))
            for (i in seq_along(gs)) {
              pages[[i]] <- build_pdf_page(gs[[i]], cols, key_raw, key_deri,
                                           ref_raw, ref_deri, mc_df, deri_df,
                                           color_var, tr, title = row_header_text(gs[[i]]))
              progress$set(value = 0.9 * i / length(gs),
                           detail = sprintf("%s (%d/%d)", gs[[i]], i, length(gs)))
            }
            progress$set(value = 0.95, detail = "Writing PDFs...")

            combined <- file.path(plot_dir, "000_all_groups.pdf")
            dev_fn(combined, width = 14, height = 6)
            tryCatch(
              for (p in pages) print(p),
              error = function(e) showNotification(paste("Some PDF pages failed:", e$message),
                                                   type = "warning", duration = 8),
              finally = grDevices::dev.off()
            )
            for (i in seq_along(gs)) {
              ggsave(file.path(plot_dir, sprintf("%03d_%s.pdf", i, sanitize_filename(gs[[i]]))),
                     plot = pages[[i]], width = 14, height = 6, device = dev_fn)
            }
            pdf_count <- length(list.files(plot_dir, pattern = "\\.pdf$"))
          }

          pdf_files <- list.files(plot_dir, pattern = "\\.pdf$")
          zip::zip(zipfile = file, files = pdf_files, root = plot_dir)
          progress$close()
          removeNotification(msg_id)
          showNotification(sprintf("PDF export ready: %d files.", pdf_count),
                           type = "message", duration = 6)
        }, error = function(e) {
          progress$close()
          removeNotification(msg_id)
          showNotification(paste("Error creating zip:", e$message),
                           type = "error", duration = 10)
        })
      },
      contentType = "application/zip"
    )
  })
}

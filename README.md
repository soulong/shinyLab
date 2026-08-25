# shinyBioTools

A standalone R Shiny application for biological data analysis. It bundles seven analysis tools in a single `shinydashboard` interface: qPCR analysis, shRNA/sgRNA primer design, score normalization, drug synergy, thermal shift (DSF), and an end-to-end RNA-seq pipeline.

## Features

- **Real-time PCR**: Analyze qPCR data with the delta-delta-CT method, with faceted / combined bar plots, melt-curve visualization, and Excel export
- **Easy shRNA**: Convert gene IDs (SYMBOL/Ensembl/Entrez), query the splashRNA database, and design overlapping-PCR primers for miR-30 based vectors
- **Easy sgRNA**: Design LentiCRISPRv2 cloning primers with customizable prefix/suffix sequences (forward and reverse)
- **Score Norm**: Rank-based inverse normal transformation (van der Waerden INT) with summary statistics, side-by-side histograms, and CSV export
- **SynergyFinder**: Drug synergy scoring with ZIP, HSA, Bliss, and Loewe models, including dose-response, 2D synergy, surface, and barometer plots
- **DSF Analysis**: Differential scanning fluorimetry — Tm (derivative or Boltzmann) and delta-Tm scatter plots, raw/derivative curves, reference-ligand handling, and optional chemical structures from SMILES
- **RNA-seq**: End-to-end DESeq2 pipeline — data upload with editable metadata, QC (PCA / sample correlation / variable-gene heatmap), Wald comparisons, volcano plots, and ORA + GSEA enrichment against MSigDB

## Quick Start

Run the app from the project root (the app sources modules via relative paths):

```r
# In RStudio: open app.R and click "Run App"
source("app.R")
```

Or from terminal:

```bash
Rscript app.R
```

The app listens on `http://127.0.0.1:5005` (configured in `app.R`); a browser window opens automatically when launched from RStudio. Open that address in your browser if it does not.

## Requirements

- R >= 4.5.0
- R packages (all directly used by the code):

**Base / UI**
- shiny, shinydashboard, tidyverse (includes ggplot2, dplyr, tidyr, readr, stringr, purrr, tibble, forcats), magrittr, DT, shinyjs, plotly

**Real-time PCR**
- readxl, writexl, janitor

**shRNA / sgRNA**
- Biostrings, httr, rvest

**Score Norm**
- (base R only)

**SynergyFinder**
- synergyfinder, patchwork, writexl, zip

**DSF Analysis**
- zip, janitor, mgcv, gratia, patchwork, DT, writexl
- Optional: rcdk (renders 2D chemical structures from SMILES in the Ligand Details view)

**RNA-seq**
- readxl, writexl, rhandsontable, DESeq2, ashr, corrplot, ggforce, ggrepel, clusterProfiler, msigdbr, enrichplot, future, furrr, matrixStats

**Gene ID conversion**
- AnnotationDbi, org.Hs.eg.db (human), org.Mm.eg.db (mouse)

Implicit dependencies (installed automatically): `scales`, `rlang`, `gtable` (via ggplot2), `xml2` (via rvest), `parallel` (base).

Install packages (`BiocManager::install` handles both CRAN and Bioconductor packages):

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "shiny", "shinydashboard", "tidyverse", "magrittr", "patchwork",
  "writexl", "readxl", "janitor", "httr", "rvest",
  "zip", "mgcv", "gratia", "DT", "shinyjs", "plotly",
  "corrplot", "ggforce", "ggrepel", "future", "furrr", "matrixStats",
  "Biostrings", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db",
  "DESeq2", "ashr", "clusterProfiler", "msigdbr", "enrichplot",
  "synergyfinder"
))
```

Note: `rhandsontable` is archived on CRAN; install it from the CRAN archive if needed:

```r
install.packages("https://cran.r-project.org/src/contrib/Archive/rhandsontable/rhandsontable_0.3.8.tar.gz", repos = NULL, type = "source")
```

## Project Layout

| Sidebar item | Tab `tabName` | Module file |
|---|---|---|
| Real-time PCR | `rtpcr` | `R/mod_rtpcr.R` |
| Easy shRNA | `shrna` | `R/mod_shrna.R` |
| Easy sgRNA | `sgrna` | `R/mod_sgrna.R` |
| Score Norm | `scorenorm` | `R/mod_scorenorm.R` |
| SynergyFinder | `synergy` | `R/mod_synergy.R` |
| DSF Analysis | `dsf` | `R/mod_dsf.R` |
| RNA-seq | `rnaseq` | `R/mod_rnaseq.R` |

Shared helpers live in `R/utils.R`; the app entry point is `app.R`. A `manifest.json` is included for deployment to Posit Connect / shinyapps.io via `rsconnect::writeManifest()`.

## Author

Hao He <haohe90@gmail.com>

## License

MIT License
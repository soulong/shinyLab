# shinyBioTools

A standalone R Shiny application for biological data analysis.

## Features

- **Real-time PCR Analysis**: Analyze qPCR data using the delta-delta-CT method with visualization
- **Gene ID Conversion**: Convert between different gene identifiers (Ensembl, Symbol, Entrez)
- **shRNA Primer Design**: Design shRNA primers using splashRNA database queries
- **Drug Synergy Analysis**: Calculate drug synergy scores using ZIP, HSA, Bliss, and Loewe models

## Quick Start

### Option 1: Run Standalone (Recommended)

Simply run the R script directly:

```r
# In RStudio or R console
source("R/app.R")
```

Or from terminal:
```bash
Rscript R/app.R
```

The application will open at `http://localhost:5001`

### Option 2: Package Installation

For package installation, you'll need to rebuild the package structure.

## Requirements

- R >= 4.0.0
- R packages: shiny, shinydashboard, tidyverse, magrittr, patchwork, writexl, readxl, stringr, Biostrings, httr, xml2, limma, tibble, synergyfinder

Install missing packages:
```r
install.packages(c(
  "shiny", "shinydashboard", "tidyverse", "magrittr", "patchwork", "writexl",
  "readxl", "stringr", "Biostrings", "httr", "xml2", "limma", "tibble"
))

# Install synergyfinder from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("synergyfinder")
```

## Usage

### Real-time PCR Analysis

1. Upload a QuantStudio export file (.xls/.xlsx)
2. Select targets and samples to analyze
3. Click "Analyze" to compute delta-delta-CT values
4. View results in Facet Plot, Combined Plot, or Data Summary tabs
5. Check Melting Curve tab for quality control

### shRNA Design

1. Enter gene symbols (one per line)
2. Select species (human/mouse) and input type
3. Click "Convert IDs" to get Entrez IDs
4. Click "Design Primers" to query splashRNA and design primers
5. Download results as Excel file

### Drug Synergy Analysis

1. Upload a formatted CSV file with drug combination data
2. Select effect type (viability or inhibition)
3. Click "Submit" to calculate synergy scores
4. View results in Output and Parameters tabs
5. Download all results as a ZIP file containing PDFs and Excel files

## Version

Current version: 1.5

## Author

Hao He <haohe90@gmail.com>

## License

MIT License

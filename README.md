# shinyBioTools

A standalone R Shiny application for biological data analysis.

## Features

- **Real-time PCR Analysis**: Analyze qPCR data using the delta-delta-CT method with visualization
- **Gene ID Conversion**: Convert between different gene identifiers (Ensembl, Symbol, Entrez)
- **shRNA Primer Design**: Design shRNA primers using splashRNA database queries

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
- R packages: shiny, shinythemes, ggplot2, dplyr, readxl, writexl, stringr, Biostrings, httr, xml2, limma, tibble

Install missing packages:
```r
install.packages(c(
  "shiny", "shinythemes", "ggplot2", "dplyr", "readxl", "writexl",
  "stringr", "Biostrings", "httr", "xml2", "limma", "tibble"
))
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

## File Structure

```
shinyBioTools/
├── R/
│   └── standalone_app.R    # Self-contained Shiny application
├── README.md               # This file
└── data/                   # Data files (for package mode)
```

## Version

Current version: 1.4

## Author

Hao He <haohe90@gmail.com>

## License

MIT License

# shinyBioTools

A standalone R Shiny application for biological data analysis.

## Features

- **Real-time PCR Analysis**: Analyze qPCR data using the delta-delta-CT method with visualization
- **shRNA Primer Design**: Design shRNA primers using splashRNA database queries
- **Drug Synergy Analysis**: Calculate drug synergy scores using ZIP, HSA, Bliss, and Loewe models

## Quick Start

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


## Requirements

- R >= 4.5.0
- R packages: shiny, shinydashboard, tidyverse, magrittr, patchwork, writexl, readxl, stringr, Biostrings, httr, xml2, limma, tibble, synergyfinder

Install packages:
```r
BiocManager::install(c(
  "shiny", "shinydashboard", "tidyverse", "magrittr", "patchwork", "writexl",
  "readxl", "stringr", "Biostrings", "httr", "xml2", "limma", "tibble", "synergyfinder"
))
```

## Author

Hao He <haohe90@gmail.com>

## License

MIT License

# shinyBioTools - Main Application
# Deploy: rsconnect::deployApp()

required_pkgs <- c("shiny", "shinydashboard", "ggplot2", "dplyr", "readxl",
                   "writexl", "Biostrings", "httr", "xml2", "AnnotationDbi",
                   "org.Hs.eg.db", "org.Mm.eg.db")

missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "))
}

for (pkg in required_pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}


options(shiny.maxRequestSize=100*1024^2)


source("R/utils.R")
source("R/mod_rtpcr.R")
source("R/mod_shrna.R")

# UI
ui <- dashboardPage(
  dashboardHeader(title = "shinyBioTools"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Real-time PCR", tabName = "rtpcr", icon = icon("chart-bar")),
      menuItem("Easy shRNA", tabName = "shrna", icon = icon("dna")),
      menuItem("Easy sgRNA", tabName = "sgrna", icon = icon("scissors"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "rtpcr", rtPCRUI("rtPCR")),
      tabItem(tabName = "shrna", shRNAUI("shRNA")),
      tabItem(tabName = "sgrna", h4("Coming soon..."))
    )
  )
)

# Server
server <- function(input, output, session) {
  rtPCRServer("rtPCR")
  shRNAServer("shRNA")
}

shinyApp(ui = ui, server = server)

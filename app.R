
# Deploy: rsconnect::deployApp()

library(shiny)
library(shinydashboard)
library(tidyverse)

source("R/utils.R")
source("R/mod_rtpcr.R")
source("R/mod_shrna.R")
source("R/mod_sgrna.R")
source("R/mod_scorenorm.R")
source("R/mod_synergy.R")
source("R/mod_dsf.R")
source("R/mod_rnaseq.R")

options(shiny.host = "0.0.0.0", 
        shiny.port = 5005, 
        shiny.launch.browser = T, 
        shiny.maxRequestSize=100*1024^2
        )


ui <- dashboardPage(
  dashboardHeader(title = "shinyBioTools"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Real-time PCR", tabName = "rtpcr", icon = icon("chart-bar")),
      menuItem("Easy shRNA", tabName = "shrna", icon = icon("dna")),
      menuItem("Easy sgRNA", tabName = "sgrna", icon = icon("scissors")),
      menuItem("Score Norm", tabName = "scorenorm", icon = icon("sort-numeric-down-alt")),
      menuItem("SynergyFinder", tabName = "synergy", icon = icon("flask")),
      menuItem("DSF Analysis", tabName = "dsf", icon = icon("thermometer-half")),
      menuItem("RNA-seq", tabName = "rnaseq", icon = icon("chart-line"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "rtpcr", rtPCRUI("rtPCR")),
      tabItem(tabName = "shrna", shRNAUI("shRNA")),
      tabItem(tabName = "sgrna", sgRNAUI("sgrna")),
      tabItem(tabName = "scorenorm", scoreNormUI("scorenorm")),
      tabItem(tabName = "synergy", synergyUI("synergy")),
      tabItem(tabName = "dsf", dsfUI("dsf")),
      tabItem(tabName = "rnaseq", rnaseqUI("rnaseq"))
    )
  )
)


server <- function(input, output, session) {
  rtPCRServer("rtPCR")
  shRNAServer("shRNA")
  synergyServer("synergy")
  dsfServer("dsf")
  sgRNAServer("sgrna")
  scoreNormServer("scorenorm")
  rnaseqServer("rnaseq")
}

runApp(shinyApp(ui = ui, server = server))





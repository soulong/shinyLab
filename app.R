
# rstudioapi::getActiveDocumentContext()$path |> dirname() |> setwd()
# rsconnect::writeManifest()
# Deploy: rsconnect::deployApp()

library(shiny)
library(shinydashboard)
library(tidyverse)

source("R/utils.R")
source("R/mod_rtpcr.R")
source("R/mod_shrna.R")
source("R/mod_synergy.R")
source("R/mod_dsf.R")

options(shiny.maxRequestSize=100*1024^2)



ui <- dashboardPage(
  dashboardHeader(title = "shinyBioTools"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Real-time PCR", tabName = "rtpcr", icon = icon("chart-bar")),
      menuItem("Easy shRNA", tabName = "shrna", icon = icon("dna")),
      menuItem("Easy sgRNA", tabName = "sgrna", icon = icon("scissors")),
      menuItem("SynergyFinder", tabName = "synergy", icon = icon("flask")),
      menuItem("DSF Analysis", tabName = "dsf", icon = icon("thermometer-half"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "rtpcr", rtPCRUI("rtPCR")),
      tabItem(tabName = "shrna", shRNAUI("shRNA")),
      tabItem(tabName = "sgrna", h4("Coming soon...")),
      tabItem(tabName = "synergy", synergyUI("synergy")),
      tabItem(tabName = "dsf", dsfUI("dsf"))
    )
  )
)


server <- function(input, output, session) {
  rtPCRServer("rtPCR")
  shRNAServer("shRNA")
  synergyServer("synergy")
  dsfServer("dsf")
}


shinyApp(ui = ui, server = server)

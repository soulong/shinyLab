# Deploy: rsconnect::deployApp()

library(shiny)
library(shinydashboard)
library(tidyverse)

source("R/utils.R")
source("R/mod_rtpcr.R")
source("R/mod_shrna.R")

options(shiny.maxRequestSize=100*1024^2)



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


server <- function(input, output, session) {
  rtPCRServer("rtPCR")
  shRNAServer("shRNA")
}


shinyApp(ui = ui, server = server)

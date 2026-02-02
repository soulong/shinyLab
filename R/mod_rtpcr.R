# shinyBioTools - rtPCR Module
library(readxl)
library(writexl)
library(janitor)
require(Hmisc)

rtPCRUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(title="Settings", width=3, status="primary", solidHeader=TRUE,
      fileInput(ns("userfile"), "Upload File", accept=c(".xls", ".xlsx")),
      selectInput(ns("target"), "Targets", choices=NULL, multiple=TRUE),
      selectInput(ns("sample"), "Samples", choices=NULL, multiple=TRUE),
      checkboxInput(ns("na.do"), "Auto handle NA", value=TRUE),
      radioButtons(ns("show_data"), "Display Data", 
                   choices=c('ct','d_ct','dd_ct','dd_ct_2n','dd_ct_2n_norm'), 
                   selected='dd_ct_2n_norm'),
      hr(),
      actionButton(ns("apply"), "Analyze", icon=icon("play"), class="btn-primary"),
      downloadButton(ns("download_data"), "Download")
    ),
    box(title="Results", width=9, status="primary",
      tabsetPanel(
        tabPanel("Facet Plot",
          fluidRow(
            column(9, plotOutput(ns("facet_plot"))),
            column(3,
              sliderInput(ns("w_facet"), "Width", 300, 1200, 700, 50),
              sliderInput(ns("h_facet"), "Height", 200, 800, 400, 50),
              sliderInput(ns("sz_facet"), "Label size", 8, 20, 12),
              sliderInput(ns("ang_facet"), "X angle", 0, 90, 45, 15)
            )
          )
        ),
        tabPanel("Combined Plot",
          fluidRow(
            column(9, plotOutput(ns("combine_plot"))),
            column(3,
              sliderInput(ns("w_comb"), "Width", 300, 1200, 700, 50),
              sliderInput(ns("h_comb"), "Height", 200, 800, 400, 50),
              sliderInput(ns("sz_comb"), "Label size", 8, 20, 12),
              sliderInput(ns("ang_comb"), "X angle", 0, 90, 45, 15)
            )
          )
        ),
        tabPanel("Data",
          tableOutput(ns("table_summary")), hr(),
        ),
        tabPanel("Melt Curve",
          plotOutput(ns("plot_mc_1"), width="100%", height="300px"), hr(),
          plotOutput(ns("plot_mc_2"), width="100%", height="400px")
        )
      )
    )
  )
}

rtPCRServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    rv <- reactiveValues(ct=NULL, mc=NULL, result=NULL)
    
    # Read uploaded file
    observeEvent(input$userfile, {
      req(input$userfile)
      
      new_path <- paste0(input$userfile$datapath, ".xlsx")
      file.rename(input$userfile$datapath, new_path)
      # new_path <- "C:\\Users\\haohe\\Desktop\\gy20260121 TEAD CIP P51 P53 TRULI CTGF CYR61_Copy_20260121_100251_Admin_Results_20260121_120831.xlsx"

       # find start row
      for(row_x in 24:50) {
        # print(row_x)
        suppressMessages(ct <- read_excel(new_path, sheet='Results', skip=row_x))
        if(colnames(ct)[1] == 'Well') break
      }
      
      # deal with CT data
      ct <- ct %>% 
        dplyr::select(2,4,5,13) %>% 
        set_names(c("well", "sample", "target", "ct")) %>% 
        janitor::remove_empty('rows') #%>% print()
      rv$ct <- ct
      
      # deal with MC data
      mc_sheet <- ifelse('Melt Curve Raw' %in% excel_sheets(new_path), 'Melt Curve Raw', 
                         ifelse('Melt Curve Raw Data' %in% excel_sheets(new_path), 'Melt Curve Raw Data', NA))
      if(!is.na(mc_sheet)) {
        print('read melting cureve data')
        mc <- read_excel(new_path, sheet=mc_sheet, skip=row_x) %>% 
          dplyr::select(2,4:7) %>% 
          set_names(c("well", "target", "temperature", "fluorescence", "derivative")) %>% 
          janitor::remove_empty('rows') %>% 
          left_join(ct[, 1:3]) #%>% print()
        rv$mc <- mc
      }

      updateSelectInput(session, "sample", choices=unique(ct$sample))
      updateSelectInput(session, "target", choices=unique(ct$target))
    })
    
    
    # Analyze
    observeEvent(input$apply, {
      req(rv$ct, input$sample, input$target)
      
      # filter sample, target
      ct_filtered <- rv$ct %>%
        dplyr::filter(sample %in% input$sample, target %in% input$target)
      # deal with NA
      if (isTRUE(input$na.do)) {
        ct_filtered <- ct_filtered %>% 
          mutate(mean=mean(ct, na.rm=T), .by=c(sample, target)) %>% 
          mutate(ct=ifelse(is.na(ct), rnorm(1, mean, sd=0.5), ct))
      }
      
      # calculate dCT
      ref_target <- ct_filtered %>% 
        # dplyr::filter(target == 'ACTB') %>% 
        dplyr::filter(target == input$target[1]) %>%
        dplyr::distinct(sample, mean) %>% 
        dplyr::rename(ref_target_mean=mean)
      d_ct <- ct_filtered %>% 
        left_join(ref_target) %>% 
        mutate(d_ct=ct - ref_target_mean)
      
      # calculate ddCT
      ref_sample <- d_ct %>% 
        # dplyr::filter(sample == "1") %>% 
        dplyr::filter(sample == input$sample[1]) %>%
        reframe(ref_sample_mean=mean(d_ct, na.rm=T), .by=c(target))
      dd_ct <- d_ct %>% 
        left_join(ref_sample) %>% 
        mutate(dd_ct=d_ct - ref_sample_mean) %>% 
        mutate(dd_ct_2n=2^dd_ct, 
               dd_ct_2n_mean=mean(dd_ct_2n, na.rn=T), .by=c(sample, target)
               )
      
      # norm ref_sample to 1
      ref_norm_factor <- dd_ct %>% 
        # dplyr::filter(sample == "1") %>% 
        dplyr::filter(sample == input$sample[1]) %>%
        reframe(norm_factor=mean(dd_ct_2n_mean, na.rm=T), .by=c(target))
      dd_ct_final <- dd_ct %>% 
        left_join(ref_norm_factor) %>% 
        mutate(dd_ct_2n_norm=dd_ct_2n / norm_factor) %>% 
        dplyr::select(well, sample, target, ct, 
                      d_ct, dd_ct, dd_ct_2n, dd_ct_2n_norm) #%>% print()
      
      # update data
      rv$result <- dd_ct_final

      # melting curve
      if(!is.null(rv$mc)) {
        mc_filtered <- rv$mc %>%
          dplyr::filter(sample %in% input$sample, target %in% input$target)
        # update data
        rv$mc <- mc_filtered
      }
    })
    
    
    # Plot
    output$facet_plot <- renderPlot({
      req(rv$result)
      
      ggplot(rv$result, aes(sample, !!as.name(input$show_data), fill=target)) +
        stat_summary(geom="errorbar", fun.data=mean_sdl, fun.args=list(mult=1),
                     width=0.4, position=position_dodge(0.75)) + 
        stat_summary(geom="bar", fun=mean,
                     width=0.7,  position=position_dodge(0.75)) +
        facet_wrap(vars(target), scales="free_y") + 
        labs(x="") +
        theme_classic() +
        theme(axis.text=element_text(size=input$sz_facet),
              axis.text.x=element_text(angle=input$ang_facet, vjust=0.5),
              legend.position="none")
    })
    
    output$combine_plot <- renderPlot({
      req(rv$result)
      
      ggplot(rv$result, aes(sample, !!as.name(input$show_data), fill=target)) +
        stat_summary(geom="errorbar", fun.data=mean_sdl, fun.args=list(mult=1),
                     width=0.4, position=position_dodge(0.75)) + 
        stat_summary(geom="bar", fun=mean, 
                     width=0.7,  position=position_dodge(0.75)) + 
        labs(x="") +
        theme_classic() +
        theme(axis.text=element_text(size=input$sz_comb),
              axis.text.x=element_text(angle=input$ang_comb, vjust=0.5),
              legend.position="top", legend.title=element_blank())
    })

    output$plot_mc_1 <- renderPlot({
      req(rv$mc)
      
      ggplot(rv$mc, aes(temperature, derivative, color=target)) +
        facet_wrap(vars(target), scales="free_y") + 
        geom_line() +
        theme_minimal() + theme(legend.position="none")
    })
    
    output$plot_mc_2 <- renderPlot({
      req(rv$mc)
      
      ggplot(rv$mc, aes(temperature, derivative, color=target)) +
        facet_wrap(vars(target, sample), scales="free_y") + 
        geom_line() +
        theme_minimal() + theme(legend.position="none")
    })
    
    output$table_summary <- renderTable({
      req(rv$result); rv$result
    }, striped=TRUE, hover=TRUE)
    
    output$download_data <- downloadHandler(
      filename = function() {
        str_split(as.character(Sys.time()), pattern='\\.', simplify=T)[1] %>% 
          str_replace_all(":", "") %>% 
          str_replace_all(" ", "_") %>% 
          str_c("_qpcr_data.xlsx")
      },
      content = function(file) {
        write_xlsx(rv$result, file)
      })
    
  })
}

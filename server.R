
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)
library(knitr)
library(xtable)
source("r/entero.r")


shinyServer(function(input, output) {

  #process data#
  results <- reactive({
    if(input$ent_process == 0)return(NULL)
    isolate({
      processor <- get(input$Assay)
      
      # Process data  
      result <- try(processor(input$data_file$datapath, input$plat, input$org))
      success <- class(result) != "try-error"
      # Store file
      if(input$db & success){
        time <- gsub(" ", "-", Sys.time())
        file.copy(input$data_file$datapath, 
                  paste0("submissions/",
                         input$Assay, "/",
                         input$org, "/",
                         input$date, "_", 
                         time, "_",
                         input$describe, ".csv"))
      }
      result
    })
  })
  
  output$goodresults <- downloadHandler(
    filename = "report.pdf",
    content = function(f){
      knit2pdf(input=paste0("templates/", input$Assay, "/report.Rtex"),
               output = paste0("templates/", input$Assay, "/report.tex"),
               compiler="xelatex")
      file.copy(paste0("templates/", input$Assay, "/report.pdf"),
                f)
    })
  
  output$result <- renderUI({
    if(class(results()) == "try-error")
      renderText(results())
    else if(is.null(results()))
      NULL 
    else
      downloadButton("goodresults", "Download Results")  
  })  
})

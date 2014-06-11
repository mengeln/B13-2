library(shiny)

shinyUI(bootstrapPage(
  titlePanel("Bight '13 qPCR Tool"),
  fluidRow(
    column(1),
    column(4, navlistPanel(tabPanel("Overview", includeMarkdown("resources/overview.Rmd")),
                           tabPanel("Resources", includeMarkdown("resources/resources.Rmd")),
                           tabPanel("Plate Setup", includeMarkdown("resources/plates.Rmd")))),
    column(3,
           radioButtons("Assay", "Choose assay:",
                        c("Entero TaqEnviron" = "ent",
                          "HF183" = "hf")),
           fileInput("data_file", "Select file"),
           conditionalPanel(condition = "input.Assay == 'hf'",
                            fileInput("sketa_file", "Select sketa file")),
           radioButtons("plat", "Select Platform:",
                        c("CFX" = "cfx",
                          "ABI" = "abi")),
           hr(),
           checkboxInput("db", "Submit data to SCCWRP", FALSE),
           conditionalPanel(condition = "input.db == 1",
                            selectInput("org", "Organization",
                                        c("SCCWRP" = "SCCWRP",
                                          "OCSD" = "OCSD")),
                            dateInput("date", "Date of Experiment"),
                            textInput("describe", "Description/PlateID")
           ),
           hr(),  
           actionButton("ent_process", "Process"),
           uiOutput("result"))
    )
))



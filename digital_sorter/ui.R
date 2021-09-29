
library(shiny)
library(shinythemes)


#Import Data
 genelist <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 
               'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z')


# Define UI
ui <-fluidPage(
  
  h1("Digital Cell Sorter"), 
  p(style = "font-family:Impact"),
  tags$hr(),
  
  
  sidebarLayout(
    
    
    sidebarPanel(
      titlePanel("Input"),
      #width = 3,
      
      
      #gene
      selectizeInput(inputId = "gene", label= "Select genes",
                     choices = genelist,
                     selected = NULL,
                     multiple = TRUE,
                     options = list(
                       placeholder = 'Type to search for gene',
                       onInitialize = I('function() { this.setValue(""); }')
                     )
      ),
     
      #cancer-lung
      checkboxGroupInput(inputId ="cancer", label= "Select cancer types", choices = "Lung", selected = "Lung"),
      #markers
      conditionalPanel('input.cancer == "Lung"', 
                       selectizeInput(inputId = "marker", 
                                      label= "master markers",
                                      multiple = F, 
                                      choices = c("PECAM1","EPCAM","PTPRC"),
                                      selected = c("PTPRC"))
                       ),
                     
      
      selectizeInput(
        inputId = "cohort", label= "Select dataset of interest", choices = c("song_2019","nsclc_primary"), 
        selected = NULL), 
      
      uiOutput("sample"),
                       
      actionButton("go", "Go")
      
    ),#side panel end
    
    mainPanel(
      navbarPage(title = "Results",
                 windowTitle = "Digital Cell Sorter",
                 theme = shinytheme("cerulean"),
                 collapsible = TRUE,
                 tabPanel(icon("home"),
                          
                          fluidRow(
                            textOutput("result1"),textOutput("result2"),
                            tags$hr(),
                            print("You chose gene(s) "),
                            textOutput("result3")
                                  )
                          ),
                 tabPanel("Violin Plot",
                          plotlyOutput("plotv")
                          ),
                 tabPanel("Dot Plot",
                 ),
                 tabPanel("tSNEs",
                 ),
                 tabPanel("Marker Dot Plot (FACS)",
                 )
      )
      
      
    )#main panel end
     
  )#side bar layout end
      
) #page end

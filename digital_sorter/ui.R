
library(shiny)
library(shinythemes)



#Import Data
genelist <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 
              'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z')


# Define UI
ui <-fluidPage(
  useShinyjs(),
  
  h1("Digital Cell Sorter"), 
  p(style = "font-family:Impact"),
  tags$hr(),
  
  
  sidebarLayout(
    
    
    sidebarPanel(
      titlePanel("Input"),
      #width = 3,
      
      
  
      #cancer
      selectizeInput(inputId ="cancer", label= "Select cancer types:", choices = c("Lung","t.b.c."), selected = "Lung"),
      
      #show master markers after selecting cancer type
      p(strong("Masters markers :")) ,                
      verbatimTextOutput("master_markers"),
      
      
      #select dataset (choices depends on cancer)                 
      uiOutput("cohort"),
      
      
      actionButton("go", "Draw main violin plots!"),
      
      tags$hr(),
      
      #gene
      uiOutput("sample_se"),
      
      uiOutput("gene"),
      
      
      checkboxInput(inputId= "testDGE", label="Automatically select differentially expressed genes", value = FALSE, width = NULL),
      actionButton("levelplot", "Draw plots for each level!")
      
    ),#side panel end
    
    mainPanel(
      navbarPage(title = "Results",
                 windowTitle = "Digital Cell Sorter",
                 theme = shinytheme("cerulean"),
                 collapsible = TRUE,
                 tabPanel(icon("home"),
                          
                          fluidRow(
                            column(width=8, offset=3, #Tree plot
                                   conditionalPanel(
                                     condition = "input.cancer == 'Lung'",
                                     img(src = "tree.PNG", width=680, height=500)  
                                   ),
                                   conditionalPanel(
                                     condition = "input.cancer == 't.b.c.'",
                                     h4("Tree plot t.b.c.")
                                   )
                            ) ),
                          fluidRow(
                            column(width=12,#violin plots
                                   h4("Level 1: CD45 (PTPRC)"),
                                   plotlyOutput("plotv_h1"),
                                   h4("Level 2: EPCAM"),
                                   plotlyOutput("plotv_h2"),
                                   h4("Level 3: PECAM1 (CD31)"),
                                   plotlyOutput("plotv_h3")
                            )
                          )
                          
                 ), #home end
                 tabPanel("Level 1",
                           fluidRow( #violin plot
                            column(width=12, plotlyOutput("plotv1"))
                            
                          ),
                          
                          fluidRow(#dot plot
                            tags$hr(),
                            h4(textOutput("dotplot_title1")),
                            column(width=8, style = "height:500px",
                                   plotOutput("plotd1_1")),
                            column(width=4,
                                  # dataTableOutput("table1_1")
                                  )
                            
                          ),
                         
                          fluidRow(#dot plot
                            h4(textOutput("dotplot_title2")),
                            column(width=8, style = "height:500px",
                                  plotOutput("plotd1_2")
                                  ),
                            column(width=4,
                                   #dataTableOutput("table1_2")
                                   )
                            
                          )
                 ),
                 tabPanel("Level 2",
                          fluidRow( #violin plot
                            plotlyOutput("plotv2")
                          ),
                          
                          fluidRow( #dot plot
                            # plotlyOutput("plotd1")
                          )
                 ),
                 tabPanel("Level 3",
                          fluidRow( #violin plot
                            plotlyOutput("plotv3")
                          ),
                          
                          fluidRow( #dot plot
                            # plotlyOutput("plotd1")
                          )
                 ),
                 tabPanel("Level 4",
                 )
      )
      
      
    )#main panel end
    
  )#side bar layout end
  
) #page end
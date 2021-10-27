
library(shiny)
library(shinythemes)
library(shinydashboard)
library(dashboardthemes)


# Define UI
ui <- dashboardPage(skin = "blue",
                   
                  
                    #1. header
                    dashboardHeader(
                      title = HTML("Digital Cell Sorter"), 
                      disable = FALSE, 
                      titleWidth  = 300,
                      dropdownMenu(type = "notifications", badgeStatus = NULL, icon = icon('comment'), headerText = NULL,
                                   messageItem("Feedback and suggestions",
                                               "",
                                               time = NULL,
                                               icon = icon("envelope"),
                                               href = "mailto:c.lee@dkfz-heidelberg.de")
                                   )
                    ), #header end
                    
                    #2. side bar
                    dashboardSidebar(
                      width = 300,
                      sidebarMenu(
                        id = 'sidebar',
                        ## 1st tab show the Main dashboard -----------
                        menuItem( "Main Dashboard", tabName = 'dashboard', icon = icon('tachometer-alt')),
                        
                        useShinyjs(),
                        
                        
                                         
                        ## 2nd tab shows results ----------
                        menuItem( "Results", tabName = "results", icon = icon('barcode'), startExpanded = F,
                                  h5("Input section for results"),
                                  uiOutput("sample_se"),
                                  uiOutput("gene"),
                                  checkboxInput(inputId= "testDGE", 
                                                label="Automatically select differentially expressed genes", 
                                                value = FALSE, width = NULL),
                                  actionBttn("levelplot", "Plots for each level! (See Results sections)", 
                                             style = "jelly", color = "success", size = "xs"),
                                  menuSubItem('Level 1', tabName = "result1"),
                                  menuSubItem('Level 2', tabName = "result2"),
                                  menuSubItem('Level 3', tabName = "result3"),
                                  menuSubItem('Level 4', tabName = "result4")
                                  
                        ),
                        ## 3rd tab Data source, definition , i.e., help ---------------
                        menuItem( "FAQs", tabName = 'help', icon = icon('question-circle') )
                      )#sidebarMenu end
                    ), #side bar end 
                    
                    dashboardBody(
                      
                      mainPanel(
                        ## 3.1 Dashboard body --------------
                        tabItems(
                          ## 3.1.1 Main dashboard ----------------------------------------------------------
                          tabItem( tabName = 'dashboard',
                                   
                                   fluidRow(
                                     box(title = "Inputs", status= "info", solidHeader = TRUE, width=3,
                                            
                                            #cancer
                                            selectizeInput(inputId ="cancer", label= "Select cancer types:", choices = c("Lung","t.b.c."), selected = "Lung"),
                                            
                                            #show master markers after selecting cancer type
                                            p(strong("Masters markers :")) ,                
                                            verbatimTextOutput("master_markers"),
                                            
                                            
                                            #select dataset (choices depends on cancer)                 
                                            uiOutput("cohort"),
                                            
                                            
                                            actionBttn("go", "Draw main violin plots!", style = "jelly", color = "primary",size = "sm")
                                     ),
                                     box(title = "Tree plot", status="warning", solidHeader = TRUE, width=9,
                                           
                                            conditionalPanel(
                                              condition = "input.cancer == 'Lung'",
                                              img(src = "tree.PNG", width=680, height=500)  
                                            ),
                                            conditionalPanel(
                                              condition = "input.cancer == 't.b.c.'",
                                              h4("Tree plot t.b.c.")
                                            )
                                     )
                                   ), #fluid row end
                                   fluidRow(
                                     tabBox(title = "Violin plots", 
                                            id="vplots", 
                                            height = "260px", 
                                            width=12, side= "left",
                                            tabPanel("Level 1: CD45 (PTPRC)", 
                                                     plotlyOutput("plotv_h1")),
                                            tabPanel("Level 2: EPCAM", 
                                                     plotlyOutput("plotv_h2")),
                                            tabPanel("Level 3: PECAM1 (CD31)", 
                                                    plotlyOutput("plotv_h3"))
                                     )
                                    
                                   )#fluid row end
                          ), #tab1 end
                          
                          ## 3.1.2 Result1 ----------------------------------------------------------
                          
                          tabItem( tabName = 'result1',
                                   
                                   fluidRow(
                                     box(title = "Violin plot", status="warning", solidHeader = TRUE, width=12,
                                     plotlyOutput("plotv1")
                                     )
                                   ),#fluid flow end
                                   fluidRow(
                                     h4(textOutput("dotplot_title1")),
                                     box(title = "Dot plot", status="warning", solidHeader = TRUE, width=8,
                                         plotlyOutput("plotd1")
                                     ),
                                     box(title = "Data table of dot plot", status="warning", solidHeader = TRUE, width=4,
                                         h4(textOutput("dotplot_title1")),
                                         plotlyOutput("table1")
                                     )
                                   )#fluid row end
                          ), #tab2 end
                          
                          ## 3.1.3 Result2 ----------------------------------------------------------
                          
                          tabItem( tabName = 'result2',
                                   
                                   
                                   fluidRow(
                                     plotlyOutput("plotv2")
                                   ),#fluid row end
                                   fluidRow(#dot plot
                                     
                                   )#fluid flow end
                          ), #tab3 end
                          
                          ## 3.1.4 Result3 ----------------------------------------------------------
                          
                          tabItem( tabName = 'result3',
                                   fluidRow(
                                     plotlyOutput("plotv3")
                                   ),#fluid flow end
                                   
                                   fluidRow(#dot plot
                                     
                                     
                                   )#fluid row end
                          ), #tab4 end
                          
                          ## 3.1.5 Result4 ----------------------------------------------------------
                          
                          tabItem( tabName = 'result4',
                                   
                                   fluidRow(
                                     
                                   ),#fluid row end
                                   fluidRow(
                                     
                                   )#fluid flow end
                          ) #tab5 end
                          
                          
                        )#tabItems end
                      )    #main Panel
                  )      #dashboard body
   )  #UI

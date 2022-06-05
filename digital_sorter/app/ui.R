library(shiny)
library(shinybusy)
library(shinythemes)
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(plotly)
library(gridExtra)
library(shinyjs)
library(shinydashboard)
library(shinyWidgets)

# Define UI
ui <-dashboardPage(
  
  skin = "blue",
  title= "Digital Sorter",
  
  
  #1. header------------
  dashboardHeader(
    
    
    title = tags$b(tags$img(src = "digital_sorter.png", width=80, height=50),
                   " Digital Sorter"), 
    
    disable = FALSE, 
    titleWidth  = 300
  ), #header end
  
  
  #2. side bar-----------
  dashboardSidebar(
    width = 300,
    useShinyjs(),
    
    sidebarMenu(
      id = 'sidebar',
      ## 1st tab show the pre-processing dashboard  
      menuItem( "Define My Own Master Markers", tabName = 'processing', icon = icon('tachometer-alt')),
      
      ## 2st tab show the Tree plot and violin plot
      menuItem( "Tree Plot & Violin Plot", tabName = 'dashboard', icon = icon('tree')),
      
      ## 3nd tab shows visualization results (dot plots of gene of interest) 
      menuItem( "Visualization (gene of interest)", tabName = "results", icon = icon('list'), startExpanded = F,
                
                
                h4("Input section for each level"),
                #select dataset (choices depends on cancer)                 
                uiOutput("cohort1"),
                
                uiOutput("sample_se1"),
              
                uiOutput("gene"),
                
                
                h5("Go to each level and check out the results!"),
                lapply(1:5, function(i){
                  menuSubItem(paste0("Level ",i), tabName = paste0("result",i))
                })
              
                
                
                
      ),
      
      ## 4th tab shows visualization results (dot plots of results of marker selection) 
      menuItem( "Visualization (marker selection)", tabName = "results2", icon = icon('brain'), startExpanded = F,
                
               
                h4("Input section for each level"),
                #select dataset (choices depends on cancer)                 
                uiOutput("cohort2"),
                
                uiOutput("sample_se2"),
                
                h6("Warning! Marker selection function takes quite some time!"),
                actionBttn("dplot2", "Do marker selection", style = "jelly", color = "warning",size = "sm"),
                
                
                fluidRow(
                  column(width=4,offset=1,
                         ## let user download the marker selection table
                         downloadButton('download_marker_csv', 'Download Marker Selection Results')
                  ),
                ),
                uiOutput("celltype"),
                
                
                h5("Go to each level and check out the results!"),
                lapply(1:5, function(i){
                  menuSubItem(paste0("Level ",i), tabName = paste0("2result",i))
                })
        
                
                
                
      ),
     
      
      
      ## 5th tab Data source, definition , i.e., help 
      menuItem( "FAQs", tabName = 'help', icon = icon('question-circle') )
    )#sidebarMenu end
    
  ), #side bar end
  
  # 3. Dashboard body --------------
  dashboardBody( 
    
    tabItems(
      ## 3.1 Pre-processing dashboard ----------------------------------------------------------
      tabItem( tabName = 'processing',
               fluidRow(
                 column(width=4,
                        box(title = "Datasets stratified with default master markers", status= "success", solidHeader = TRUE, width = NULL,
                            #cancer
                            selectizeInput(inputId ="cancer", label= "Select cancer types:", choices = c("Lung","t.b.c."), selected = "Lung"),
                            #show master markers after selecting cancer type
                            uiOutput("master_markers"),
                            
                            
                            actionBttn("go", "Confirm & Go!", style = "jelly", color = "success",size = "sm")
                        ),# box end
                        box(title = "Define your own master markers", status= "warning", solidHeader = TRUE, width = NULL,
                            uiOutput("own_markers"),
                            textInput(inputId="own_markers_name", label="Name your own master markers (e.g. Lung1)", value = "", width = NULL, placeholder = NULL),
                            actionBttn("update_marker", "Update marker list and tree plot!", style = "jelly", color = "warning",size = "sm"),
                            h6(span("Function not yet finished", style = "color:orange")),
                            h5("This is the panel for user to define their own master markers."),
                            h5("If the cancer type is not",span(" Lung", style = "color:blue") ,", user should upload the reference expression data to define the levels.")
                        ),# box end
                       
                        ),#column end
                 column(width=8,
                        box(title = "The structure of the embedded data", status="success", solidHeader = TRUE, width = NULL,
                            h5("The results shown in this app will be based on the cell groups/levels shown here."),
                            
                            h5("If you want to change the stratifying structure, please use the",
                            span(" yellow ", style = "color:orange") , "panel. Otherwise, you can go to the ",
                            span("Tree Plot & Violin Plot", style = "color:green")," for the following visualization."),
                            
                            #Tree plot
                            conditionalPanel(
                              condition = "input.cancer == 'Lung'",
                              img(src = "tree.PNG", width=550, height=480)  
                            ),
                            conditionalPanel(
                              condition = "input.cancer == 't.b.c.'",
                              h4("Tree plot t.b.c.")
                              
                            )
                        ),#box end
                 ),#column end
               ), #fluidRow end
              
               
                
            
               
      ), #tab1 end
      
      ## 3.2 Tree plot and violin plot ----------------------------------------------------------
      tabItem( tabName = 'dashboard',
               
               fluidRow(
                 
                 box(title = "Updating Tree Plot", status="warning", solidHeader = TRUE, width=8,
                        h5("The results shown in this app were the analysis done within the cell groups shown here."),
                        #Tree plot
                        h3("Should be the updating tree plot")
                     ##automatic tree plot
                     ######################
                  )
                 ), #fluid row end
               
               fluidRow(
                 tabBox(title = "Violin plots of all cell types", 
                        id="vplots", 
                        width=12, side= "left",
                        tabPanel("Description", 
                                 #select dataset (choices depends on cancer)                 
                                 uiOutput("cohort_h"),
                                 
                                 h4("Showing cancer type:",textOutput("cancer_type")),
                                
                                 
                                 h5("If you want to change the cancer type, go to panel ", span("Define My Own Master Markers.", style = "color:orange") ),
                                 h4("Otherwise, press",span(" Plot! ", style = "color:green")),
                                 actionBttn("plotv", "Plot!", style = "jelly", color = "success",size = "sm")
                                 
                                
                                 ),
                        tabPanel("Marker 1", 
                                 withSpinner(plotOutput("plotv_h1"))),
                        tabPanel("Marker 2", 
                                 withSpinner(plotOutput("plotv_h2") )),
                        tabPanel("Marker 3", 
                                 withSpinner(plotOutput("plotv_h3") )),
                        tabPanel("Marker 4", 
                                 withSpinner(plotOutput("plotv_h4") ))
                        
                        
                 )
               )#fluid row end
      ), #tab end
      
      ## 3.3 Result (gene of interest)----------------------------------------------------------
       ### 3.3.1 Level 1 (gene of interest)----------------------------------------------------------
 
      tabItem( tabName = 'result1',
      
               
               fluidRow(#dot plot
      
                 h3("Level 1 Dot plot with genes of interest"),
                 h4(textOutput("dotplot_title1")),
                 column(width=6, style = "height:200px;",
                        withSpinner(plotOutput("plotd1") )
                        
              
                 ),
                 column(width=4,
                        withSpinner(dataTableOutput("table1") )
                        
                 )
                 
               )#fluid row end
            ), #tab end
      
      ### 3.3.2 Level 2 (gene of interest)----------------------------------------------------------
      
      tabItem( tabName = 'result2',
              
               fluidRow(#dot plot
                 h3("Level 2 Dot plot with genes of interest"),
                 h4(textOutput("dotplot_title2")),
                 column(width=6, style = "height:200px;",
                        withSpinner(plotOutput("plotd2") )
                        
                   ),
                 column(width=4,
                        withSpinner(dataTableOutput("table2") )
                        
                 )
                 
               )#fluid row end
      ), #tab end
      
      ### 3.3.3 Level 3 (gene of interest)----------------------------------------------------------
      
      tabItem( tabName = 'result3',
              
               
               fluidRow(#dot plot
                 h3("Level 3 Dot plot with genes of interest"),
                 h4(textOutput("dotplot_title3")),
                 
                 column(width=6, style = "height:200px;",
                        
                        withSpinner( plotOutput("plotd3") )
                
                 ),
                 column(width=4,
                        withSpinner(dataTableOutput("table3") )
                        
                 )
                        
                 )#fluid row end
      ), #tab end
      
      ### 3.3.4 Level 4 (gene of interest)----------------------------------------------------------
      
      tabItem( tabName = 'result4',
               
               
               fluidRow(#dot plot
                 h3("Level 4 Dot plot with genes of interest"),
                 h4(textOutput("dotplot_title4")),
                 
                 column(width=6, style = "height:200px;",
                        withSpinner(plotOutput("plotd4")  )
                        
                 ),
                 column(width=4,
                        withSpinner(dataTableOutput("table4") )
                        
                
                 )      
                 )#fluid row end
      ), #tab end
      
      ### 3.3.5 Level 5 (gene of interest)----------------------------------------------------------
      tabItem( tabName = 'result5',
               
              
             fluidRow(#dot plot
               h3("Level 5 Dot plot with genes of interest"),
                 h4(textOutput("dotplot_title5")),
                 
                 column(width=6, style = "height:200px;",
                        withSpinner(plotOutput("plotd5") )
                        
                 ),
                 column(width=4,
                        withSpinner(dataTableOutput("table5") )
                       
                        
                 )      
               )#fluid row end
      ), #tab end
      
      ## 3.4 Result2 (marker selection)----------------------------------------------------------
      ### 3.4.1 Level 1 (marker selection)----------------------------------------------------------
      tabItem( tabName = '2result1',
               
               fluidRow(#dot plot
                 h3("Level 1 Dot plot with specific cell surface markers"),
                 h4(textOutput("dotplot2_title1")),
                 column(width=6, style = "height:200px;",
                        withSpinner(plotOutput("plotd_cellMark1") )
                        
                        
                 ),
                 column(width=4,
                        withSpinner(dataTableOutput("table_cellMark1") )
                        
                 )
                 
               )#fluid row end
      ), #tab end
      
      ### 3.4.2 Level 2 (marker selection)----------------------------------------------------------
      
      tabItem( tabName = '2result2',
              
               fluidRow(#dot plot
                h3("Level 2 Dot plot with specific cell surface markers"),
                 h4(textOutput("dotplot2_title2")),
                 column(width=6, style = "height:200px;",
                        withSpinner(plotOutput("plotd_cellMark2") )
                        
                 ),
                 column(width=4,
                        withSpinner(dataTableOutput("table_cellMark2") )
                        
                 )
                 
               )#fluid row end
      ), #tab end
      
      ### 3.4.3 Level 3 (marker selection)----------------------------------------------------------
      
      tabItem( tabName = '2result3',
               
               fluidRow(#dot plot
                 
                 h3("Level 3 Dot plot with specific cell surface markers"),
                 h4(textOutput("dotplot2_title3")),
                 
                 column(width=6, style = "height:200px;",
                        
                        withSpinner( plotOutput("plotd_cellMark3") )
                        
                 ),
                 column(width=4,
                        withSpinner(dataTableOutput("table_cellMark3") )
                        
                 )
                 
               )#fluid row end
      ), #tab end
      
      ### 3.4.4 Level 4 (marker selection)----------------------------------------------------------
      
      tabItem( tabName = '2result4',
               
               fluidRow(#dot plot
                 
                 h3("Level 4 Dot plot with specific cell surface markers"),
                 h4(textOutput("dotplot2_title4")),
                 
                 column(width=6, style = "height:200px;",
                        withSpinner(plotOutput("plotd_cellMark4")  )
                        
                 ),
                 column(width=4,
                        withSpinner(dataTableOutput("table_cellMark4") )
                        
                        
                 )      
               )#fluid row end
      ), #tab end
      ### 3.4.5 Level 5 (marker selection)----------------------------------------------------------
      tabItem( tabName = '2result5',
               
            
               fluidRow(#dot plot
                
                 h3("Level 5 Dot plot with specific cell surface markers"),
                 h4(textOutput("dotplot2_title5")),
                 
                 column(width=6, style = "height:200px;",
                        withSpinner(plotOutput("plotd_cellMark5") )
                        
                 ),
                 column(width=4,
                        withSpinner(dataTableOutput("table_cellMark5") )
                        
                        
                 )      
               )#fluid row end
      ), #tab end
      
      ## 3.5 FAQs ----------------------------------------------------------
      
      tabItem( tabName = 'help',
               
               fluidRow(
                 column( width=11.5,offset = 0.5,
                         h3("Tutorial and FAQs"),
                         h4("Summary of the integrated dataset"),
                         p("The integrated dataset imbedded in this shiny app consists of 7 datasets which were extracted from the publications from 2019 to 2021", style = "font-family: 'times'; font-si16pt"),
                         strong("There are:"),
                         p(em("lambrechts_2018"),"includes",span("adjacent normal, LUAD, LUSC,", style = "color:blue"),"and",span("NSCLC primary", style = "color:blue"),"sample types."),
                         p(em("song_2019"),"includes",span("adjacent normal", style = "color:blue"),"and",span("NSCLC primary", style = "color:blue"),"sample types."),
                         p(em("zilionis_2019"),"includes",span("blood", style = "color:blue"),"and",span("NSCLC primary", style = "color:blue"),"sample types."),
                         p(em("kim_2020"),"includes",span("NSCLC primary, NSCLC meta, pleural fluids, lymphNode normal, lymphNode meta", style = "color:blue"),"and",span("adjacent normal", style = "color:blue"),"sample types."),
                         p(em("travaglini_2020"),"includes only",span("general normal", style = "color:blue"),"sample types."),
                         p(em("bischoff_2021"),"includes",span("LUAD", style = "color:blue"),"and",span("adjacent normal", style = "color:blue"),"sample types."),
                         p(em("wu_2021"),"includes",span("LUAD and LUSC", style = "color:blue"),"and",span("NSCLC primary", style = "color:blue"),"sample types.")
                         
                         )
               
                
               )
      ) #tab end
    )#tabItems end
  )#dashboard body end
) #page end

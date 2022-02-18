
library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinycssloaders)


# Define UI
ui <-dashboardPage(
  
  skin = "blue",
  
  #1. header
  dashboardHeader(
    title = tags$b(tags$img(src = "digital_sorter.png", width=80, height=50),
                   " Digital Sorter"), 
    
    disable = FALSE, 
    titleWidth  = 300,
    dropdownMenu( type = "notifications", badgeStatus = NULL, icon = icon('comment'), headerText = NULL,
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
    useShinyjs(),
    sidebarMenu(
      id = 'sidebar',
      ## 1st tab show the Main dashboard -----------
      menuItem( "Main Dashboard", tabName = 'dashboard', icon = icon('tachometer-alt')),
      ## 2nd tab shows results ----------
      menuItem( "Results", tabName = "results", icon = icon('barcode'), startExpanded = F,
                
                lapply(1:5, function(i){
                  menuSubItem(paste0("Level ",i), tabName = paste0("result",i))
                }),
                #uiOutput("ui_menuSubItem"), #format would be different
                
                h4("Input section for each level"),
                #select dataset (choices depends on cancer)                 
                uiOutput("cohort2"),
                
                uiOutput("sample_se"),
              
                uiOutput("gene"),
                
                checkboxInput(inputId= "cellMark", 
                              label="Automatically select top markers of specific cell type!", 
                              value = FALSE, width = NULL),
                downloadButton('download_marker_csv', 'Download Marker Selection Results'),
                uiOutput("celltype"),
                
                actionBttn("dplot", "Plot/ Renew plots!", 
                           style = "jelly", color = "warning", size = "sm"), 
                h5("Go to each level and check out the results!")
                
                
      ),
     
      
      
      ## 3rd tab Data source, definition , i.e., help ---------------
      menuItem( "FAQs", tabName = 'help', icon = icon('question-circle') )
    )#sidebarMenu end
  ), #side bar end
  
  #3. body
  dashboardBody( 
    
    ## 3.1 Dashboard body --------------
    tabItems(
      ## 3.1.1 Main dashboard ----------------------------------------------------------
      tabItem( tabName = 'dashboard',
               
               fluidRow(
                 box(title = "Inputs", status= "success", solidHeader = TRUE, width=4,
                        #cancer
                        selectizeInput(inputId ="cancer", label= "Select cancer types:", choices = c("Lung","t.b.c."), selected = "Lung"),
                        
                        #show master markers after selecting cancer type
                        uiOutput("master_markers"),
                        
                        
                        #select dataset (choices depends on cancer)                 
                        uiOutput("cohort"),
                        
                        
                        actionBttn("go", "Select & Go!", style = "jelly", color = "success",size = "sm")
                  ),
                 box(title = "Tree plot", status="warning", solidHeader = TRUE, width=8,
                        
                        #Tree plot
                        conditionalPanel(
                          condition = "input.cancer == 'Lung'",
                          img(src = "tree.PNG", width=550, height=480)  
                        ),
                        conditionalPanel(
                          condition = "input.cancer == 't.b.c.'",
                          h4("Tree plot t.b.c.")
                        )
                  )
                 ), #fluid row end
               
               fluidRow(
                 tabBox(title = "Violin plots of all cell types", 
                        id="vplots", 
                        width=12, side= "left",
                        tabPanel("Marker 1", 
                                 withSpinner(plotOutput("plotv_h1"))),
                        tabPanel("Marker 2", 
                                 withSpinner(plotOutput("plotv_h2") )),
                        tabPanel("Marker 3", 
                                 withSpinner(plotOutput("plotv_h3") )),
                        tabPanel("Marker 4", 
                                 withSpinner(plotOutput("plotv_h4") ))
                        # error will appear
                        #lapply(1:4, function(i){ 
                        #  tabPanel(paste0("Marker ",i), 
                        #           withSpinner(plotOutput(paste0("plotv_h",i))))
                        #})
                        
                 )
               )#fluid row end
      ), #tab1 end
      
      ## 3.1.2 Result1 ----------------------------------------------------------
 
      tabItem( tabName = 'result1',
      
               fluidRow(#violin plot
                 h3("Level 1 Violin plot"),
                 tags$hr(),
                 withSpinner(plotlyOutput("plotv1") )
                        
                 ),#fluid row end
               fluidRow(#dot plot
                 tags$hr(),
                 h3("Level 1 Dot plot"),
                 h4(textOutput("dotplot_title1")),
                 column(width=6, style = "height:200px;",
                        withSpinner(plotOutput("plotd1") )
                        
              
                 ),
                 column(width=4,
                        withSpinner(dataTableOutput("table1") )
                        
                 )
                 
               )#fluid row end
            ), #tab2 end
      
      ## 3.1.3 Result2 ----------------------------------------------------------
      
      tabItem( tabName = 'result2',
               
               
               fluidRow(
                 h3("Level 2 Violin plot"),
                 withSpinner(plotlyOutput("plotv2") ) 
                 
                 
               ),#fluid row end
               fluidRow(#dot plot
                 tags$hr(),
                 h3("Level 2 Dot plot"),
                 h4(textOutput("dotplot_title2")),
                 column(width=6, style = "height:200px;",
                        withSpinner(plotOutput("plotd2") )
                        
                   ),
                 column(width=4,
                        withSpinner(dataTableOutput("table2") )
                        
                 )
                 
               )#fluid row end
      ), #tab3 end
      
      ## 3.1.4 Result3 ----------------------------------------------------------
      
      tabItem( tabName = 'result3',
               fluidRow(
                 h3("Level 3 Violin plot"),
                 withSpinner(plotlyOutput("plotv3") ) 
                 
                 
               ),#fluid row end
               
               fluidRow(#dot plot
                 tags$hr(),
                 h3("Level 3 Dot plot"),
                 h4(textOutput("dotplot_title3")),
                 
                 column(width=6, style = "height:200px;",
                        
                        withSpinner( plotOutput("plotd3") )
                
                 ),
                 column(width=4,
                        withSpinner(dataTableOutput("table3") )
                        
                 )
                        
                 )#fluid row end
      ), #tab4 end
      
      ## 3.1.5 Result4 ----------------------------------------------------------
      
      tabItem( tabName = 'result4',
               
               fluidRow(
                 h3("Level 4 Violin plot"),
                 withSpinner(plotlyOutput("plotv4") )
                 
          
               ),#fluid row end
               fluidRow(#dot plot
                 tags$hr(),
                 h3("Level 4 Dot plot"),
                 h4(textOutput("dotplot_title4")),
                 
                 column(width=6, style = "height:200px;",
                        withSpinner(plotOutput("plotd4")  )
                        
                 ),
                 column(width=4,
                        withSpinner(dataTableOutput("table4") )
                        
                
                 )      
                 )#fluid row end
      ), #tab5 end
      ## 3.1.6 Result5 ----------------------------------------------------------
      tabItem( tabName = 'result5',
               
               fluidRow(style = "height:200px;",
                        h3("Level 5 Violin plot"),
                        tags$hr(),
                        h3("Level 5 Violin plot is the same as Level 4"),
                        #plotlyOutput("plotv4")
                        
               ),#fluid row end
               fluidRow(#dot plot
                 tags$hr(),
                 h3("Level 5 Dot plot"),
                 h4(textOutput("dotplot_title5")),
                 
                 column(width=6, style = "height:200px;",
                        withSpinner(plotOutput("plotd5") )
                        
                 ),
                 column(width=4,
                        withSpinner(dataTableOutput("table5") )
                       
                        
                 )      
               )#fluid row end
      ), #tab5 end
      ## 3.1.6 FAQs ----------------------------------------------------------
      
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
      ) #tab6 end
    )#tabItems end
  )#dashboard body end
) #page end

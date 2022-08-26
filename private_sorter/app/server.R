library(shiny)
library(shinybusy)
library(shinythemes)
library(shinyjs)
library(shinydashboard)
library(shinyWidgets)
library(shinycssloaders)

library(Seurat)
#library(biomaRt)
library(dplyr)
library(reshape2)
library(ggplot2)
library(plotly)
library(gridExtra)
library(RColorBrewer)

#library(shinyauthr)
library(DT)
library(sodium)
# Main login screen
loginpage <- div(id = "loginpage", style = "width: 500px; max-width: 100%; margin: 0 auto; padding: 20px;",
                 wellPanel(
                   tags$h2("LOG IN", class = "text-center", style = "padding-top: 0;color:#333; font-weight:600;"),
                   textInput("userName", placeholder="Username", label = tagList(icon("user"), "Username")),
                   passwordInput("passwd", placeholder="Password", label = tagList(icon("unlock-alt"), "Password")),
                   br(),
                   div(
                     style = "text-align: center;",
                     actionButton("login", "SIGN IN", style = "color: white; background-color:#3c8dbc;
                                 padding: 10px 15px; width: 150px; cursor: pointer;
                                 font-size: 18px; font-weight: 600;"),
                     shinyjs::hidden(
                       div(id = "nomatch",
                           tags$p("Oops! Incorrect username or password!",
                                  style = "color: red; font-weight: 600; 
                                            padding-top: 5px;font-size:16px;", 
                                  class = "text-center")))
                   ))
)

# dataframe that holds usernames, passwords and other user data
credentials = data.frame(
  username_id = c("Lung", "Jessie"),
  passod   = sapply(c("2022", "Jessie"),password_store),
  permission  = c("basic", "advanced"), 
  stringsAsFactors = F
)
ReadRDSFiles("data")
#server.R
server <- function(input, output, session){
  
  login = FALSE
  USER <- reactiveValues(login = login)
  
  observe({ 
    if (USER$login == FALSE) {
      if (!is.null(input$login)) {
        if (input$login > 0) {
          Username <- isolate(input$userName)
          Password <- isolate(input$passwd)
          if(length(which(credentials$username_id==Username))==1) { 
            pasmatch  <- credentials["passod"][which(credentials$username_id==Username),]
            pasverify <- password_verify(pasmatch, Password)
            if(pasverify) {
              USER$login <- TRUE
            } else {
              shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade")
              shinyjs::delay(3000, shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade"))
            }
          } else {
            shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade")
            shinyjs::delay(3000, shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade"))
          }
        } 
      }
    }    
  })
  
  output$logoutbtn <- renderUI({
    req(USER$login)
    tags$li(a(icon("door-open"), "Logout", 
              href="javascript:window.location.reload(true)"),
            class = "dropdown", 
            style = "background-color: #eee !important; border: 0;
                    font-weight: bold; margin:5px; padding: 5px;")
  })
  
  output$menu <- renderUI({
    sidebarMenu(
      id = 'sidebar',
      ## 1st tab show the pre-processing dashboard  
      menuItem( "Define My Own Master Markers", tabName = 'processing', icon = icon('tachometer-alt')),
      
      ## 2st tab show the Tree plot and violin plot
      menuItem( "Tree Plot & Violin Plot", tabName = 'dashboard', icon = icon('tree')),
      
      ## 3nd tab shows visualization results (dot plots of gene of interest) 
      menuItem( "Visualization (gene of interest)", expandedName = "results", icon = icon('list'), startExpanded = F,
                
                
                h4("Input section for each level"),
                #select dataset (choices depends on cancer)                 
                uiOutput("cohort1"),
                
                uiOutput("sample_se1"),
                
                uiOutput("gene"),
                
                
                lapply(1:5, function(i){
                  menuSubItem(paste0("Level ",i), tabName = paste0("result",i))
                })
                
                
                
                
      ),
      
      ## 4th tab shows visualization results (dot plots of results of marker selection) 
      menuItem( "Visualization (marker selection)", expandedName = "2results", icon = icon('brain'), startExpanded = F,
                
                
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
                  )
                ),#fluidflow
                uiOutput("celltype"),
                sliderInput("n_genes", "Number of genes on the dotplot",
                            min = 5, max = 20,
                            value = 10),
                
                lapply(1:5, function(i){
                  menuSubItem(paste0("Level ",i), tabName = paste0("2result",i))
                })
                
                
      ),
      
      
      
      ## 5th tab Data source, definition , i.e., help 
      menuItem( "FAQs", tabName = 'help', icon = icon('question-circle') )
    )#sidebarMenu end
  })
    
  output$body <- renderUI({
    if (USER$login == TRUE ) {  
      
      ## 3.1 Dashboard body --------------
      tabItems(
        ## 3.1 Pre-processing dashboard ----------------------------------------------------------
        tabItem( tabName = 'processing',
                 fluidRow(
                   column(width=4,
                          box(title = "Datasets stratified with default master markers", status= "success", solidHeader = TRUE, width = NULL,
                              #genes list
                              selectizeInput(inputId ="genelists", 
                                             label= "Select type of gene list:", 
                                             choices = c("Protein coding genes","Cell surface markers"), 
                                             selected = "Protein coding genes"),
                              
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
                          )# box end
                          
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
                          )#box end
                   )#column end
                 )#fluidRow end
                 
                 
        ), #tab1 end
        
        ## 3.2 Tree plot and violin plot ----------------------------------------------------------
        tabItem( tabName = 'dashboard',
                 
                 fluidRow(
                   
                   box(title = "Updating Tree Plot", status="warning", solidHeader = TRUE, width=12,
                       h5("The results shown in this app were the analysis done within the cell groups shown here."),
                       #Tree plot
                       conditionalPanel(
                         condition = "input.go == 1",
                         img(src = "tree.PNG", width=550, height=480)  
                       ),
                       conditionalPanel(
                         condition = "input.update_marker == 1",
                         h3("Should be the updating tree plot")
                         ##automatic tree plot
                         ######################
                       )
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
                 ), #fluid row end
                 
                 fluidRow(
                   box(title = "What's next?", status="primary", solidHeader = TRUE, width=12,
                       h5("The violin plots here showed which level each cell type belongs to. "),
                       h5("For example, cell types which express marker 1 belong to Level 1. Cell types which don't express marker 1 would be used to draw violin plots for marker 2. 
                    Among those cell types, the ones which express marker 2 belong to Level 2."),
                       h5("To check the expression of ", strong("specific genes")," in cell types of different levels, please ", strong("expand Visualization (gene of interest) ") ,"panel on the sidebar menu."),
                       h5("Select the ",strong("dataset and genes of interest")," from the sidebar menu."),
                       actionBttn("go_to_r1", "Go to Visualization (gene of interest)", style = "jelly", color = "primary",size = "sm"),
                       
                       h5("To check the expression of ", strong("cell surface markers")," (defined by the app) in cell types of different levels, please ", strong("expand Visualization (marker selection) ") ,"panel on the sidebar menu."),
                       h5("Select the ",strong("dataset and cell type of interest")," from the sidebar menu."),
                       actionBttn("go_to_r2", "Go to Visualization (marker selection)", style = "jelly", color = "primary",size = "sm")
                   )#box end
                 )#fluid row end
        ), #tab end
        
        ## 3.3 Result (gene of interest)----------------------------------------------------------
        
        ### 3.3.1 Level 1 (gene of interest)----------------------------------------------------------
        
        tabItem( tabName = 'result1',
                 
                 
                 fluidRow(#dot plot
                   
                   h3("Level 1 Dot plot with genes of interest"),
                   h4(strong("Expand the siderbar menu on the left!")),
                   img(src = "tab1.PNG", width=240, height=33),
                   h4(code("Error messgage? Because you have not select the genes of interest!")),
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
                   h4(strong("Expand the siderbar menu on the left!")),
                   img(src = "tab1.PNG", width=240, height=33),
                   h4(code("Error messgage? Because you have not select the genes of interest!")),
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
                   h4(strong("Expand the siderbar menu on the left!")),
                   img(src = "tab1.PNG", width=240, height=33),
                   h4(code("Error messgage? Because you have not select the genes of interest!")),
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
                   h4(strong("Expand the siderbar menu on the left!")),
                   img(src = "tab1.PNG", width=240, height=33),
                   h4(code("Error messgage? Because you have not select the genes of interest!")),
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
                   h4(strong("Expand the siderbar menu on the left!")),
                   img(src = "tab1.PNG", width=240, height=33),
                   h4(code("Error messgage? Because you have not select the genes of interest!")),
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
                 fluidRow(
                   h4(strong("Expand the siderbar menu on the left, select the dataset and press 'Do marker selection'!")),
                   img(src = "tab2.PNG", width=240, height=33),
                   
                   box(title = "Marker selection results (Table)", status="primary", solidHeader = TRUE, width=12,
                       dataTableOutput("marker_gene1")
                   )
                 ),
                 
                 fluidRow(#dot plot
                   
                   h3("Level 1 Dot plot with specific cell surface markers"),
                   h4(code("Error messgage? Because you have not select the cell type of interest!")),
                   
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
                 fluidRow(
                   h4(strong("Expand the siderbar menu on the left, select the dataset and press 'Do marker selection'!")),
                   img(src = "tab2.PNG", width=240, height=33),
                   
                   box(title = "Marker selection results (Table)", status="primary", solidHeader = TRUE, width=12,
                       dataTableOutput("marker_gene2")
                   )
                 ),
                 fluidRow(#dot plot
                   h3("Level 2 Dot plot with specific cell surface markers"),
                   h4(code("Error messgage? Because you have not select the cell type of interest!")),
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
                 fluidRow(
                   h4(strong("Expand the siderbar menu on the left, select the dataset and press 'Do marker selection'!")),
                   img(src = "tab2.PNG", width=240, height=33),
                   
                   box(title = "Marker selection results (Table)", status="primary", solidHeader = TRUE, width=12,
                       dataTableOutput("marker_gene3")
                   )
                 ),
                 fluidRow(#dot plot
                   
                   h3("Level 3 Dot plot with specific cell surface markers"),
                   h4(code("Error messgage? Because you have not select the cell type of interest!")),
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
                 fluidRow(
                   h4(strong("Expand the siderbar menu on the left, select the dataset and press 'Do marker selection'!")),
                   img(src = "tab2.PNG", width=240, height=33),
                   
                   box(title = "Marker selection results (Table)", status="primary", solidHeader = TRUE, width=12,
                       dataTableOutput("marker_gene4")
                   )
                 ),
                 fluidRow(#dot plot
                   
                   h3("Level 4 Dot plot with specific cell surface markers"),
                   h4(code("Error messgage? Because you have not select the cell type of interest!")),
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
                 fluidRow(
                   h4(strong("Expand the siderbar menu on the left, select the dataset and press 'Do marker selection'!")),
                   img(src = "tab2.PNG", width=240, height=33),
                   
                   box(title = "Marker selection results (Table)", status="primary", solidHeader = TRUE, width=12,
                       dataTableOutput("marker_gene5")
                   )
                 ),
                 
                 fluidRow(#dot plot
                   
                   h3("Level 5 Dot plot with specific cell surface markers"),
                   h4(code("Error messgage? Because you have not select the cell type of interest!")),
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
                           h4("For the tutorial and FAQs, please visit",a("our GitHub", href="https://github.com/CompEpigen/digital_sorter", target="_blank" ),"." ),
                           
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
    }else {
      loginpage
    }
  })
  
  shinybusy::show_modal_spinner(text = "Loading datasets...") # show the modal window
  
  if(F){ 
    # read me
    # the dataset put in folder data/should be a SplitObject of seurat.
    # Create split object (split by datasets)
    # use separated R script 20211208_merge_seuratobject_addsplits_inspect.R
    dataset <- readRDS("data/NEW_datasets.rds")
    datasets <- SplitObject(dataset, split.by = "dataset_origin")
    saveRDS(datasets,"data/datasets.rds")
  }
  if(F){ 
    # read me
    # The markerlist.rds in folder data/ should be a list including master markers of interested cancer types.
    markerlist <- list("Lung"=c("PTPRC","EPCAM","PECAM1","MME"),
                       "t.b.c."=c("t.b.c."))
    saveRDS(markerlist,"./data/markerlist.rds")
  }
  ## read datasets ####
  
  markerlist <- markerlist
  
  genelist <- reactive({
    if(input$genelists == "Protein coding genes"){ 
      return(prot.coding.genes)
    }else if(input$genelists == "Cell surface markers"){return(Union.cell.surface.marker)} 
  })
  
  list_samples_disease <- reactive({
    list_samples_disease <- list()
    for(i in 1: length(datasets)){
      new_list <- list(as.character(unique(sort(datasets[[i]]@meta.data[["disease"]]))))
      list_samples_disease <- c(list_samples_disease,new_list)
    }
    names(list_samples_disease) <- names(datasets)
    return(list_samples_disease)})
  
  
  
  
  ## settings for plots ####
  cols <- reactive(RColorBrewer::brewer.pal(length(levels(datasets[[1]]@meta.data[["disease"]])), name="Paired"))
  
  aes_list <- reactive({theme(legend.text=element_text(size=12),
                              axis.text=element_text(size=12),
                              axis.title=element_text(size=14),
                              legend.title=element_text(size=12))
  })
  
  shinybusy::remove_modal_spinner() # remove it when done
  
  
  ## interactive markers #### 
  
  output$master_markers <- renderUI({
    myMarkers <- markerlist[[input$cancer]]
    selectizeInput(inputId = "marker", label= "Select master markers (default selected):",
                   choices = myMarkers,
                   selected = myMarkers,
                   multiple = TRUE
    )
  })
  
  ## interactive UI layout ####
  marker_num <- reactive({
    length(input$marker)
  })
  
  ### home violinplot tabPanel ####
  output$ui_tabBox <- renderUI({
    lapply(1:(marker_num()), function(i){
      tabPanel(paste0("Marker ",i), 
               withSpinner(plotOutput(paste0("plotv_h",i))))
    }) 
    
  })
  
  ### sidebar Results- menuSubItem ####
  # output$ui_menuSubItem <- renderUI({
  #    lapply(1:(marker_num()+1), function(i){
  #     menuSubItem(paste0("Level ",i), tabName = paste0("result",i))
  #    }) 
  #format would be different
  # })
  
  ## select cohort #### 
  output$cohort_h <- renderUI({
    if(input$cancer == "Lung"){ 
      myCohort <- as.character(names(datasets))
    } 
    else{myCohort <- c("t.b.c.")} 
    
    selectizeInput(inputId = "cohort", label= "Select dataset of interest:", choices = myCohort, 
                   selected = NULL)
  })
  
  output$cohort1 <- renderUI({
    if(input$cancer == "Lung"){ 
      myDataset <-as.character(names(datasets))
    } 
    else{myDataset <- c("t.b.c.")} 
    
    selectizeInput(inputId = "cohort1", label= "Select dataset of interest:", choices = myDataset, 
                   selected = NULL)
  })
  
  output$cohort2 <- renderUI({
    if(input$cancer == "Lung"){ 
      myDataset2 <-as.character(names(datasets))
    } 
    else{myDataset <- c("t.b.c.")} 
    
    selectizeInput(inputId = "cohort2", label= "Select dataset of interest:", choices = myDataset2, 
                   selected = NULL)
  })
  
  ## Event reactive objects after cohort selected #### 
  #for plots in home section
  dataset_reactive_h <- reactive({
    paste(input$cohort[1:length(input$cohort)])
  })
  
  output$dataset_h <-  renderText({input$cohort })
  output$cancer_type <-  renderText({input$cancer })
  
  #for plots in home section
  ob_reactive_h <- reactive({
    ob_reactive <- datasets[[dataset_reactive_h()]]
    ob_reactive@meta.data$annotation.l1=factor(ob_reactive@meta.data$annotation.l1,levels=unique(sort(datasets[[1]]@meta.data$annotation.l1)))
    ob_reactive@meta.data$disease=factor(ob_reactive@meta.data$disease,levels=unique(sort(ob_reactive@meta.data$disease)))
    
    return(ob_reactive)
  })
  
  
  #for plots in results section
  ob_reactive1 <- reactive({
    ob_reactive1 <- datasets[[input$cohort1]]
    ob_reactive1@meta.data$annotation.l1=factor(ob_reactive1@meta.data$annotation.l1,levels=unique(sort(datasets[[1]]@meta.data$annotation.l1)))
    return(ob_reactive1)
  })
  
  ob_reactive_split1<- reactive({
    ob_reactive_split1 <- SplitObject(ob_reactive1(), split.by = "disease")
    
    return(ob_reactive_split1)
  })
  ob_selected1 <- reactive({
    ob_selected1 <- ob_reactive_split1()[[input$sample1]]
    
    return(ob_selected1)
  })
  
  
  #for plots in results2 section
  ob_reactive2 <- reactive({
    ob_reactive2 <- datasets[[input$cohort2]]
    ob_reactive2@meta.data$annotation.l1=factor(ob_reactive2@meta.data$annotation.l1,levels=unique(sort(datasets[[1]]@meta.data$annotation.l1)))
    return(ob_reactive2)
  })
  ob_reactive_split2<- reactive({
    ob_reactive_split2 <- SplitObject(ob_reactive2(), split.by = "disease")
    
    return(ob_reactive_split2)
  })
  
  ob_selected2<- reactive({
    ob_selected2 <- ob_reactive_split2()[[input$sample2]]
    
    return(ob_selected2)
  })
  
  
  
  ## select sample ####  
  
  output$sample_se1 <- renderUI({
    mySample <- list_samples_disease()[[input$cohort1]]
    
    selectizeInput(inputId = "sample1","Select sample types:",
                   choices <- mySample, selected = NULL,
                   multiple = F)
  })
  
  output$sample_se2 <- renderUI({
    mySample2 <- list_samples_disease()[[input$cohort2]]
    
    selectizeInput(inputId = "sample2","Select sample types:",
                   choices <- mySample2, selected = NULL,
                   multiple = F)
  })
  
  ## Select Section-genes ####  
  output$gene <- renderUI({
    selectizeInput(inputId = "gene", label= "Select cell surface marker genes for the dot plots:",
                   choices = genelist(),
                   selected = NULL,
                   multiple = TRUE,
                   options = list(
                     placeholder = 'Type to search for gene',
                     onInitialize = I('function() { this.setValue(""); }')
                   )
    )
  })
  
  
  ## Select self defined master markers ####  
  output$own_markers <- renderUI({
    selectizeInput(inputId = "own_markers", label= "Select cell surface marker genes for stratification:",
                   choices = genelist(),
                   selected = NULL,
                   multiple = TRUE,
                   options = list(
                     placeholder = 'Type to search for gene',
                     onInitialize = I('function() { this.setValue(""); }')
                   )
    )
  })
  
  
  
  ## Select cell types for selecting markers ####
  output$celltype <- renderUI({
    if (!is.null(marker_gene_table())){  
      celltypes <- as.character(unique(sort(marker_gene_table()$annotation.l2)))
      selectizeInput(inputId = "celltype", label= "Select cell types for the features of dot plots:",
                     choices = celltypes,
                     selected = celltypes[1],
                     multiple = F,
                     options = list(
                       placeholder = 'Type to search for cell types whose markers were found',
                       onInitialize = I('function() { this.setValue(""); }')
                     )
      )
    }
  })
  
  
  ## Event reactive objects #### 
  
  marker_list <- reactive({ 
    paste(input$marker[1:length(input$marker)])
  })
  
  plot_marker <- eventReactive(input$cancer == "Lung", {
    plot_marker <- c()
    for(i in 1:length(markerlist[[input$cancer]])){
      markers <- c(paste0(markerlist[[input$cancer]][i],"+"),  paste0(markerlist[[input$cancer]][i],"-"))
      plot_marker <- c(plot_marker,markers)
    }
    
    return(plot_marker)}
  )
  
  
  
  
  genes <- reactive({ 
    paste(input$gene[1:length(input$gene)])
  })
  
  cell_chosen <- reactive({ 
    paste(input$celltype[1:length(input$celltype)])
  })
  
  #could pre-run this function for all the datasets/splits to save user's time
  ## marker selection function----------------
  marker_gene_table <- eventReactive(input$dplot2, { 
    shinybusy::show_modal_progress_line(text = "Get marker genes...") # show the modal window
    
    #level 1 
    ob1 <- SplitObject(ob_selected2(), split.by = "split1") 
    ob_split1 = ob1[[plot_marker()[1]]]
    seurat1 <- getMarkers(
      ob_split1,
      assay = 'RNA',
      organism = 'hg',
      groups = c('annotation.l2'),
      name = 'cerebro_seurat',
      only_pos = F,
      min_pct = 0.5,
      thresh_logFC = 1,
      thresh_p_val = 0.05,
      test = 'wilcox', 
      verbose = TRUE
    )
    marker_gene_table <- seurat1@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]] #3708
    
    subset_marker_gene1 <- marker_gene_table %>% 
      filter(marker_gene_table$gene %in% genelist())%>%
      arrange(annotation.l2,desc(abs(avg_log2FC)),p_val)
    update_modal_progress(
      0.25,
      text = "Level 1 completed",
      session = shiny::getDefaultReactiveDomain()
    )
    #level2
    ob2 = ob1[[plot_marker()[2]]]
    ob2 <- SplitObject(ob2, split.by = "split2") 
    ob_split2 = ob2[[plot_marker()[3]]]
    
    seurat2 <- getMarkers(
      ob_split2,
      assay = 'RNA',
      organism = 'hg',
      groups = c('annotation.l2'),
      name = 'cerebro_seurat',
      only_pos = F,
      min_pct = 0.5,
      thresh_logFC = 1,
      thresh_p_val = 0.05,
      test = 'wilcox', 
      verbose = TRUE
    )
    marker_gene_table <- seurat2@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]] #3708
    
    subset_marker_gene2 <- marker_gene_table %>% 
      filter(marker_gene_table$gene %in% genelist())%>%
      arrange(annotation.l2,desc(abs(avg_log2FC)),p_val)
    update_modal_progress(
      0.50,
      text = "Level 2 completed",
      session = shiny::getDefaultReactiveDomain()
    )
    #level3
    ob3 = ob2[[plot_marker()[4]]]
    ob3 <- SplitObject(ob3, split.by = "split3") 
    ob_split3 = ob3[[plot_marker()[5]]]
    
    seurat3 <- getMarkers(
      ob_split3,
      assay = 'RNA',
      organism = 'hg',
      groups = c('annotation.l2'),
      name = 'cerebro_seurat',
      only_pos = F,
      min_pct = 0.5,
      thresh_logFC = 1,
      thresh_p_val = 0.05,
      test = 'wilcox', 
      verbose = TRUE
    )
    marker_gene_table <- seurat3@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]] #3708
    
    subset_marker_gene3 <- marker_gene_table %>% 
      filter(marker_gene_table$gene %in% genelist())%>%
      arrange(annotation.l2,desc(abs(avg_log2FC)),p_val)
    update_modal_progress(
      0.75,
      text = "Level 3 completed",
      session = shiny::getDefaultReactiveDomain()
    )
    #level4
    ob_split4 = ob3[[plot_marker()[6]]]
    
    seurat4 <- getMarkers(
      ob_split4,
      assay = 'RNA',
      organism = 'hg',
      groups = c('annotation.l2'),
      name = 'cerebro_seurat',
      only_pos = F,
      min_pct = 0.5,
      thresh_logFC = 1,
      thresh_p_val = 0.05,
      test = 'wilcox', 
      verbose = TRUE
    )
    marker_gene_table <- seurat4@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]] #3708
    
    subset_marker_gene4 <- marker_gene_table %>% 
      filter(marker_gene_table$gene %in% genelist())%>%
      arrange(annotation.l2,desc(abs(avg_log2FC)),p_val)
    update_modal_progress(
      1,
      text = "Level 4 completed",
      session = shiny::getDefaultReactiveDomain()
    )
    #combine
    subset_marker_gene <- rbind(subset_marker_gene1,subset_marker_gene2,subset_marker_gene3,subset_marker_gene4)
    
    shinybusy::remove_modal_progress()
    return(subset_marker_gene)
  })
  
  lapply(1:5, function(i) {
    outputId <- paste0("dotplot_title", i)
    output[[outputId]] <- renderText({
      title <- list(1,
                    c(2,3),
                    c(2,4,5),
                    c(2,4,6,7),
                    c(2,4,6,8))
      if(input$cancer == "Lung"){plot_marker()[title[[i]]]}else {"to be continued"}
    })
  })
  
  lapply(1:5, function(i) {
    outputId <- paste0("dotplot2_title", i)
    output[[outputId]] <- renderText({
      title <- list(1,
                    c(2,3),
                    c(2,4,5),
                    c(2,4,6,7),
                    c(2,4,6,8))
      if(input$cancer == "Lung"){plot_marker()[title[[i]]]}else {"to be continued"}
    })
  })
  
  output$download_marker_csv <- downloadHandler(
    filename = "Marker_selection_table.csv",
    content = function(file) {
      write.csv(marker_gene_table(), file)
    })
  ## output marker gene table ####
  marker_gene_table_output <- reactive({
    marker_gene_table <- marker_gene_table()
    marker_gene_table_output <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
    return(marker_gene_table_output)
  })
  
  lapply(1:5, function(i) {
    outputId <- paste0("marker_gene", i)
    output[[outputId]] <- renderDataTable(marker_gene_table_output(), options = list(pageLength = 10))
  })
  
  
  ##Change the selected tab on the client ####
  observeEvent(input$go, {
    updateTabItems( session = getDefaultReactiveDomain(), "sidebar",
                    selected = "dashboard"
    )
  })
  
  observeEvent(input$plotv, {
    updateTabsetPanel(session, "vplots",
                      selected = "Marker 1"
    )
  })
  
  observeEvent(input$dplot2, {
    updateTabItems( session = getDefaultReactiveDomain(), "sidebar",
                    selected = "2result1"
    )
  })
  
  observeEvent(input$update_marker, {
    updateTabItems( session = getDefaultReactiveDomain(), "sidebar",
                    selected = "dashboard"
    )
  })
  
  observeEvent(input$go_to_r1, {
    updateTabItems( session = getDefaultReactiveDomain(), "sidebar",
                    selected = "result1"
    )
  })
  observeEvent(input$go_to_r1, {
    shinyjs::alert("Expand the sidebar menu 'Visualization (gene of interest)' on the left!")
  })
  
  
  observeEvent(input$go_to_r2, {
    updateTabItems( session = getDefaultReactiveDomain(), "sidebar",
                    selected = "2result1"
    )
  })
  observeEvent(input$go_to_r2, {
    shinyjs::alert("Expand the sidebar menu 'Visualization (marker selection)' on the left!")
  })
  
  
  
  
  ## Violin Plot #### 
  ### home ####
  
  v_h1 <-eventReactive(input$plotv,{
    #ob_selected <- subset(x = ob_reactive_h(), subset = dataset_origin == input$cohort)
    
    VlnPlot(ob_reactive_h(),marker_list()[1], split.by = "disease", group.by = "annotation.l2", cols=cols(),
            sort = TRUE,pt.size = 0, combine = FALSE)
  })
  v_h2 <-eventReactive(input$plotv,{
    #ob_selected <- subset(x = ob_reactive_h(), subset = dataset_origin == input$cohort)
    
    VlnPlot(ob_reactive_h(),marker_list()[2], split.by = "disease", group.by = "annotation.l2", cols=cols(),
            sort = TRUE,pt.size = 0, combine = FALSE)
  })
  v_h3 <-eventReactive(input$plotv,{
    #ob_selected <- subset(x = ob_reactive_h(), subset = dataset_origin == input$cohort)
    
    VlnPlot(ob_reactive_h(),marker_list()[3], split.by = "disease", group.by = "annotation.l2", cols=cols(),
            sort = TRUE,pt.size = 0, combine = FALSE)
  })
  v_h4 <-eventReactive(input$plotv,{
    #ob_selected <- subset(x = ob_reactive_h(), subset = dataset_origin == input$cohort)
    
    VlnPlot(ob_reactive_h(),marker_list()[4], split.by = "disease", group.by = "annotation.l2", cols=cols(),
            sort = TRUE,pt.size = 0, combine = FALSE)
  })
  
  hvlist <- reactive(list(v_h1(),v_h2(),v_h3(),v_h4()))
  lapply(1:4, function(i) {
    outputId <- paste0("plotv_h", i)
    output[[outputId]] <- renderPlot(hvlist()[[i]])  
  })
  
  
  
  ## Dot plot ####
  
  d1 <- reactive({
    ob1_1 <- subset(x = ob_selected1(), subset = split1 == plot_marker()[1])
    
    DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
      RotatedAxis()+ aes_list()
    
  })
  
  d2 <- reactive({
    ob1_1 <- subset(x = ob_selected1(), subset = split2 == plot_marker()[3])
    
    DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
      RotatedAxis()+ aes_list()
    
  })
  
  d3 <- reactive({
    ob1_1 <- subset(x = ob_selected1(), subset = split3 == plot_marker()[5])
    
    DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
      RotatedAxis()+ aes_list()
    
  })
  d4 <- reactive({
    ob1_1 <- subset(x = ob_selected1(),subset = split4 == plot_marker()[7])
    
    
    DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
      RotatedAxis()+ aes_list()
    
  })
  
  d5 <- reactive({
    ob1_1 <- subset(x = ob_selected1(), subset = split4 == plot_marker()[8])
    
    
    DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
      RotatedAxis()+aes_list()
    
  })
  
  dlist <- reactive(list(d1(),d2(),d3(),d4(),d5()))
  
  lapply(1:5, function(i) {
    outputId <- paste0("plotd", i)
    output[[outputId]] <- renderPlot(dlist()[[i]])
  })
  
  
  
  ## table next to the dot plot ####
  
  d1_table <- reactive({
    table_d <- dlist()[[1]][["data"]][, c(3,4,1,2,5)]
    table_d <- table_d%>% arrange(desc(table_d$avg.exp))
    return(table_d)
  })
  d2_table <- reactive({
    table_d <- dlist()[[2]][["data"]][, c(3,4,1,2,5)]
    table_d <- table_d%>% arrange(desc(table_d$avg.exp))
    return(table_d)
  })
  d3_table <- reactive({
    table_d <- dlist()[[3]][["data"]][, c(3,4,1,2,5)]
    table_d <- table_d%>% arrange(desc(table_d$avg.exp))
    return(table_d)
  })
  d4_table <- reactive({
    table_d <- dlist()[[4]][["data"]][, c(3,4,1,2,5)]
    table_d <- table_d%>% arrange(desc(table_d$avg.exp))
    return(table_d)
  })
  d5_table <- reactive({
    table_d <- dlist()[[5]][["data"]][, c(3,4,1,2,5)]
    table_d <- table_d%>% arrange(desc(table_d$avg.exp))
    return(table_d)
  })
  
  dlist_table <- reactive(list(d1_table(),d2_table(),d3_table(),d4_table(),d5_table()))
  lapply(1:5, function(i) {
    outputId <- paste0("table", i)
    output[[outputId]] <- renderDataTable(dlist_table()[[i]] , 
                                          options = list(pageLength = 5))
  })
  
  ## Dot plot2 for marker selection ####
  
  d1_cellMark <- reactive({
    ob1_1 <- subset(x = ob_selected2(), subset = split1 == plot_marker()[1])
    
    #cell_types1 <- unique(sort(ob1_1@meta.data[["annotation.l2"]]))
    #marker_gene_table <- marker_gene_table()
    #marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_types1,]
    #gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
    
    marker_gene_table <- marker_gene_table()
    marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
    gene_dge <- unique(marker_gene_table_sub$gene)[1:input$n_genes]
    
    DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
      RotatedAxis()+ aes_list() 
    
  })
  
  d2_cellMark <-reactive({
    ob1_1 <- subset(x = ob_selected2(), subset = split2 == plot_marker()[3])
    
    marker_gene_table <- marker_gene_table()
    marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
    gene_dge <- unique(marker_gene_table_sub$gene)[1:input$n_genes]
    
    DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
      RotatedAxis()+ aes_list()
    
  })
  
  d3_cellMark <- reactive({
    ob1_1 <- subset(x = ob_selected2(), subset = split3 == plot_marker()[5])
    
    marker_gene_table <- marker_gene_table()
    marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
    gene_dge <- unique(marker_gene_table_sub$gene)[1:input$n_genes]
    
    DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
      RotatedAxis()+ aes_list()
    
  })
  
  d4_cellMark <- reactive({
    ob1_1 <- subset(x = ob_selected2(),subset = split4 == plot_marker()[7])
    
    
    marker_gene_table <- marker_gene_table()
    marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
    gene_dge <- unique(marker_gene_table_sub$gene)[1:input$n_genes]
    
    DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
      RotatedAxis()+aes_list()
    
  })
  
  d5_cellMark <- reactive({
    ob1_1 <- subset(x = ob_selected2(), subset = split4 == plot_marker()[8])
    
    
    marker_gene_table <- marker_gene_table()
    marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
    gene_dge <- unique(marker_gene_table_sub$gene)[1:input$n_genes]
    
    DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
      RotatedAxis()+ aes_list()
    
  })
  
  dlist_cellMark <- reactive(list(d1_cellMark(),d2_cellMark(),d3_cellMark(),d4_cellMark(),d5_cellMark()))
  
  lapply(1:5, function(i) {
    outputId <- paste0("plotd_cellMark", i)
    output[[outputId]] <- renderPlot(dlist_cellMark()[[i]])
  })
  
  
  
  ## table next to the dot plot ####
  d1_cellMark_table <- reactive({
    table <- dlist_cellMark()[[1]][["data"]][, c(3,4,1,2,5)]
    table <- table %>% arrange(desc(table$avg.exp))
    return(table)
  })
  d2_cellMark_table <- reactive({
    table <- dlist_cellMark()[[2]][["data"]][, c(3,4,1,2,5)]
    table <- table %>% arrange(desc(table$avg.exp))
    return(table)
  })
  d3_cellMark_table <- reactive({
    table <- dlist_cellMark()[[3]][["data"]][, c(3,4,1,2,5)]
    table <- table %>% arrange(desc(table$avg.exp))
    return(table)
  })
  d4_cellMark_table <- reactive({
    table <- dlist_cellMark()[[4]][["data"]][, c(3,4,1,2,5)]
    table <- table %>% arrange(desc(table$avg.exp))
    return(table)
  })
  d5_cellMark_table <- reactive({
    table <- dlist_cellMark()[[5]][["data"]][, c(3,4,1,2,5)]
    table <- table %>% arrange(desc(table$avg.exp))
    return(table)
  })
  
  dlist_cellMark_table <- reactive(list(d1_cellMark_table(),d2_cellMark_table(),
                                        d3_cellMark_table(),d4_cellMark_table(),
                                        d5_cellMark_table()))
  
  lapply(1:5, function(i) {
    outputId <- paste0("table_cellMark", i)
    output[[outputId]] <- renderDataTable(dlist_cellMark_table()[[i]], 
                                          options = list(pageLength = 5))
  })
}


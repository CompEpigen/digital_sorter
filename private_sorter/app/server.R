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
      ## 1st tab show the preprocessing dashboard  Main dashboard -----------
      menuItem( "Define My Own Master Markers", tabName = 'processing', icon = icon('tachometer-alt')),
      ## 2st tab show the Main dashboard-----------
      menuItem( "Tree Plot & Violin Plot", tabName = 'dashboard', icon = icon('tree')),
      ## 3nd tab shows results ----------
      menuItem( "Visualization (gene of interest)", tabName = "results", icon = icon('list'), startExpanded = F,
                
                lapply(1:5, function(i){
                  menuSubItem(paste0("Level ",i), tabName = paste0("result",i))
                }),
                #uiOutput("ui_menuSubItem"), #format would be different
                
                h4("Input section for each level"),
                #select dataset (choices depends on cancer)                 
                uiOutput("cohort1"),
                
                uiOutput("sample_se1"),
                
                uiOutput("gene"),
                
                
                
                actionBttn("dplot", "Plot/ Renew Dot Plots!", 
                           style = "jelly", color = "warning", size = "sm"), 
                h5("Go to each level and check out the results!")
                
                
      ),
      
      menuItem( "Visualization (marker selection)", tabName = "results2", icon = icon('brain'), startExpanded = F,
                
                lapply(1:5, function(i){
                  menuSubItem(paste0("Level ",i), tabName = paste0("2result",i))
                }),
                #uiOutput("ui_menuSubItem"), #format would be different
                
                h4("Input section for each level"),
                #select dataset (choices depends on cancer)                 
                uiOutput("cohort2"),
                
                uiOutput("sample_se2"),
                
                fluidRow(
                  column(width=4,offset=1,
                         downloadButton('download_marker_csv', 'Download Marker Selection Results')
                  ),
                ),
                uiOutput("celltype"),
                
                actionBttn("dplot2", "Plot/ Renew Dot Plots!", 
                           style = "jelly", color = "warning", size = "sm"), 
                h5("Go to each level and check out the results!")
                
                
      ),
      
      
      
      ## 3rd tab Data source, definition , i.e., help ---------------
      menuItem( "FAQs", tabName = 'help', icon = icon('question-circle') )
    )
  })
    
  output$body <- renderUI({
    if (USER$login == TRUE ) {  
      
      ## 3.1 Dashboard body --------------
      tabItems(
        ## 3.1.0 Pre-processing dashboard ----------------------------------------------------------
        tabItem( tabName = 'processing',
                 fluidRow(
                   column(width=4,
                          box(title = "Datasets stratified with default master markers", status= "success", solidHeader = TRUE, width = NULL,
                              #cancer
                              selectizeInput(inputId ="cancer", label= "Select cancer types:", choices = c("Lung","t.b.c."), selected = "Lung"),
                              #show master markers after selecting cancer type
                              uiOutput("master_markers"),
                              #select dataset (choices depends on cancer)                 
                              uiOutput("cohort_h"),
                              
                              actionBttn("go", "Confirm & Go!", style = "jelly", color = "success",size = "sm")
                          ),# box end
                          box(title = "Define your own master markers", status= "warning", solidHeader = TRUE, width = NULL,
                              uiOutput("own_markers"),
                              textInput(inputId="own_markers_name", label="Name your own master markers (e.g. Lung1)", value = "", width = NULL, placeholder = NULL),
                              actionBttn("update_marker", "Update marker list and tree plot!", style = "jelly", color = "warning",size = "sm"),
                              h4("This is the panel for user to define their own master markers."),
                              h4("If the cancer type is not",span(" Lung", style = "color:blue") ,", user should upload the reference expression data to define the levels.")
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
        
        ## 3.1.1 Main dashboard ----------------------------------------------------------
        tabItem( tabName = 'dashboard',
                 
                 fluidRow(
                   
                   box(title = "Updating Tree Plot", status="warning", solidHeader = TRUE, width=8,
                       h5("The results shown in this app were the analysis done within the cell groups shown here."),
                       #Tree plot
                       h3("Should be the updating tree plot")
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
        
        ## 3.2.1 Result1 ----------------------------------------------------------
        tabItem( tabName = '2result1',
                 
                 fluidRow(#dot plot
                   tags$hr(),
                   h3("Level 1 Dot plot"),
                   h4(textOutput("dotplot2_title1")),
                   column(width=6, style = "height:200px;",
                          withSpinner(plotOutput("plotd_cellMark1") )
                          
                          
                   ),
                   column(width=4,
                          withSpinner(dataTableOutput("table_cellMark1") )
                          
                   )
                   
                 )#fluid row end
        ), #tab2 end
        
        ## 3.2.2 Result2 ----------------------------------------------------------
        
        tabItem( tabName = '2result2',
                 
                 fluidRow(#dot plot
                   tags$hr(),
                   h3("Level 2 Dot plot"),
                   h4(textOutput("dotplot2_title2")),
                   column(width=6, style = "height:200px;",
                          withSpinner(plotOutput("plotd_cellMark2") )
                          
                   ),
                   column(width=4,
                          withSpinner(dataTableOutput("table_cellMark2") )
                          
                   )
                   
                 )#fluid row end
        ), #tab3 end
        
        ## 3.2.3 Result3 ----------------------------------------------------------
        
        tabItem( tabName = '2result3',
                 
                 fluidRow(#dot plot
                   tags$hr(),
                   h3("Level 3 Dot plot"),
                   h4(textOutput("dotplot2_title3")),
                   
                   column(width=6, style = "height:200px;",
                          
                          withSpinner( plotOutput("plotd_cellMark3") )
                          
                   ),
                   column(width=4,
                          withSpinner(dataTableOutput("table_cellMark3") )
                          
                   )
                   
                 )#fluid row end
        ), #tab4 end
        
        ## 3.2.4 Result4 ----------------------------------------------------------
        
        tabItem( tabName = '2result4',
                 
                 fluidRow(#dot plot
                   tags$hr(),
                   h3("Level 4 Dot plot"),
                   h4(textOutput("dotplot2_title4")),
                   
                   column(width=6, style = "height:200px;",
                          withSpinner(plotOutput("plotd_cellMark4")  )
                          
                   ),
                   column(width=4,
                          withSpinner(dataTableOutput("table_cellMark4") )
                          
                          
                   )      
                 )#fluid row end
        ), #tab5 end
        ## 3.2.5 Result5 ----------------------------------------------------------
        tabItem( tabName = '2result5',
                 
                 
                 fluidRow(#dot plot
                   tags$hr(),
                   h3("Level 5 Dot plot"),
                   h4(textOutput("dotplot2_title5")),
                   
                   column(width=6, style = "height:200px;",
                          withSpinner(plotOutput("plotd_cellMark5") )
                          
                   ),
                   column(width=4,
                          withSpinner(dataTableOutput("table_cellMark5") )
                          
                          
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
    }else {
      loginpage
    }
  })
  
  shinybusy::show_modal_spinner(text = "Loading datasets...") # show the modal window
  
  if(F){ 
    # read me
    # the dataset put infolder data/should be a SplitObject of seurat.
    # Create split object (split by datasets)
    # use separated R script 20211208_merge_seuratobject_addsplits_inspect.R
    dataset <- readRDS("private_sorter/data/NEW_datasets.rds")
    split <- as.character(unique(sort(dataset@meta.data$dataset_origin)))
    datasets <- SplitObject(dataset, split.by = "dataset_origin")
    saveRDS(datasets,"private_sorter/data/datasets.rds")
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
  genelist <- reactive({return(Union.cell.surface.marker)})
  
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
  
  #for plots in home section
  ob_reactive_h <- reactive({
    ob_reactive <- datasets[[dataset_reactive_h()]]
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
    celltypes <- as.character(unique(sort(datasets[[1]]@meta.data[["annotation.l2"]])))
    selectizeInput(inputId = "celltype", label= "Select cell types for automatically choose the features of dot plots:",
                   choices = celltypes,
                   selected = NULL,
                   multiple = F,
                   options = list(
                     placeholder = 'Type to search for cell types',
                     onInitialize = I('function() { this.setValue(""); }')
                   )
    )
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
    
    return(plot_marker)})
  
  
  
  
  genes <- reactive({ 
    paste(input$gene[1:length(input$gene)])
  })
  
  cell_chosen <- eventReactive(input$dplot2, { 
    paste(input$celltype[1:length(input$celltype)])
  })
  
  #could pre-run this function for all the datasets/splits to save user's time
  marker_gene_table <- eventReactive(input$cellMark,{
    shinybusy::show_modal_spinner(text = "Get marker genes...") # show the modal window
    seurat <- getMarkers(
      ob_reactive2(),
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
    marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]] #3708
    
    subset_marker_gene <- marker_gene_table %>% 
      filter(marker_gene_table$gene %in% genelist())
    subset_marker_gene <- subset_marker_gene %>%
      filter(subset_marker_gene$p_val_adj<1)%>% 
      arrange(annotation.l2,p_val)
    
    shinybusy::remove_modal_spinner()
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
  
  ##Change the selected tab on the client ####
  observeEvent(input$dplot, {
    updateTabItems( session = getDefaultReactiveDomain(), "sidebar",
                    selected = "result2"
    )
  })
  
  observeEvent(input$dplot2, {
    updateTabItems( session = getDefaultReactiveDomain(), "sidebar",
                    selected = "2result2"
    )
  })
  
  observeEvent(input$go, {
    updateTabItems( session = getDefaultReactiveDomain(), "sidebar",
                    selected = "dashboard"
    )
  })
  
  observeEvent(input$update_marker, {
    updateTabItems( session = getDefaultReactiveDomain(), "sidebar",
                    selected = "dashboard"
    )
  })
  
  ## Violin Plot #### 
  ### home ####
  
  v_h1 <-eventReactive(input$go,{
    VlnPlot(ob_reactive_h(),marker_list()[1], split.by = "disease", group.by = "annotation.l2", cols=cols(),
            sort = TRUE,pt.size = 0, combine = FALSE)
  })
  v_h2 <-eventReactive(input$go,{
    VlnPlot(ob_reactive_h(),marker_list()[2], split.by = "disease", group.by = "annotation.l2", cols=cols(),
            sort = TRUE,pt.size = 0, combine = FALSE)
  })
  v_h3 <-eventReactive(input$go,{
    VlnPlot(ob_reactive_h(),marker_list()[3], split.by = "disease", group.by = "annotation.l2", cols=cols(),
            sort = TRUE,pt.size = 0, combine = FALSE)
  })
  v_h4 <-eventReactive(input$go,{
    VlnPlot(ob_reactive_h(),marker_list()[4], split.by = "disease", group.by = "annotation.l2", cols=cols(),
            sort = TRUE,pt.size = 0, combine = FALSE)
  })
  
  hvlist <- reactive(list(v_h1(),v_h2(),v_h3(),v_h4()))
  lapply(1:4, function(i) {
    outputId <- paste0("plotv_h", i)
    output[[outputId]] <- renderPlot(hvlist()[[i]])  
  })
  
  #Level1: CD45 (PTPRC) 
  output$plotv1 <- renderPlotly({
    ob_selected <- ob_reactive_split1()[[input$sample1]]
    VlnPlot(ob_selected, marker_list()[1], split.by = "disease", group.by = "annotation.l2", cols=cols(),
            sort = TRUE,pt.size = 0, combine = FALSE)
    ggplotly(ggplot2::last_plot())
  }) 
  #Level2: EPCAM
  output$plotv2 <- renderPlotly({
    ob_selected <- ob_reactive_split1()[[input$sample1]]
    ob_selected_1 <- subset(x = ob_selected, subset = split1 == plot_marker()[2])
    
    VlnPlot(ob_selected_1, marker_list()[2], split.by = "disease", group.by = "annotation.l2", cols=cols(),
            sort = TRUE, pt.size = 0, combine = FALSE)
    ggplotly(ggplot2::last_plot())
  }) 
  
  
  #Level3: PECAM1 (CD31)
  output$plotv3 <- renderPlotly({
    ob_selected <- ob_reactive_split1()[[input$sample1]]
    ob_selected_1 <- subset(x = ob_selected, subset = split2 == plot_marker()[4])
    
    VlnPlot(ob_selected_1, marker_list()[3], split.by = "disease", group.by = "annotation.l2", cols=cols(),
            sort = TRUE,pt.size = 0, combine = FALSE)
    ggplotly(ggplot2::last_plot()) 
  }) 
  
  #Level4: MME (CD10)
  output$plotv4 <- renderPlotly({
    ob_selected <- ob_reactive_split1()[[input$sample1]]
    ob_selected_1 <- subset(x = ob_selected, subset = split3 == plot_marker()[6])
    
    VlnPlot(ob_selected_1, marker_list()[4], split.by = "disease", group.by = "annotation.l2", cols=cols(),
            sort = TRUE,pt.size = 0, combine = FALSE)
    ggplotly(ggplot2::last_plot()) 
    
  }) 
  
  ## Dot plot ####
  
  d1 <- eventReactive(input$dplot,ignoreInit = T,{
    ob_selected <- ob_reactive_split1()[[input$sample1]]
    
    ob1_1 <- subset(x = ob_selected, subset = split1 == plot_marker()[1])
    
    DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
      RotatedAxis()+ aes_list()
    
  })
  
  d2 <- eventReactive(input$dplot,{
    ob_selected <- ob_reactive_split1()[[input$sample1]]
    
    ob1_1 <- subset(x = ob_selected, subset = split2 == plot_marker()[3])
    
    DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
      RotatedAxis()+ aes_list()
    
  })
  
  d3 <- eventReactive(input$dplot,{
    ob_selected <- ob_reactive_split1()[[input$sample1]]
    
    ob1_1 <- subset(x = ob_selected, subset = split3 == plot_marker()[5])
    
    DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
      RotatedAxis()+ aes_list()
    
  })
  d4 <- eventReactive(input$dplot,{
    ob_selected <- ob_reactive_split1()[[input$sample1]]
    
    ob1_1 <- subset(x = ob_selected,subset = split4 == plot_marker()[7])
    
    
    DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
      RotatedAxis()+ aes_list()
    
  })
  
  d5 <- eventReactive(input$dplot,{
    ob_selected <- ob_reactive_split1()[[input$sample1]]
    
    ob1_1 <- subset(x = ob_selected, subset = split4 == plot_marker()[8])
    
    
    DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
      RotatedAxis()+aes_list()
    
  })
  
  dlist <- reactive(list(d1(),d2(),d3(),d4(),d5()))
  
  lapply(1:5, function(i) {
    outputId <- paste0("plotd", i)
    output[[outputId]] <- renderPlot(dlist()[[i]])
  })
  
  
  
  ## table next to the dot plot ####
  lapply(1:5, function(i) {
    outputId <- paste0("table", i)
    output[[outputId]] <- renderDataTable(dlist()[[i]][["data"]][, c(3,4,1,2,5)], 
                                          options = list(pageLength = 5))
  })
  
  ## Dot plot2 for marker selection ####
  
  d1_cellMark <- eventReactive(input$dplot2,ignoreInit = T,{
    ob_selected <- ob_reactive_split2()[[input$sample2]]
    
    ob1_1 <- subset(x = ob_selected, subset = split1 == plot_marker()[1])
    
    #cell_types1 <- unique(sort(ob1_1@meta.data[["annotation.l2"]]))
    #marker_gene_table <- marker_gene_table()
    #marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_types1,]
    #gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
    
    marker_gene_table <- marker_gene_table()
    marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
    gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
    
    DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
      RotatedAxis()+ aes_list() 
    
  })
  
  d2_cellMark <- eventReactive(input$dplot2,{
    ob_selected <- ob_reactive_split2()[[input$sample2]]
    
    ob1_1 <- subset(x = ob_selected, subset = split2 == plot_marker()[3])
    
    marker_gene_table <- marker_gene_table()
    marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
    gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
    
    DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
      RotatedAxis()+ aes_list()
    
  })
  
  d3_cellMark <- eventReactive(input$dplot2,{
    ob_selected <- ob_reactive_split2()[[input$sample2]]
    
    ob1_1 <- subset(x = ob_selected, subset = split3 == plot_marker()[5])
    
    marker_gene_table <- marker_gene_table()
    marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
    gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
    
    DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
      RotatedAxis()+ aes_list()
    
  })
  
  d4_cellMark <- eventReactive(input$dplot2,{
    ob_selected <- ob_reactive_split2()[[input$sample2]]
    
    ob1_1 <- subset(x = ob_selected,subset = split4 == plot_marker()[7])
    
    
    marker_gene_table <- marker_gene_table()
    marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
    gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
    
    DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
      RotatedAxis()+aes_list()
    
  })
  
  d5_cellMark <- eventReactive(input$dplot2,{
    ob_selected <- ob_reactive_split2()[[input$sample2]]
    
    ob1_1 <- subset(x = ob_selected, subset = split4 == plot_marker()[8])
    
    
    marker_gene_table <- marker_gene_table()
    marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
    gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
    
    DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
      RotatedAxis()+ aes_list()
    
  })
  
  dlist_cellMark <- reactive(list(d1_cellMark(),d2_cellMark(),d3_cellMark(),d4_cellMark(),d5_cellMark()))
  
  lapply(1:5, function(i) {
    outputId <- paste0("plotd_cellMark", i)
    output[[outputId]] <- renderPlot(dlist_cellMark()[[i]])
  })
  
  
  
  ## table next to the dot plot ####
  lapply(1:5, function(i) {
    outputId <- paste0("table_cellMark", i)
    output[[outputId]] <- renderDataTable(dlist_cellMark()[[i]][["data"]][, c(3,4,1,2,5)], 
                                          options = list(pageLength = 5))
  })
}


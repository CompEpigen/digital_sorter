library(shiny)
library(shinybusy)
library(Seurat)
#library(SeuratDisk) #for saving file
library(dplyr)
library(reshape2)
library(ggplot2)
library(plotly)
library(gridExtra)
library(shinyjs)
library(RColorBrewer)
ReadRDSFiles("data")

#server.R
server <- function(input, output, update_gene_list=F) {
 
  shinybusy::show_modal_spinner(text = "Loading datasets...") # show the modal window
  
  if(F){ 
    # read me
    # the dataset put infolder data/should be a SplitObject of seurat.
    # Create split object (split by datasets)
    # use separated R script 20211208_merge_seuratobject_addsplits_inspect.R
    dataset <- readRDS("data/NEW_datasets.rds")
    split <- as.character(unique(sort(dataset@meta.data$dataset_origin)))
    datasets <- SplitObject(dataset, split.by = "dataset_origin")
    saveRDS(datasets,"data/datasets.rds")
  }
  if(F){ 
    # read me
    # The markerlist.rds in folder data/ should be a list including master markers of interested cancer types.
    markerlist <- list("Lung"=c("PTPRC","EPCAM","PECAM1","MME"),
                       "t.b.c."=c("t.b.c."))
    saveRDS(markerlist,"digital_sorter/data/markerlist.rds")
  }
  ## read datasets ####
  
  markerlist <- markerlist
  genelist <- Union.cell.surface.marker
  
  list_samples_disease <- reactive({
    list_samples_disease <- list()
    for(i in 1: length(datasets)){
      new_list <- list(as.character(unique(sort(datasets[[i]]@meta.data[["disease"]]))))
      list_samples_disease <- c(list_samples_disease,new_list)
    }
    names(list_samples_disease) <- names(datasets)
    return(list_samples_disease)})

  

  
  ## update cell surface marker genes or not
  if(update_gene_list){
    cell.surface.marker <- get_cell_markers()
    new_cell_marker_gene <- SurfaceGenie.cell.surface.marker
    genelist <- union(new_cell_marker_gene,cell.surface.marker) 
  }
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
  output$ui_menuSubItem <- renderUI({
    lapply(1:(marker_num()+1), function(i){
      menuSubItem(paste0("Level ",i), tabName = paste0("result",i))
    }) 
    #format would be different
  })
   
  ## select cohort #### 
  output$cohort <- renderUI({
    if(input$cancer == "Lung"){ 
      myCohort <- as.character(names(datasets))
      } 
    else{myCohort <- c("t.b.c.")} 
   
    selectizeInput(inputId = "cohort", label= "Select dataset of interest:", choices = myCohort, 
                   selected = NULL)
  })

  output$cohort2 <- renderUI({
    if(input$cancer == "Lung"){ 
      myDataset <-as.character(names(datasets))
      } 
    else{myDataset <- c("t.b.c.")} 
    
    selectizeInput(inputId = "cohort2", label= "Select dataset of interest:", choices = myDataset, 
                   selected = NULL)
  })

  ## Event reactive objects after cohort selected #### 
  #for plots in home section
  dataset_reactive <- reactive({
    paste(input$cohort[1:length(input$cohort)])
  })
  #for plots in results section
  dataset_reactive2 <- reactive({
    paste(input$cohort2[1:length(input$cohort2)])
  })
  #for plots in home section
  ob_reactive <- reactive({
    ob_reactive <- datasets[[dataset_reactive()]]
    return(ob_reactive)
  })
  #for plots in results section
  ob_reactive2 <- reactive({
    ob_reactive2 <- datasets[[dataset_reactive2()]]
    ob_reactive2@meta.data$annotation.l1=factor(ob_reactive2@meta.data$annotation.l1,levels=unique(sort(datasets[[1]]@meta.data$annotation.l1)))
    return(ob_reactive2)
  })

  #for plots in results section
  ob_reactive_split<- reactive({
    ob_reactive_split <- SplitObject(ob_reactive2(), split.by = "disease")
    
    return(ob_reactive_split)
  })
  ## select sample ####  
  output$sample_se <- renderUI({
    mySample <- list_samples_disease()[[dataset_reactive2()]]
    
    selectizeInput(inputId = "sample","Select sample types:",
                   choices <- mySample, selected = NULL,
                   multiple = F)
  })
  
  ## Select Section-genes ####  
  output$gene <- renderUI({
    selectizeInput(inputId = "gene", label= "Select cell surface marker genes for the dot plots:",
                   choices = genelist,
                   selected = NULL,
                   multiple = TRUE,
                   options = list(
                     placeholder = 'Type to search for gene',
                     onInitialize = I('function() { this.setValue(""); }')
                   )
    )
  })
  
  observeEvent(input$cellMark,{
    if (input$cellMark){
      shinyjs::disable(id="gene")  
      shinyjs::enable(id="celltype")  
    } 
    else{
      shinyjs::enable(id="gene") 
      shinyjs::disable(id="celltype") 
    } 
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
  
  
  sample_types <- reactive({input$sample})
  
  genes <- reactive({ 
    paste(input$gene[1:length(input$gene)])
  })
  
  cell_chosen <- eventReactive(input$dplot, { 
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
      filter(marker_gene_table$gene %in% genelist)
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
  
  
  ## Violin Plot #### 
  ### home ####
  
  v_h1 <-eventReactive(input$go,{
      VlnPlot(ob_reactive(),marker_list()[1], split.by = "disease", group.by = "annotation.l2", cols=cols(),
              sort = TRUE,pt.size = 0, combine = FALSE)
  })
  v_h2 <-eventReactive(input$go,{
     VlnPlot(ob_reactive(),marker_list()[2], split.by = "disease", group.by = "annotation.l2", cols=cols(),
              sort = TRUE,pt.size = 0, combine = FALSE)
  })
  v_h3 <-eventReactive(input$go,{
    VlnPlot(ob_reactive(),marker_list()[3], split.by = "disease", group.by = "annotation.l2", cols=cols(),
              sort = TRUE,pt.size = 0, combine = FALSE)
  })
  v_h4 <-eventReactive(input$go,{
    VlnPlot(ob_reactive(),marker_list()[4], split.by = "disease", group.by = "annotation.l2", cols=cols(),
              sort = TRUE,pt.size = 0, combine = FALSE)
  })
 
  hvlist <- reactive(list(v_h1(),v_h2(),v_h3(),v_h4()))
  lapply(1:4, function(i) {
    outputId <- paste0("plotv_h", i)
    output[[outputId]] <- renderPlot(hvlist()[[i]])  
  })
  
  #Level1: CD45 (PTPRC) 
  output$plotv1 <- renderPlotly({
    ob_selected <- ob_reactive_split()[[sample_types()]]
    VlnPlot(ob_selected, marker_list()[1], split.by = "disease", group.by = "annotation.l2", cols=cols(),
              sort = TRUE,pt.size = 0, combine = FALSE)
    ggplotly(ggplot2::last_plot())
  }) 
  #Level2: EPCAM
  output$plotv2 <- renderPlotly({
    ob_selected <- ob_reactive_split()[[sample_types()]]
    ob_selected_1 <- subset(x = ob_selected, subset = split1 == plot_marker()[2])
      
      VlnPlot(ob_selected_1, marker_list()[2], split.by = "disease", group.by = "annotation.l2", cols=cols(),
              sort = TRUE, pt.size = 0, combine = FALSE)
      ggplotly(ggplot2::last_plot())
  }) 
  
  
  #Level3: PECAM1 (CD31)
  output$plotv3 <- renderPlotly({
    ob_selected <- ob_reactive_split()[[sample_types()]]
    ob_selected_1 <- subset(x = ob_selected, subset = split2 == plot_marker()[4])
      
      VlnPlot(ob_selected_1, marker_list()[3], split.by = "disease", group.by = "annotation.l2", cols=cols(),
              sort = TRUE,pt.size = 0, combine = FALSE)
      ggplotly(ggplot2::last_plot()) 
  }) 
  
  #Level4: MME (CD10)
  output$plotv4 <- renderPlotly({
    ob_selected <- ob_reactive_split()[[sample_types()]]
    ob_selected_1 <- subset(x = ob_selected, subset = split3 == plot_marker()[6])
      
      VlnPlot(ob_selected_1, marker_list()[4], split.by = "disease", group.by = "annotation.l2", cols=cols(),
              sort = TRUE,pt.size = 0, combine = FALSE)
      ggplotly(ggplot2::last_plot()) 
   
  }) 
  
  ## Dot plot ####
  
  d1 <- eventReactive(input$dplot,ignoreInit = T,{
    ob_selected <- ob_reactive_split()[[sample_types()]]
    
    ob1_1 <- subset(x = ob_selected, subset = split1 == plot_marker()[1])
      if(input$cellMark){
        #cell_types1 <- unique(sort(ob1_1@meta.data[["annotation.l2"]]))
        #marker_gene_table <- marker_gene_table()
        #marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_types1,]
        #gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
        
        marker_gene_table <- marker_gene_table()
        marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
        gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
                                                                                            
      DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
        RotatedAxis()+ aes_list() 
      }else{
        DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
          RotatedAxis()+ aes_list()
        }
  })
  
  d2 <- eventReactive(input$dplot,{
    ob_selected <- ob_reactive_split()[[sample_types()]]
    
    ob1_1 <- subset(x = ob_selected, subset = split2 == plot_marker()[3])
      if(input$cellMark){ 
        marker_gene_table <- marker_gene_table()
        marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
        gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
     
      DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
        RotatedAxis()+ aes_list()
      }else{
        DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
          RotatedAxis()+ aes_list()
        }
  })
  
  d3 <- eventReactive(input$dplot,{
    ob_selected <- ob_reactive_split()[[sample_types()]]
    
    ob1_1 <- subset(x = ob_selected, subset = split3 == plot_marker()[5])
      if(input$cellMark){
        marker_gene_table <- marker_gene_table()
        marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
        gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
       
        DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
          RotatedAxis()+ aes_list()
      }else{
        DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
          RotatedAxis()+ aes_list()
        }
  })
  d4 <- eventReactive(input$dplot,{
    ob_selected <- ob_reactive_split()[[sample_types()]]
    
    ob1_1 <- subset(x = ob_selected,subset = split4 == plot_marker()[7])
      
      if(input$cellMark){
        marker_gene_table <- marker_gene_table()
        marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
        gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
        
         DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
           RotatedAxis()+aes_list()
      }else{
        DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
          RotatedAxis()+ aes_list()
      }
  })
  
  d5 <- eventReactive(input$dplot,{
    ob_selected <- ob_reactive_split()[[sample_types()]]
    
    ob1_1 <- subset(x = ob_selected, subset = split4 == plot_marker()[8])
      
      if(input$cellMark){
        marker_gene_table <- marker_gene_table()
        marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
        gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
        
        DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
          RotatedAxis()+ aes_list()
      }else{
        DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
          RotatedAxis()+aes_list()
      }
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
}


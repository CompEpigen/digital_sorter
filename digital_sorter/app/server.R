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
server <- function(input, output, session){

  
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
                   selected = "song_2019")
  })

  output$cohort1 <- renderUI({
    if(input$cancer == "Lung"){ 
      myDataset <-as.character(names(datasets))
      } 
    else{myDataset <- c("t.b.c.")} 
    
    selectizeInput(inputId = "cohort1", label= "Select dataset of interest:", choices = myDataset, 
                   selected = myDataset[1])
  })
  
  output$cohort2 <- renderUI({
    if(input$cancer == "Lung"){ 
      myDataset2 <-as.character(names(datasets))
    } 
    else{myDataset <- c("t.b.c.")} 
    
    selectizeInput(inputId = "cohort2", label= "Select dataset of interest:", choices = myDataset2, 
                   selected = myDataset2[1])
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
                   choices <- mySample, selected = mySample[1],
                   multiple = F)
  })
  
  output$sample_se2 <- renderUI({
    mySample2 <- list_samples_disease()[[input$cohort2]]
    
    selectizeInput(inputId = "sample2","Select sample types:",
                   choices <- mySample2, selected = mySample2[1],
                   multiple = F)
  })
  
  ## Select Section-genes ####  
  output$gene <- renderUI({
    selectizeInput(inputId = "gene", label= "Select cell surface marker genes for the dot plots:",
                   choices = genelist(),
                   selected = "ABCA3",
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


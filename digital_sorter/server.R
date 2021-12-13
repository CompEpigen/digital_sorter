library(shiny)
library(shinybusy)
library(Seurat)
#library(SeuratDisk) #for saving file
library(dplyr)
library(reshape2)
library(ggplot2)
library(plotly)
library(grid)
library(gridExtra)
library(shinyjs)
library(cerebroApp)


source("R/util.R")


#server.R
server <- function(input, output, update_gene_list=F) {
  ## update cell surface marker genes or not
  if(update_gene_list){
    genelist <- get_cell_markers()
  }
  shinybusy::show_modal_spinner(text = "Loading datasets...") # show the modal window
  
  if(F){ 
    ## Create split object (split by datasets)
    # use separated R script 20211208_merge_seuratobject_inspect.R
    dataset <- readRDS("rawdata/lung_mapped_cellxgene_fixed_addsplits.rds")
    split <- as.character(unique(sort(dataset@meta.data$dataset_origin)))
    ob_all <- SplitObject(dataset, split.by = "dataset_origin")
    saveRDS(ob_all,"rawdata/lung_mapped_cellxgene_seuratsplit_fixed_7datasets.rds")
  }
  ob_all <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/lung_mapped_cellxgene_seuratsplit_fixed_7datasets.rds")
  
  list_samples_disease <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/LIST_lung_seuratsplit_7datasets_samples_disease.rds")
  
  
  genelist <-readRDS("R/cell.surface.marker.rds")
  #cell-surface marker list 
  cols = c("steelblue","darkred","gold","violet","blue","yellowgreen","darkgreen","coral2","pink","orange1")
  
  shinybusy::remove_modal_spinner() # remove it when done
  
  ##interactive markers####  
  output$master_markers <- renderUI({
    if(input$cancer == "Lung"){ myMarkers <-  c("PTPRC","EPCAM","PECAM1","MME")} 
    else{myMarkers <- c("t.b.c.")} 
    selectizeInput(inputId = "marker", label= "Select master markers (default selected):",
                   choices = myMarkers,
                   selected = myMarkers,
                   multiple = TRUE
    )
  })
  
  ## select cohort #### 
  output$cohort <- renderUI({
    if(input$cancer == "Lung"){ 
      #myCohort <-as.character(unique(sort(ob_all@meta.data$dataset_origin)))
      myCohort <-c("bischoff_2021","kim_2020","lambrechts_2018","song_2019","travaglini_2020","wu_2021","zilionis_2019")
      } 
    else{myCohort <- c("t.b.c.")} 
   
    selectizeInput(inputId = "cohort", label= "Select dataset of interest:", choices = myCohort, 
                   selected = NULL)
  })
  
  output$cohort2 <- renderUI({
    if(input$cancer == "Lung"){ 
      myDataset <-c("bischoff_2021","kim_2020","lambrechts_2018","song_2019","travaglini_2020","wu_2021","zilionis_2019")
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
    ob_reactive <- ob_all[[dataset_reactive()]]
    return(ob_reactive)
  })
  #for plots in results section
  ob_reactive2 <- reactive({
    ob_reactive2 <- ob_all[[dataset_reactive2()]]
    
    return(ob_reactive2)
  })
  
  #for plots in results section
  ob_reactive_split<- reactive({
    ob_reactive_split <- SplitObject(ob_reactive2(), split.by = "disease")
    
    return(ob_reactive_split)
  })
  ## select sample ####  
  output$sample_se <- renderUI({
    mySample <- list_samples_disease[[dataset_reactive2()]]
    
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
  
  celltypes <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/digital_sorter/R/cell_types.rds")
  
  output$celltype <- renderUI({
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
    #should also be reactive (by create a text file and the following vectors could be extracted from the file)
    #or should be user defined (haven't work on it)
    #marker_list <- eventReactive(input$cancer == "Lung", c("PTPRC","EPCAM","PECAM1","MME"))
      paste(input$marker[1:length(input$marker)])
  })
  
  plot_marker <- eventReactive(input$cancer == "Lung", c("PTPRC+","PTPRC-","EPCAM+","EPCAM-","PECAM1+","PECAM1-","MME+","MME-"))
  
  
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
    
    subset_marker_gene <- readRDS(paste0("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata_addsplits/lung_mapped_cellxgene_fixed_addsplits/",dataset_reactive2(),"/",sample_types(),"_marker_gene_table_subset_arranged50_exp.rds"))
    shinybusy::remove_modal_spinner()
    return(subset_marker_gene)
  })
  
  
  output$dotplot_title1 <- renderText({
    if(input$cancer == "Lung"){plot_marker()[1]}else {"to be continued"}
  })
  output$dotplot_title2 <- renderText({
    if(input$cancer == "Lung"){plot_marker()[c(2,3)]}else {"to be continued"}
  })
  output$dotplot_title3 <- renderText({
    if(input$cancer == "Lung"){plot_marker()[c(2,4,5)]}else {"to be continued"}
  })
  output$dotplot_title4 <- renderText({
    if(input$cancer == "Lung"){plot_marker()[c(2,4,6,7)]}else {"to be continued"}
  })
  output$dotplot_title5 <- renderText({
    if(input$cancer == "Lung"){plot_marker()[c(2,4,6,8)]}else {"to be continued"}
  })
  
  ## Violin Plot #### 
  ### home ####
  
  v_h1<-eventReactive(input$go,{
      VlnPlot(ob_reactive(),marker_list()[1], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
  })
  v_h2 <-eventReactive(input$go,{
     VlnPlot(ob_reactive(),marker_list()[2], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
  })
  v_h3 <-eventReactive(input$go,{
    VlnPlot(ob_reactive(),marker_list()[3], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
  })
  v_h4 <-eventReactive(input$go,{
    VlnPlot(ob_reactive(),marker_list()[4], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
  })
  output$plotv_h1 <- renderPlot(v_h1())  
  output$plotv_h2 <- renderPlot(v_h2()) 
  output$plotv_h3 <- renderPlot(v_h3())
  output$plotv_h4 <- renderPlot(v_h4()) 
  
  
  #Level1: CD45 (PTPRC) 
  output$plotv1 <- renderPlotly({
    ob_selected <- ob_reactive_split()[[sample_types()]]
    VlnPlot(ob_selected, marker_list()[1], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
    ggplotly(ggplot2::last_plot())
  }) 
  #Level2: EPCAM
  output$plotv2 <- renderPlotly({
    ob_selected <- ob_reactive_split()[[sample_types()]]
    ob_selected_1 <- subset(x = ob_selected, subset = split1 == plot_marker()[2])
      
      VlnPlot(ob_selected_1, marker_list()[2], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE, pt.size = 0, combine = FALSE)
      ggplotly(ggplot2::last_plot())
  }) 
  
  
  #Level3: PECAM1 (CD31)
  output$plotv3 <- renderPlotly({
    ob_selected <- ob_reactive_split()[[sample_types()]]
    ob_selected_1 <- subset(x = ob_selected, subset = split2 == plot_marker()[4])
      
      VlnPlot(ob_selected_1, marker_list()[3], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      ggplotly(ggplot2::last_plot()) 
  }) 
  
  #Level4: MME (CD10)
  output$plotv4 <- renderPlotly({
    ob_selected <- ob_reactive_split()[[sample_types()]]
    ob_selected_1 <- subset(x = ob_selected, subset = split3 == plot_marker()[6])
      
      VlnPlot(ob_selected_1, marker_list()[4], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      ggplotly(ggplot2::last_plot()) 
   
  }) 
  
  ## Dot plot ####
  
  ## draw dot plot and save as png?
  
  d1 <- eventReactive(input$dplot,ignoreInit = T,{
    ob_selected <- ob_reactive_split()[[sample_types()]]
    ob_selected@meta.data[["annotation.l2"]] <- factor(ob_selected@meta.data[["annotation.l2"]], 
                                                        levels=c("Natural Killer","Natural Killer T","Proliferating NK/T","CD8+ Memory/Effector T","CD8+ Naive T",
                                                                 "CD4+ Memory/Effector T","CD4+ Naive T","B", "Basophil/Mast 1","Basophil/Mast 2", "Plasmacytoid Dendritic",
                                                                 "Classical Monocyte","OLR1+ Classical Monocyte","Intermediate Monocyte","Nonclassical Monocyte",
                                                                 "TREM2+ Dendritic","EREG+ Dendritic","IGSF21+ Dendritic","Myeloid Dendritic Type 1","Myeloid Dendritic Type 2",
                                                                 "Macrophage","Proliferating Macrophage"))
    ob1_1 <- subset(x = ob_selected, subset = split1 == plot_marker()[1])
      if(input$cellMark){
        cell_types1 <- unique(sort(ob1_1@meta.data[["annotation.l2"]]))
        marker_gene_table <- marker_gene_table()
        marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_types1,]
        gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
        
          marker_gene_table <- marker_gene_table()
          marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
          gene_dge <- unique(marker_gene_table_sub$gene)
                                                                                            
      DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
        RotatedAxis()+ theme(legend.text=element_text(size=12),
                        axis.text=element_text(size=12),
                        axis.title=element_text(size=14),
                        legend.title=element_text(size=12))
      }else{
        DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
          RotatedAxis()+ theme(legend.text=element_text(size=12),
                               axis.text=element_text(size=12),
                               axis.title=element_text(size=14),
                               legend.title=element_text(size=12))
        }
  })
  
  d2 <- eventReactive(input$dplot,{
    ob_selected <- ob_reactive_split()[[sample_types()]]
    ob_selected@meta.data[["annotation.l2"]] <- factor(ob_selected@meta.data[["annotation.l2"]], 
                                                       levels=c("Signaling Alveolar Epithelial Type 2", "Alveolar Epithelial Type 2", 
                                                                "Club" ,"Alveolar Epithelial Type 1",
                                                                "Mucous" , "Goblet" ,"Basal", "Differentiating Basal", "Proximal Basal" ,"Proliferating Basal", 
                                                                "Ionocyte" ,"Neuroendocrine" ,
                                                                "Ciliated", "Proximal Ciliated" ,"Serous"))
    ob1_1 <- subset(x = ob_selected, subset = split2 == plot_marker()[3])
      if(input$cellMark){ 
      cell_types2 <- unique(sort(ob1_1@meta.data[["annotation.l2"]]))
      marker_gene_table <- marker_gene_table()
      marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_types2,]
      gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
      
        marker_gene_table <- marker_gene_table()
        marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
        gene_dge <- unique(marker_gene_table_sub$gene)
     
      DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
        RotatedAxis()+ theme(legend.text=element_text(size=12),
                             axis.text=element_text(size=12),
                             axis.title=element_text(size=14),
                             legend.title=element_text(size=12))
      }else{
        DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
          RotatedAxis()+ theme(legend.text=element_text(size=12),
                               axis.text=element_text(size=12),
                               axis.title=element_text(size=14),
                               legend.title=element_text(size=12))
        }
  })
  
  d3 <- eventReactive(input$dplot,{
    ob_selected <- ob_reactive_split()[[sample_types()]]
    ob_selected@meta.data[["annotation.l2"]] <- factor(ob_selected@meta.data[["annotation.l2"]], 
                                                       levels=c(#Level3
                                                                "Vein","Bronchial Vessel 1","Bronchial Vessel 2","Artery","Capillary",
                                                                "Capillary Intermediate 1","Capillary Intermediate 2","Capillary Aerocyte","Lymphatic"))
    ob1_1 <- subset(x = ob_selected, subset = split3 == plot_marker()[5])
      if(input$cellMark){
        cell_types3 <- unique(sort(ob1_1@meta.data[["annotation.l2"]]))
        marker_gene_table <- marker_gene_table()
        marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_types3,]
        gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
        
          marker_gene_table <- marker_gene_table()
          marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
          gene_dge <- unique(marker_gene_table_sub$gene)
       
        DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
          RotatedAxis()+ theme(legend.text=element_text(size=12),
                               axis.text=element_text(size=12),
                               axis.title=element_text(size=14),
                               legend.title=element_text(size=12))
      }else{
        DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
          RotatedAxis()+ theme(legend.text=element_text(size=12),
                               axis.text=element_text(size=12),
                               axis.title=element_text(size=14),
                               legend.title=element_text(size=12))
        }
  })
  d4 <- eventReactive(input$dplot,{
    ob_selected <- ob_reactive_split()[[sample_types()]]
    ob_selected@meta.data[["annotation.l2"]] <- factor(ob_selected@meta.data[["annotation.l2"]], 
                                                       levels=c(#Level 4&5
                                                                "Platelet/Megakaryocyte","Adventitial Fibroblast","Alveolar Fibroblast"))
    ob1_1 <- subset(x = ob_selected,subset = split4 == plot_marker()[7])
      
      if(input$cellMark){
        cell_types4 <- unique(sort(ob1_1@meta.data[["annotation.l2"]]))
        marker_gene_table <- marker_gene_table()
        marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_types4,]
        
        gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
        
          marker_gene_table <- marker_gene_table()
          marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
          gene_dge <- unique(marker_gene_table_sub$gene)
        
         DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
           RotatedAxis()+ theme(legend.text=element_text(size=12),
                                axis.text=element_text(size=12),
                                axis.title=element_text(size=14),
                                legend.title=element_text(size=12))
      }else{
        DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
          RotatedAxis()+ theme(legend.text=element_text(size=12),
                               axis.text=element_text(size=12),
                               axis.title=element_text(size=14),
                               legend.title=element_text(size=12))
      }
  })
  
  d5 <- eventReactive(input$dplot,{
    ob_selected <- ob_reactive_split()[[sample_types()]]
    ob_selected@meta.data[["annotation.l2"]] <- factor(ob_selected@meta.data[["annotation.l2"]], 
                                                       levels=c("Mesothelial","Lipofibroblast",
                                                                "Myofibroblast","Fibromyocyte","Airway Smooth Muscle","Vascular Smooth Muscle","Pericyte","Plasma"))
    ob1_1 <- subset(x = ob_selected, subset = split4 == plot_marker()[8])
      
      if(input$cellMark){
        cell_types5 <- unique(sort(ob1_1@meta.data[["annotation.l2"]]))
        marker_gene_table <- marker_gene_table()
        marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_types5,]
        gene_dge <- unique(marker_gene_table_sub$gene)[1:10]
        
          marker_gene_table <- marker_gene_table()
          marker_gene_table_sub <- marker_gene_table[marker_gene_table$annotation.l2 %in% cell_chosen(),]
          gene_dge <- unique(marker_gene_table_sub$gene)
       
        DotPlot(ob1_1, features = gene_dge, group.by = "annotation.l2") + 
          RotatedAxis()+ theme(legend.text=element_text(size=12),
                               axis.text=element_text(size=12),
                               axis.title=element_text(size=14),
                               legend.title=element_text(size=12))
      }else{
        DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + 
          RotatedAxis()+ theme(legend.text=element_text(size=12),
                               axis.text=element_text(size=12),
                               axis.title=element_text(size=14),
                               legend.title=element_text(size=12))
      }
  })
  
  output$plotd1 <- renderPlot(d1())
  output$plotd2 <- renderPlot(d2())
  output$plotd3 <- renderPlot(d3())
  output$plotd4 <- renderPlot(d4())
  output$plotd5 <- renderPlot(d5())

  
  ## table next to the dot plot #### 
  output$table1 <- renderDataTable(d1()[["data"]][, c(3,4,1,2,5)], 
                                     options = list(pageLength = 5)
  )
  output$table2 <- renderDataTable(d2()[["data"]][, c(3,4,1,2,5)], 
                                   options = list(pageLength = 5)
  )
  output$table3 <- renderDataTable(d3()[["data"]][, c(3,4,1,2,5)], 
                                   options = list(pageLength = 5)
  )
  output$table4 <- renderDataTable(d4()[["data"]][, c(3,4,1,2,5)], 
                                   options = list(pageLength = 5)
  )
  output$table5 <- renderDataTable(d5()[["data"]][, c(3,4,1,2,5)], 
                                   options = list(pageLength = 5)
  )
  
  
}


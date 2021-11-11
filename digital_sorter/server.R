library(shiny)
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


source("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/util.R")
#To-do: loading a data automatically (for loop?)
song_2019 <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/song_2019_addsplits.rds")
#nsclc <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/nsclc_primary.rds")

song_split <- SplitObject(song_2019, split.by = "disease")

genelist <-readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/cell.surface.marker.rds") 
#cell-surface marker list #different from cerebroApp????????????????????

# CD45 (PTPRC)
# CD31 (PECAM1)
cols = c("steelblue","darkred","gold","coral2")

#server.R
server <- function(input, output, update_gene_list=F) {
  if(update_gene_list){
    genelist <- get_cell_markers()
  }
  ##interactive markers####  
  output$master_markers <- renderUI({
    if(input$cancer == "Lung") myMarkers <-  c("PTPRC","EPCAM","PECAM1","MME")
    
    else myMarkers <- c("t.b.c.")
    selectizeInput(inputId = "marker", label= "Select master markers (default selected):",
                   choices = myMarkers,
                   selected = myMarkers,
                   multiple = TRUE
    )
  })
  
  ## select cohort #### 
  output$cohort <- renderUI({
    if(input$cancer == "Lung") myCohort <-  c("song_2019","nsclc_primary") #should detect the sample itself
    
    else myCohort <- c("t.b.c.")
    
    selectizeInput(inputId = "cohort", label= "Select dataset of interest:", choices = myCohort, 
                   selected = NULL)
  })
  
  ## select sample ####  
  output$sample_se <- renderUI({
    if(input$cohort == "song_2019") myChoices <-  c("adjacent normal","nsclc") #should detect the sample itself
    else if(input$cohort == "nsclc_primary") myChoices <-  c("nsclc")
    else myChoices <- c("sample1","sample2","sample3")
    
    selectizeInput(inputId = "sample","Select sample types:",
                   choices <- myChoices, selected = myChoices,
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
  
  observeEvent(input$testDGE,{
    if (input$testDGE) shinyjs::disable(id="gene")  
    else shinyjs::enable(id="gene") 
  })
  
  marker_gene_table <- eventReactive(input$testDGE,{
    
      seurat <- getMarkerGenes(
        song_2019,
        assay = 'RNA',
        organism = 'hg',
        groups = c('annotation.l2'),
        name = 'cerebro_seurat',
        only_pos = TRUE,
        min_pct = 0.7,
        thresh_logFC = 0.25,
        thresh_p_val = 0.01,
        test = 'wilcox',
        verbose = TRUE
      )
      marker_gene <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]][["gene"]]
      marker_gene_table <- seurat@misc[["marker_genes"]][["cerebro_seurat"]][["annotation.l2"]]
      marker_gene_table <- marker_gene_table %>% filter(marker_gene_table$on_cell_surface =="TRUE")%>%
        arrange(desc(avg_log2FC))
      return(marker_gene_table)
  })
  

  
  
  ## Event reactive objects #### 
  datasets <- eventReactive(input$go, {input$cohort})
  
  #should also be reactive (by create a text file and the following vectors could be extracted from the file)
  #or should be user defined (haven't work on it)
  marker_list <- eventReactive(input$cancer == "Lung", c("PTPRC","EPCAM","PECAM1","MME"))
  plot_marker <- eventReactive(input$cancer == "Lung", c("PTPRC+","PTPRC-","EPCAM+","EPCAM-","PECAM1+","PECAM1-","MME+","MME-"))
  
  
  sample_types <- eventReactive(input$levelplot, {input$sample})
  
  genes <- eventReactive(input$levelplot, { 
    if (input$testDGE){
      #marker_gene_table <- marker_gene_table()
      #gene_dge <- unique(marker_gene_table$gene)[1:10]
      gene_dge <- c("SPARC","LY6D","CLU","TFPI","AREG","EGFL7","HYAL2","VAMP5", "C3","HLA-E") 
      return(gene_dge)
    }else    paste(input$gene[1:length(input$gene)])
  })
  
 
  
  output$dotplot_title1 <- renderText({
    if(input$cancer == "Lung"){plot_marker()[1]}else {"to be continued"}
  })
  output$dotplot_title2 <- renderText({
    if(input$cancer == "Lung"){plot_marker()[2]}else {"to be continued"}
  })
  output$dotplot_title3 <- renderText({
    if(input$cancer == "Lung"){plot_marker()[c(2,4)]}else {"to be continued"}
  })
  output$dotplot_title4 <- renderText({
    if(input$cancer == "Lung"){plot_marker()[c(2,4,6)]}else {"to be continued"}
  })
  output$dotplot_title5 <- renderText({
    if(input$cancer == "Lung"){plot_marker()[c(2,4,6,8)]}else {"to be continued"}
  })
  
  ## Violin Plot #### 
  
  
  #home
  v_h1<-eventReactive(input$go,{
    if(datasets() == "song_2019"){
      
      VlnPlot(song_2019,marker_list()[1], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      
      
    }else if(datasets() == "nsclc_primary") {#not yet finished
      VlnPlot(nsclc, marker_list()[1], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
     
    }
  })
  
  
  
  v_h2 <-eventReactive(input$go,{
    if(datasets() == "song_2019"){
      
      VlnPlot(song_2019,marker_list()[2], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      
     
    }else if(datasets() == "nsclc_primary") {#not yet finished
      VlnPlot(nsclc, marker_list()[2], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)+ theme(legend.text=element_text(size=10), 
                                                               axis.text=element_text(size=10), 
                                                               legend.title=element_text(size=10))
     
    }
  })
  v_h3 <-eventReactive(input$go,{
    if(datasets() == "song_2019"){
      
      VlnPlot(song_2019,marker_list()[3], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      
      
    }else if(datasets() == "nsclc_primary") {#not yet finished
      VlnPlot(nsclc, marker_list()[3], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)+ theme(legend.text=element_text(size=10), 
                                                               axis.text=element_text(size=10), 
                                                               legend.title=element_text(size=10))
    } 
  })
  v_h4 <-eventReactive(input$go,{
    if(datasets() == "song_2019"){
      
      VlnPlot(song_2019,marker_list()[4], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      
      
    }else if(datasets() == "nsclc_primary") {#not yet finished
      VlnPlot(nsclc, marker_list()[4], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)+ theme(legend.text=element_text(size=10), 
                                                               axis.text=element_text(size=10), 
                                                               legend.title=element_text(size=10))
    } 
  })
  output$plotv_h1 <- renderPlot(v_h1())  
  output$plotv_h2 <- renderPlot(v_h2()) 
  output$plotv_h3 <- renderPlot(v_h3())
  output$plotv_h4 <- renderPlot(v_h4()) 
  
  #Level1: CD45 (PTPRC) 
  output$plotv1 <- renderPlotly({
    if(datasets() == "song_2019"){
      if(sample_types()[1] =="nsclc"){
        song_2019_selected <- song_split[["nsclc"]]
      }else if(sample_types()[1] =="adjacent normal"){
        song_2019_selected <- song_split[["adjacent normal"]]
      }
      
      VlnPlot(song_2019_selected,marker_list()[1], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      
      ggplotly(ggplot2::last_plot())
    } 
    
    else if(datasets() == "nsclc_primary") {#not yet finished
      VlnPlot(nsclc, marker_list()[1], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      ggplotly(ggplot2::last_plot())
    }
  }) 
  #Level2: EPCAM
  output$plotv2 <- renderPlotly({
    if(datasets() == "song_2019"){
      if(sample_types()[1] =="nsclc"){
        song_2019_selected <- song_split[["nsclc"]]
      }else if(sample_types()[1] =="adjacent normal"){
        song_2019_selected <- song_split[["adjacent normal"]]
      }
      
      song_2019_selected_1 <- subset(x = song_2019_selected, subset = split1 == plot_marker()[2])
      
      VlnPlot(song_2019_selected_1,marker_list()[2], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      
      ggplotly(ggplot2::last_plot())
    } else if(datasets() == "nsclc_primary") {#not yet finished
      VlnPlot(nsclc, marker_list()[2], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      
      ggplotly(ggplot2::last_plot())
    }
  }) 
  
  
  #Level3: PECAM1 (CD31)
  output$plotv3 <- renderPlotly({
    if(datasets() == "song_2019"){
      if(sample_types()[1] =="nsclc"){
        song_2019_selected <- song_split[["nsclc"]]
      }else if(sample_types()[1] =="adjacent normal"){
        song_2019_selected <- song_split[["adjacent normal"]]
      }
      
      song_2019_selected_1 <- subset(x = song_2019_selected, subset = split2 == plot_marker()[4])
      
      VlnPlot(song_2019_selected_1,marker_list()[3], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      
      ggplotly(ggplot2::last_plot()) 
    } 
    
    else if(datasets() == "nsclc_primary") {#not yet finished
      VlnPlot(nsclc, marker_list()[1], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      ggplotly(ggplot2::last_plot()) 
    }
  }) 
  
  #Level4: MME (CD10)
  output$plotv4 <- renderPlotly({
    if(datasets() == "song_2019"){
      if(sample_types()[1] =="nsclc"){
        song_2019_selected <- song_split[["nsclc"]]
      }else if(sample_types()[1] =="adjacent normal"){
        song_2019_selected <- song_split[["adjacent normal"]]
      }
      
      song_2019_selected_1 <- subset(x = song_2019_selected, subset = split3 == plot_marker()[6])
      
      VlnPlot(song_2019_selected_1,marker_list()[4], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      
      ggplotly(ggplot2::last_plot()) 
    } 
    
    else if(datasets() == "nsclc_primary") {#not yet finished
      VlnPlot(nsclc, marker_list()[1], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      ggplotly(ggplot2::last_plot()) 
    }
  }) 
  
  ## Dot plot ####
  
  ## draw dot plot and save as png
  d1 <- eventReactive(input$levelplot,{
    if(datasets() == "song_2019"){
      
        song_2019_selected <- song_split[[sample_types()]]
        ob1_1 <- subset(x = song_2019_selected, subset = split1 == plot_marker()[1])
        DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + RotatedAxis()+ theme(legend.text=element_text(size=12),
                                                                                              axis.text=element_text(size=12),
                                                                                              axis.title=element_text(size=14),
                                                                                              legend.title=element_text(size=12))
       
    }else if(datasets() == "nsclc_primary") {#not yet finished
      
      DotPlot(nsclc$`PTPRC+`, features = genes(), group.by = "annotation.l2") + RotatedAxis()+ theme(axis.text=element_text(size=12))
      
      }
  })
  d2 <- eventReactive(input$levelplot,{
    if(datasets() == "song_2019"){
      
      song_2019_selected <- song_split[[sample_types()]]
      ob1_1 <- subset(x = song_2019_selected, subset = split1 == plot_marker()[2])
      DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + RotatedAxis()+ theme(legend.text=element_text(size=12),
                                                                                            axis.text=element_text(size=12),
                                                                                            axis.title=element_text(size=14),
                                                                                            legend.title=element_text(size=12))
      
    }else if(datasets() == "nsclc_primary") {#not yet finished
      
      DotPlot(nsclc$`PTPRC+`, features = genes(), group.by = "annotation.l2") + RotatedAxis()+ theme(axis.text=element_text(size=12))
      
    }
  })
  d3 <- eventReactive(input$levelplot,{
    if(datasets() == "song_2019"){
      
      song_2019_selected <- song_split[[sample_types()]]
      ob1_1 <- subset(x = song_2019_selected, subset = split2 == plot_marker()[4])
      DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + RotatedAxis()+ theme(legend.text=element_text(size=12),
                                                                                            axis.text=element_text(size=12),
                                                                                            axis.title=element_text(size=14),
                                                                                            legend.title=element_text(size=12))
      
    }else if(datasets() == "nsclc_primary") {#not yet finished
      
      DotPlot(nsclc$`PTPRC+`, features = genes(), group.by = "annotation.l2") + RotatedAxis()+ theme(axis.text=element_text(size=12))
      
    }
  })
  d4 <- eventReactive(input$levelplot,{
    if(datasets() == "song_2019"){
      
      song_2019_selected <- song_split[[sample_types()]]
      ob1_1 <- subset(x = song_2019_selected, subset = split3 == plot_marker()[6])
      DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + RotatedAxis()+ theme(legend.text=element_text(size=12),
                                                                                            axis.text=element_text(size=12),
                                                                                            axis.title=element_text(size=14),
                                                                                            legend.title=element_text(size=12))
      
    }else if(datasets() == "nsclc_primary") {#not yet finished
      
      DotPlot(nsclc$`PTPRC+`, features = genes(), group.by = "annotation.l2") + RotatedAxis()+ theme(axis.text=element_text(size=12))
      
    }
  })
  
  output$plotd1 <- renderPlot(d1())
  output$plotd2 <- renderPlot(d2())
  output$plotd3 <- renderPlot(d3())
  output$plotd4 <- renderPlot(d4())
  
  
  ## table next to the dot plot #### 
  output$table1 <- renderDataTable(d1()[["data"]], 
                                     options = list(pageLength = 5)
  )
  output$table2 <- renderDataTable(d2()[["data"]], 
                                   options = list(pageLength = 5)
  )
  output$table3 <- renderDataTable(d3()[["data"]], 
                                   options = list(pageLength = 5)
  )
  output$table4 <- renderDataTable(d4()[["data"]], 
                                   options = list(pageLength = 5)
  )
  
  
}



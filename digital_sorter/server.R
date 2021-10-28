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

#To-do: loading a data automatically (for loop?)
song_2019 <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/song_2019.rds")
#nsclc <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/nsclc_primary.rds")

song_split <- SplitObject(song_2019, split.by = "disease")
genelist <- c("CD8A","CD8B","CD4", "MRC1","CD14","SIGLEC1") #to be changed by cell-surface marker list

# CD45 (PTPRC)
# CD31 (PECAM1)
cols = c("steelblue","darkred","gold","coral2")

#server.R
server <- function(input, output) {
  
  ##interactive markers####  
  output$master_markers <- renderText({
    if(input$cancer == "Lung"){
      
      paste("You selected lung.","Level 1: PTPRC (CD45) +","Level 2: EPCAM +","Level 3: PECAM1 (CD31) +","Level 4: PECAM1 (CD31) -", sep="\n")
      
    }
    else {"to be continued"}
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
    selectizeInput(inputId = "gene", label= "Select genes of interest for the dot plot:",
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
  
  
  ## Event reactive objects #### 
  datasets <- eventReactive(input$go, {input$cohort})
  
  marker_list <- eventReactive(input$cancer == "Lung", c("PTPRC","EPCAM","PECAM1"))
  plot_marker <- eventReactive(input$cancer == "Lung", c("PTPRC+","EPCAM+","PECAM1+","PECAM1-"))
  
  
  sample_types <- eventReactive(input$levelplot, {input$sample})
  
  genes <- eventReactive(input$dotplot, { 
    paste(input$gene[1:length(input$gene)])
  })
  
  output$dotplot_title1 <- renderText({
    if(input$cancer == "Lung"){plot_marker()[1]}else {"to be continued"}
  })
  output$dotplot_title2 <- renderText({
    if(input$cancer == "Lung"){plot_marker()[2]}else {"to be continued"}
  })
  output$dotplot_title3 <- renderText({
    if(input$cancer == "Lung"){plot_marker()[3]}else {"to be continued"}
  })
  output$dotplot_title4 <- renderText({
    if(input$cancer == "Lung"){plot_marker()[4]}else {"to be continued"}
  })
  
  ## Violin Plot #### 
  
  
  #home
  v_h1<-eventReactive(input$go,{
    if(datasets() == "song_2019"){
      
      VlnPlot(song_2019,marker_list()[1], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      
      ggplotly(ggplot2::last_plot())
      
    }else if(datasets() == "nsclc_primary") {#not yet finished
      VlnPlot(nsclc, marker_list()[1], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      ggplotly(ggplot2::last_plot())
    }
  })
  
  
  
  v_h2 <-eventReactive(input$go,{
    if(datasets() == "song_2019"){
      
      VlnPlot(song_2019,marker_list()[2], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      
      ggplotly(ggplot2::last_plot())
      
    }else if(datasets() == "nsclc_primary") {#not yet finished
      VlnPlot(nsclc, marker_list()[2], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)+ theme(legend.text=element_text(size=10), 
                                                               axis.text=element_text(size=10), 
                                                               legend.title=element_text(size=10))
      ggplotly(ggplot2::last_plot())
    }
  })
  v_h3 <-eventReactive(input$go,{
    if(datasets() == "song_2019"){
      
      VlnPlot(song_2019,marker_list()[3], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)
      
      ggplotly(ggplot2::last_plot())
      
    }else if(datasets() == "nsclc_primary") {#not yet finished
      VlnPlot(nsclc, marker_list()[3], split.by = "disease", group.by = "annotation.l2", cols=cols,
              sort = TRUE,pt.size = 0, combine = FALSE)+ theme(legend.text=element_text(size=10), 
                                                               axis.text=element_text(size=10), 
                                                               legend.title=element_text(size=10))
      ggplotly(ggplot2::last_plot())
    } 
  })
  output$plotv_h1 <- renderPlotly(v_h1())  
  output$plotv_h2 <- renderPlotly(v_h2()) 
  output$plotv_h3 <- renderPlotly(v_h3()) 
  
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
  
  
  #Level3&4: PECAM1 (CD31)
  output$plotv3 <- renderPlotly({
    if(datasets() == "song_2019"){
      if(sample_types()[1] =="nsclc"){
        song_2019_selected <- song_split[["nsclc"]]
      }else if(sample_types()[1] =="adjacent normal"){
        song_2019_selected <- song_split[["adjacent normal"]]
      }
      
      song_2019_selected_1 <- subset(x = song_2019_selected, subset = split1 == plot_marker()[2])
      
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
  
  
  
  ## Dot plot ####
  
  #Level1
  d1 <- eventReactive(input$levelplot,{
    if(datasets() == "song_2019"){
      
      song_2019_selected <- song_split[[sample_types()]]
      ob1_1 <- subset(x = song_2019_selected, subset = split1 == plot_marker()[1])
      DotPlot(ob1_1, features = genes(), group.by = "annotation.l2") + RotatedAxis()+ theme(legend.text=element_text(size=10),
                                                                                            axis.text=element_text(size=10),
                                                                                            axis.title=element_text(size=12),
                                                                                            legend.title=element_text(size=10))
    }else if(datasets() == "nsclc_primary") {#not yet finished
      DotPlot(nsclc$`PTPRC+`, features = genes(), group.by = "annotation.l2") + RotatedAxis()+ theme(axis.text=element_text(size=10))
    }
  })
  
  output$plotd1 <- renderPlot(d1())
  
  
  
  
  ## table next to the dot plot #### 
  output$table1 <- renderDataTable(d1()[["data"]], 
                                   options = list(pageLength = 5)
  )
  
  
}



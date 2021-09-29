library(shiny)
library(Seurat)
#library(SeuratDisk) #for saving file
library(dplyr)
library(reshape2)
library(ggplot2)
library(plotly)
library(grid)
library(gridExtra)
#read data 
#To-do: loading a data automatically (for loop?)
song_2019 <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/song_2019.rds")
#nsclc <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/nsclc_primary.rds")

song_split <- SplitObject(song_2019, split.by = "disease")


# CD45 (PTPRC)
# CD31 (PECAM1)
cols = c("steelblue","darkred","gold","coral2")

#server.R
server <- function(input, output) {
  
  output$sample <- renderUI({
    if(input$cohort == "song_2019") myChoices <-  c("adjacent normal","nsclc") #should detect the sample itself
    else if(input$cohort == "nsclc_primary") myChoices <-  c("nsclc")
    else myChoices <- c("sample1","sample2","sample3")
    
    selectizeInput(inputId = "sample","Select sample types",
                 choices <- myChoices, selected = myChoices,
                 multiple = TRUE)
  })
  
  features <- eventReactive(input$go, {input$marker})
  datasets <- eventReactive(input$go, {input$cohort})
  sample_types <- eventReactive(input$go, {input$sample})

  
  data <- eventReactive(input$go, { 
    paste("You chose", input$cancer,"-", input$cohort,"-")
  })
  data2 <- eventReactive(input$go, { 
    paste(input$sample[1:length(input$sample)], " ")
  })
  
  genes <- eventReactive(input$go, { 
    paste(input$gene[1:length(input$gene)])
  })
 
  output$result1 <- renderText({
  data()
  })
  output$result2 <- renderText({
    data2()
  })
  output$result3 <- renderText({
    genes()
  })
  
  
  output$plotv <- renderPlotly({
    if(datasets() == "song_2019"){
      if(length(sample_types())==1 & sample_types()[1] =="nsclc"){
        song_2019_selected <- song_split[["nsclc"]]
      }else if(length(sample_types())==1 & sample_types()[1] =="adjacent normal"){
        song_2019_selected <- song_split[["adjacent normal"]]
      }else{
        song_2019_selected <- song_2019
      }
      
      VlnPlot(song_2019_selected,features(), split.by = "disease", group.by = "annotation.l2", cols=cols,
                                            sort = TRUE,pt.size = 0, combine = FALSE)
     
      ggplotly(ggplot2::last_plot())
    } 
    
    else if(datasets() == "nsclc_primary") {#not yet finished
      VlnPlot(nsclc, features(), split.by = "disease", group.by = "annotation.l2", cols=cols,
                                                       sort = TRUE,pt.size = 0, combine = FALSE)
      ggplotly(ggplot2::last_plot())
    }
   
    
  }) 
 
  
  output$plotd1 <- renderPlotly({
    if(datasets() == "song_2019"){
      
        ob1 <- SplitObject(song_split[[1]], split.by = "split1")
       
        DotPlot(ob1$`PTPRC+`, features = features(), group.by = "annotation.l2") + RotatedAxis()
        ggplotly(ggplot2::last_plot())
      
    }else if(datasets() == "nsclc_primary") {#not yet finished
      DotPlot(nsclc$`PTPRC+`, features = features(), group.by = "annotation.l2") + RotatedAxis()
      ggplotly(ggplot2::last_plot())
      }
  })
  output$plotd2 <- renderPlotly({
    if(datasets() == "song_2019"){
      
      ob1 <- SplitObject(song_split[[2]], split.by = "split1")
      
      DotPlot(ob1$`PTPRC+`, features = features(), group.by = "annotation.l2") + RotatedAxis()
      ggplotly(ggplot2::last_plot())
      
    }else if(datasets() == "nsclc_primary") {#not yet finished
      DotPlot(nsclc$`PTPRC+`, features = features(), group.by = "annotation.l2") + RotatedAxis()
      ggplotly(ggplot2::last_plot())
    }
  })
    
  
}

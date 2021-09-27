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
song_2019 <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/song_2019.rds")
nsclc <- readRDS("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter/rawdata/nsclc_primary.rds")
kim_2020 <- subset(x = nsclc, subset = dataset_origin == "kim_2020")

#print(table(song_2019@meta.data$dataset_origin))
#print(table(song_2019@meta.data$disease))
#print(table(nsclc@meta.data$dataset_origin))
#print(table(nsclc@meta.data$disease))
#print(table(kim_2020@meta.data$dataset_origin))
#print(table(kim_2020@meta.data$disease))


# CD45 (PTPRC)
# CD31 (PECAM1)
cols = c("steelblue","darkred","gold","coral2")

#server.R
server <- function(input, output) {
  
  output$sample <- renderUI({
    if(input$cohort == "song_2019") myChoices <-  c("adjacent normal","nsclc") 
    else if(input$cohort == "kim_2020") myChoices <-  c("nsclc")
    else myChoices <- c("sample1","sample2","sample3")
    
    selectizeInput(inputId = "sample","Select disease",
                 choices <- myChoices, selected = myChoices,
                 multiple = TRUE)
  })
  
  features <- eventReactive(input$go, {input$marker})
  datasets <- eventReactive(input$go, {input$cohort})
  
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
      VlnPlot(song_2019, features(), split.by = "disease", group.by = "annotation.l2", cols=cols,
                                            sort = TRUE,pt.size = 0, combine = FALSE)
     
      ggplotly(ggplot2::last_plot())
    } 
    
    else if(datasets() == "travaglini_2020") {
      VlnPlot(nsclc, features(), split.by = "disease", group.by = "annotation.l2", cols=cols,
                                                       sort = TRUE,pt.size = 0, combine = FALSE)
      ggplotly(ggplot2::last_plot())
    }
    
    else if(datasets() == "kim_2020"){
      VlnPlot(kim_2020, features(), split.by = "disease", group.by = "annotation.l2", cols=cols,
                             sort = TRUE,pt.size = 0, combine = FALSE)
      ggplotly(ggplot2::last_plot())        
      
    } 
    
  }) 
  
  
}

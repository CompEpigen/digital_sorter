library(shiny)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(reshape2)
library(ggplot2)
library(plotly)
library(grid)
library(gridExtra)
#read data
song_2019 <- readRDS("D:/Uni Heidleberg/DKFZ/semester 5/thesis/digital_cell_sorter/rawdata/song_2019.rds")

table(song_2019@meta.data$disease)
table(song_2019@meta.data$donor)
# CD45 (PTPRC)
#features=c("PTPRC")  
cols = c("steelblue","darkred","gold","coral2")

#server.R
server <- function(input, output) {
  
  output$sample <- renderUI({
    if(input$cohort == "song_2019") myChoices <-  c("adjacent normal","nsclc") 
    else if(input$cohort == "travaglini_2020") myChoices <-  c("sample1","sample2")
    else myChoices <- c("sample1","sample2","sample3")
    
    selectizeInput(inputId = "sample","Select disease",
                 choices <- myChoices, selected = myChoices,
                 multiple = TRUE)
  })
  
  features <- eventReactive(input$go, {input$marker})
 
  
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
  
  output$violin <- renderPlot(
    if(input$cohort == "song_2019") ggplotly(VlnPlot(song_2019, features(), split.by = "disease", group.by = "annotation.l2", cols=cols,
                                            sort = TRUE,pt.size = 0, combine = FALSE)
                                            )
    else if(input$cohort == "travaglini_2020") ggplotly(VlnPlot(travaglini_2020, features(), split.by = "disease", group.by = "annotation.l2", cols=cols,
                                                       sort = TRUE,pt.size = 0, combine = FALSE))
    else  ggplotly(VlnPlot(kim_2020, features(), split.by = "disease", group.by = "annotation.l2", cols=cols,
                 sort = TRUE,pt.size = 0, combine = FALSE))
    
  )
}

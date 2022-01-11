library(shiny)
library(shinybusy)
library(shinythemes)
library(Seurat)
#library(SeuratDisk) #for saving file
library(dplyr)
library(reshape2)
library(ggplot2)
library(plotly)
library(grid)
library(gridExtra)
library(shinyjs)
library(shinydashboard)
library(shinyWidgets)
library(cerebroApp)


library(profvis)

setwd("/omics/groups/OE0219/internal/Jessie_2021/P01.digital_sorter")


profvis({
runApp("./digital_sorter")
  
})


library(shiny)
library(shinythemes)
library(Seurat)
#library(SeuratDisk) #for saving file
library(dplyr)
library(reshape2)
library(ggplot2)
library(plotly)
library(grid)
library(gridExtra)

library(profvis)

profvis({
runApp("./digital_sorter")
  
})


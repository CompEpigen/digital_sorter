library(shiny)
library(shinybusy)
library(shinythemes)
library(Seurat)
#library(SeuratDisk) #for saving file
library(dplyr)
library(reshape2)
library(ggplot2)
library(plotly)
library(gridExtra)
library(shinyjs)
library(shinydashboard)
library(shinyWidgets)


#library(profvis)

setwd("Path to the folder of digital sorter")


#profvis({

  runApp("./digital_sorter")
  
#})


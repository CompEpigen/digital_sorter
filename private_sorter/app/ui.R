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
#library(shinyauthr)

library(DT)






header <- dashboardHeader( title =tags$b(tags$img(src = "digital_sorter.png", width=80, height=50),
                                         " Private Sorter"),
                           disable = FALSE, 
                           titleWidth  = 300,
                           uiOutput("logoutbtn"))

sidebar <- dashboardSidebar(width = 300,
                            useShinyjs(),
                            uiOutput("menu")
                            
  ) 


body <- dashboardBody( uiOutput("body"))
ui<-dashboardPage(header, sidebar, body, skin = "blue",
                  title= "Digital Sorter")

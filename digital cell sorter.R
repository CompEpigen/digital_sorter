
library(shiny)


#Import Data
 genelist <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 
               'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z')


# Define UI

ui <- fluidPage(
  
  h1("Digital Cell Sorter"), 
  p(style = "font-family:Impact"),
  tags$hr(),
  
  
  sidebarLayout(
    
    
    sidebarPanel(
      titlePanel("Input"),
      
      #gene
      selectizeInput(inputId = "gene", label= "Select genes",
                     choices = genelist,
                     selected = NULL,
                     multiple = TRUE,
                     options = list(
                       placeholder = 'Type to search for gene',
                       onInitialize = I('function() { this.setValue(""); }')
                     )
      ),
     
      #cancer-lung
      checkboxGroupInput(inputId ="cancer", label= "Select cancer types", choices = "Lung", selected = "Lung"),
      #markers
      conditionalPanel('input.cancer == "Lung"', 
                       selectizeInput(inputId = "marker", 
                                      label= "master markers",
                                      multiple = TRUE, 
                                      choices = c("CD45","EPCAM","PECAM1"),
                                      selected = c("CD45","EPCAM","PECAM1"))
                       ),
                     
      
      selectizeInput(
        inputId = "cohort", label= "Select dataset of interest", choices = c("cohort1","cohort2","cohort3"), 
        selected = NULL), 
      
      uiOutput("sample"),
                       
      actionButton("go", "Go")
      
    ),#side panel end
    
    mainPanel(
      fluidRow(
        textOutput("result1"),textOutput("result2"),
        tags$hr(),
        print("You chose gene(s) "),
        textOutput("result3")
      )
      
    )#main panel end
     
  )#side bar layout end
      
) #page end

  


server <- function(input, output) {
  
  output$sample <- renderUI({
    if(input$cohort == "cohort1") myChoices <-  c("sample1","sample2","sample3","sample4")
    else if(input$cohort == "cohort2") myChoices <-  c("sample1","sample2")
    else myChoices <- c("sample1","sample2","sample3")
    
    selectizeInput(inputId = "sample","Select sample of interest",
                 choices <- myChoices, selected = NULL,
                 multiple = TRUE)
  })
  

  
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
  

}

shinyApp(ui = ui, server = server)
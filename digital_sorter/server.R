
library(shiny)

#server.R
function(input, output) {
  
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
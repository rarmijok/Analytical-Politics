
library(shiny)
load("cep702.rda")
namescep702 <- names(cep702)


shinyUI(pageWithSidebar(
  headerPanel(" "),
  sidebarPanel(
    selectInput('xcol', 'Variable Demográfica', namescep702[-1])
               ),
  mainPanel(
    plotOutput('newHist'),
    h6('Fuente: Encuesta Cep Septiembre-Octubre 2013')
  )
))
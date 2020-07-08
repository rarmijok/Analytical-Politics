
library(shiny)
load("cep61.Rda")
namescep611 <- names(cep611)


shinyUI(pageWithSidebar(
  headerPanel(" "),
  sidebarPanel(
    selectInput('xcol', 'Pregunta', namescep611[-c(1,2)])
  ),
  mainPanel(
    plotOutput('newHist'),
    h6('Fuente: Encuesta CEP Octubre 2009')
  )
))


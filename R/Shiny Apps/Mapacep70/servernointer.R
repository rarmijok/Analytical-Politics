
library(shiny)
load("cep70.rda")
library(ggplot2)

shinyServer(
  function(input, output) {
    
    output$newHist <- reactivePlot( function() {
      plot <- ggplot(cep70, aes(x=te2p05, fill=ddp31)) + geom_density(alpha=.3) + theme_bw()
      print(plot)
    }
    )      

      
  }
)


 
    
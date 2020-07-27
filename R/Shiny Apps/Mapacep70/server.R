
library(shiny)
load("cep702.rda")
library(ggplot2)

shinyServer(
  function(input, output) {
   
    
#     
#     output$newHist <- reactivePlot( function() {
#       plot <- ggplot(cep70, aes(x=te2p05, fill=get(input$xcol) )) + geom_density(alpha=.3) + facet_grid( ddp31 ~ .) + theme_bw()
#       print(plot)
#     }
#     )      

environment<-environment()    
plot <- ggplot(cep702, aes(x=get('<- Izquierda - Derecha ->'), fill=get(input$xcol) ), environment = environment)+ geom_density(alpha=.3) + labs(x="<- Izquierda - Derecha ->", y="densidad") + guides(fill=guide_legend(title=NULL)) + scale_x_continuous(breaks=seq(0, 0, 0)) + theme_bw()
  
    output$newHist <- reactivePlot( function() {
      print(plot)
    }
    )          
      
  }
)


 
    
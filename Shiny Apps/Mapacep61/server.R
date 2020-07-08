
library(shiny)
load("cep61.Rda")
library(ggplot2)

shinyServer(
  function(input, output) {
    
    
    
    
    environment<-environment()    
    
    # plot <- ggplot(cep702, aes(x=get('<- Izquierda - Derecha ->'), fill=get(input$xcol) ), 
    #environment = environment)+ geom_density(alpha=.3) 
    #+ labs(x="<- Izquierda - Derecha ->", y="densidad") + guides(fill=guide_legend(title=NULL)) + 
    #scale_x_continuous(breaks=seq(0, 0, 0)) + theme_bw()
    
    plot <- ggplot(cep611,aes(x=X,y=Y,color=get(input$xcol)), environment = environment) + geom_point(shape=1,position=position_jitter(width=1,height=.5))+ labs(x="<- Izquierda - Derecha ->", y="Afección Política")  + guides(color=guide_legend(title=NULL))+ coord_cartesian(xlim = c(-3, 3))+ coord_cartesian(ylim = c(-3, 3))  + guides(fill=guide_legend(title=NULL)) + theme_bw()
    
    
    output$newHist <- reactivePlot( function() {
      print(plot)
    }
    )          
    
  }
)


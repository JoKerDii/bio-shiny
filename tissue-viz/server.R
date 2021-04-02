require(dplyr)
require(ggplot2)
require(shiny)
# install.packages("shinyWidgets")
require(shinyWidgets)
library(shinyWidgets)
library(shinycssloaders)
# library(highcharter)
library(SummarizedExperiment)
library(stringr)
require(plotly)
library(httr)
library(rdrop2)
library(tidyr)
library(DT)


shinyServer(function(input, output, session) {

  
  output$PCA <- renderPlotly({
    
    if(input$MeasurePicker == "tissue") {
      mynewpc <- read.csv("data/myPC.csv")
      
    } else {
      mynewpc <- read.csv("data/myPC.csv")
    }
    
    fig <- plot_ly(mynewpc, x = ~PC1, y = ~PC2, z = ~PC3, color = ~tissue, 
                   colors = c('#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67','#F79F1F','#A3CB38'))
    fig <- fig %>% add_markers()
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                       yaxis = list(title = 'PC2'),
                                       zaxis = list(title = 'PC3')))
    fig
    
  })
  
  
  
})

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
    
    if(input$MeasurePicker == "Microarray") {
      mynewpc <- read.csv("data/tissue_pc.csv")
      
      fig <- plot_ly(mynewpc, x = ~PC1, y = ~PC2, z = ~PC3, color = ~tissue, 
                     colors = c('#ff0000','#ff4da6','#ffcc00','#47d147','#33ccff','#0066ff','#6600ff'), alpha = 0.7)
      fig <- fig %>% add_markers()
      fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                         yaxis = list(title = 'PC2'),
                                         zaxis = list(title = 'PC3')))
      fig
      
    }  else if (input$MeasurePicker == "RNAseq") {
      mynewpc <- read.csv("data/cuffdiff_pc.csv")
      
      fig <- plot_ly(mynewpc, x = ~PC1, y = ~PC2, z = ~PC3, color = ~group, 
                     colors = c('#ff0000','#ff6666','#cc0000','#00aaff','#80d4ff','#0077b3'))
      fig <- fig %>% add_markers()
      fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                         yaxis = list(title = 'PC2'),
                                         zaxis = list(title = 'PC3')))
      fig
    } else if (input$MeasurePicker == "MeRIPseq") {
      mynewpc <- read.csv("data/MeRIPseq_pc.csv")
      
      fig <- plot_ly(mynewpc, x = ~PC1, y = ~PC2, z = ~PC3, color = ~group, 
                     colors = c('#ff0000','#ff6666','#cc0000','#00aaff','#80d4ff','#0077b3'))
      fig <- fig %>% add_markers()
      fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                         yaxis = list(title = 'PC2'),
                                         zaxis = list(title = 'PC3')))
      fig
    }
      else {
      mynewpc <- read.csv("data/tissue_pc.csv")
      
      fig <- plot_ly(mynewpc, x = ~PC1, y = ~PC2, z = ~PC3, color = ~tissue, 
                     colors = c('#ff0000','#ff4da6','#ffcc00','#47d147','#33ccff','#0066ff','#6600ff'))
      fig <- fig %>% add_markers()
      fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                         yaxis = list(title = 'PC2'),
                                         zaxis = list(title = 'PC3')))
      fig
    }
    
    
    
  })
  
  
  
})

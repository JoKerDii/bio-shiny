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
library(Rtsne)
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
      # mynewpc <- read.csv("data/cuffdiff_pc.csv")
      # mynewpc <- read.csv("data/cuffdiff_multi_pc.csv")
      mynewpc <- read.csv("data/cuffdiff_three_pc.csv")
      
      
      fig <- plot_ly(mynewpc, x = ~PC1, y = ~PC2, z = ~PC3, color = ~group, 
                     colors = c('#ff0000','#ff6666','#cc0000',"#c94f62",
                                '#00aaff','#80d4ff','#0077b3',"#4fc9b6",
                                # "#770e48","#bc4d8a", "#9868be", "#6b4688",
                                "#0e773d","#4dbc7f","#3bcf1a", "#81c146"))
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
  
  output$MDS <- renderPlotly({
    if (input$MeasurePicker == "Microarray"){
      mds_tissue <- read.csv("data/tissue_mds.csv")
      fig <- plot_ly(mds_tissue, x = ~V1, y = ~V2, z = ~V3, color = ~tissue, 
                     colors = c('#ff0000','#ff4da6','#ffcc00','#47d147','#33ccff','#0066ff','#6600ff'), alpha = 0.7)
      fig <- fig %>% add_markers()
      fig <- fig %>% layout(scene = list(xaxis = list(title = 'V1'),
                                         yaxis = list(title = 'V2'),
                                         zaxis = list(title = 'V3')))
      fig
    } else {
      mds_tissue <- read.csv("data/tissue_mds.csv")
      fig <- plot_ly(mds_tissue, x = ~V1, y = ~V2, z = ~V3, color = ~tissue, 
                     colors = c('#ff0000','#ff4da6','#ffcc00','#47d147','#33ccff','#0066ff','#6600ff'), alpha = 0.7)
      fig <- fig %>% add_markers()
      fig <- fig %>% layout(scene = list(xaxis = list(title = 'V1'),
                                         yaxis = list(title = 'V2'),
                                         zaxis = list(title = 'V3')))
      fig
    }
  })
  
  output$tsne <- renderPlotly({
    if (input$MeasurePicker == "Microarray"){
      # data
      allpcs <- read.csv("data/allpcs_tissue.csv")
      tissue <- read.csv("data/tissue.csv") 
      tissue <- tissue[[1]]
      tsne.results <- Rtsne(allpcs, check_duplicates = FALSE, dims=3, theta=0.0, pca=FALSE, verbose=FALSE, max_iter=input$Max_Iteration, perplexity=input$Perplexity)
      colnames(tsne.results$Y) <- c("V1", "V2", "V3")
      mytsnedata <- data.frame(tsne.results$Y)
      mytsnedata1 <- mytsnedata %>% mutate(tissue = tissue)
      head(mytsnedata1)
      # plot
      fig <- plot_ly(mytsnedata1, x = ~V1, y = ~V2, z = ~V3, color = ~tissue, 
                     colors = c('#ff0000','#ff4da6','#ffcc00','#47d147','#33ccff','#0066ff','#6600ff'), alpha = 0.7)
      fig <- fig %>% add_markers()
      fig <- fig %>% layout(scene = list(xaxis = list(title = 'V1'),
                                         yaxis = list(title = 'V2'),
                                         zaxis = list(title = 'V3')))
      fig
    } else {
      # data
      allpcs <- read.csv("data/allpcs_tissue.csv")
      tissue <- read.csv("data/tissue.csv") 
      tsne.results <- Rtsne(allpcs, check_duplicates = FALSE, dims=3, theta=0.0, pca=FALSE, verbose=FALSE, max_iter=input$Max_Iteration, perplexity=input$Perplexity)
      colnames(tsne.results$Y) <- c("V1", "V2", "V3")
      mytsnedata <- data.frame(tsne.results$Y)
      mytsnedata <- mytsnedata %>% mutate(tissue = tissue)
      
      # plot
      fig <- plot_ly(mytsnedata, x = ~V1, y = ~V2, z = ~V3, color = ~tissue, 
                     colors = c('#ff0000','#ff4da6','#ffcc00','#47d147','#33ccff','#0066ff','#6600ff'), alpha = 0.7)
      fig <- fig %>% add_markers()
      fig <- fig %>% layout(scene = list(xaxis = list(title = 'V1'),
                                         yaxis = list(title = 'V2'),
                                         zaxis = list(title = 'V3')))
      fig
    }
  })
  

  
  
  
})

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
    
    if(input$MeasurePicker1 == "Microarray") {
      mynewpc <- read.csv("data/tissue_pc.csv")
      
      fig <- plot_ly(mynewpc, x = ~PC1, y = ~PC2, z = ~PC3, color = ~tissue, 
                     colors = c('#ff0000','#ff4da6','#ffcc00','#47d147','#33ccff','#0066ff','#6600ff'), alpha = 0.7)
      fig <- fig %>% add_markers()
      fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                         yaxis = list(title = 'PC2'),
                                         zaxis = list(title = 'PC3')))
      fig
      
    }  else if (input$MeasurePicker1 == "RNAseq") {
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
    } else if (input$MeasurePicker1 == "MeRIPseq") {
      mynewpc <- read.csv("data/MeRIPseq_pc.csv")
      
      fig <- plot_ly(mynewpc, x = ~PC1, y = ~PC2, z = ~PC3, color = ~group, 
                     colors = c('#ff0000','#ff6666','#cc0000','#00aaff','#80d4ff','#0077b3'))
      fig <- fig %>% add_markers()
      fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                         yaxis = list(title = 'PC2'),
                                         zaxis = list(title = 'PC3')))
      fig
    } else if (input$MeasurePicker1 == "Single-Cell-RNAseq") {
      mynewpc <- read.csv("data/scRNAseq_pc.csv")
      # mynewpc <- mynewpc %>% mutate(label = as.factor(label))
      fig <- plot_ly(mynewpc, x = ~PC1, y = ~PC2, z = ~PC3, color = ~as.factor(label), 
                     colors = c('#FE642E','#7cc652','#2ECCFA','#8181F7','#F7819F'), alpha = 0.7)
      fig <- fig %>% add_markers()
      fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1'),
                                         yaxis = list(title = 'PC2'),
                                         zaxis = list(title = 'PC3')))
      fig
    }
  })
  
  output$MDS <- renderPlotly({
    if (input$MeasurePicker2 == "Microarray"){
      mds_tissue <- read.csv("data/tissue_mds.csv")
      fig <- plot_ly(mds_tissue, x = ~V1, y = ~V2, z = ~V3, color = ~tissue, 
                     colors = c('#ff0000','#ff4da6','#ffcc00','#47d147','#33ccff','#0066ff','#6600ff'), alpha = 0.7)
      fig <- fig %>% add_markers()
      fig <- fig %>% layout(scene = list(xaxis = list(title = 'V1'),
                                         yaxis = list(title = 'V2'),
                                         zaxis = list(title = 'V3')))
      fig
    } else if (input$MeasurePicker2 == "Single-Cell-RNAseq"){
      mds_scRNAseq <- read.csv("data/scRNAseq_mds.csv")
      # mds_scRNAseq <- mds_scRNAseq %>% mutate(label = as.factor(label))
      fig1 <- plot_ly(mds_scRNAseq, x = ~V1, y = ~V2, z = ~V3, color = ~as.factor(label), 
                     colors = c('#FE642E','#7cc652','#2ECCFA','#8181F7','#F7819F'), alpha = 0.7)
      fig1 <- fig1 %>% add_markers()
      fig1 <- fig1 %>% layout(scene = list(xaxis = list(title = 'V1'),
                                         yaxis = list(title = 'V2'),
                                         zaxis = list(title = 'V3')))
      fig1
    } 
  })
  
  output$tsne <- renderPlotly({
    if (input$MeasurePicker3 == "Microarray"){
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
    } else if (input$MeasurePicker3 == "Single-Cell-RNAseq"){
      # data
      allpcs_scRNAseq <- read.csv("data/scRNAseq_allpc.csv")
      
      tsne.results_scRNAseq <- Rtsne(allpcs_scRNAseq[, 1:50], check_duplicates = FALSE, dims=3, theta=0.0, pca=FALSE, verbose=FALSE, max_iter=input$Max_Iteration, perplexity=input$Perplexity)
      colnames(tsne.results_scRNAseq$Y) <- c("V1", "V2", "V3")
      mytsnedata_scRNAseq <- data.frame(tsne.results_scRNAseq$Y)
      mytsnedata_scRNAseq <- mytsnedata_scRNAseq %>% mutate(label = as.factor(allpcs_scRNAseq$label))
      
      # plot
      fig1 <- plot_ly(mytsnedata_scRNAseq, x = ~V1, y = ~V2, z = ~V3, color = ~label, 
                     colors = c('#FE642E','#7cc652','#2ECCFA','#8181F7','#F7819F'), alpha = 0.7)
      fig1 <- fig1 %>% add_markers()
      fig1 <- fig1 %>% layout(scene = list(xaxis = list(title = 'V1'),
                                         yaxis = list(title = 'V2'),
                                         zaxis = list(title = 'V3')))
      fig1
    }
  })
  
})

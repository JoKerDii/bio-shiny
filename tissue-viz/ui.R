
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

shinyUI(fluidPage(theme = "style1.css",
                  titlePanel(HTML("<h3>Clustering - visualizer</h3>"),windowTitle = "Clustering - visualizer"),
                  navbarPage(
                    
                    # One tab for each plot/table.
                    tabsetPanel(
                      
                      type = "tabs",
                      
                      # Bar plot of enriched GO terms
                      tabPanel(
                        
                        "Principal Components Analysis (PCA)",
                        
                        # Sidebar panel for controls.
                        sidebarPanel(
                          pickerInput(
                            "MeasurePicker", "Choose a dataset",
                            
                            choices = c("Microarray","RNAseq", "MeRIPseq"),
                            selected = "Microarray",
                            multiple = F
                          ),
                          
                          tags$p(HTML("PCA is a demensionality reduction technique, which can be used for visualizing high dimensional data. Three coordinates corresponds to three principal components, which capture most of the variance in the data.")),
                          tags$p(HTML("'Microarray' is a tabular data containing expression levels obtained from a set of microarray experiments. There are 22215 genes in 189 samples from 7 human tissues.")),
                          tags$p(HTML("'RNAseq' is a summarized result of differential expression analysis with CuffDiff. It contains FPKM values for 26260 genes in two experimental conditions.")),
                          tags$p(HTML("'MeRIPseq' is a summarized result of differential expression analysis with Ballgown. It contains FPKM values for 87841 genes in two experimental conditions."))
                          
                        ),
                        
                        # Main panel with plot.
                        mainPanel(
                          plotlyOutput("PCA",width = "950px", height = "600px")
                        ),
                        hr(),
                        div(#class = "footer",
                          tags$br(),
                          hr(),
                          p("App created by Di Zhen in April 2021", HTML("&bull;"), 
                            "Find the code on Github:", tags$a(href = "https://github.com/JoKerDii/", "JoKerDii", tags$i(class = 'fa fa-github', style = 'color:#5000a5'), target = '_blank'), style = "font-size: 80%"),
                          p(tags$em("Last updated: April 2021"), style = 'font-size:75%'))
                        
                      ),
                      
                      # Table for Gene specific functions.
                      tabPanel(

                        "What's More",

                      #   # Sidebar panel for controls.
                        sidebarPanel(
                          pickerInput(
                            "MeasurePicker", "Choose a dataset",
                            
                            choices = c("tissue"),
                            selected = "tissue",
                            multiple = F
                          ),
                          tags$p(HTML("More algorithms are coming soon, such as t-SNE and MDS if possible. And maybe use other interesting data such as scRNAseq data.")),
                          tags$p(HTML("jajaja..."))
                        ),
                      #   
                      #   # Main panel with table.
                        mainPanel(
                          # dataTableOutput("genefunctiontable")
                        ),
                        hr(),
                        div(#class = "footer",
                          tags$br(),
                          hr(),
                          p("App created by Di Zhen in Feb 2021", HTML("&bull;"),
                            "Find the code on Github:", tags$a(href = "https://github.com/JoKerDii/", "JoKerDii", tags$i(class = 'fa fa-github', style = 'color:#5000a5'), target = '_blank'), style = "font-size: 80%"),
                          p(tags$em("Last updated: Feb 2021"), style = 'font-size:75%'))
                      #   
                      )
                      
                    )
                    
                  ),
                  
)
)

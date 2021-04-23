
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
                  titlePanel(HTML("<h3>High Dimensional Data - visualizer</h3>"),windowTitle = "High Dimensional Data - visualizer"),
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
                            "MeasurePicker1", tags$h5("Choose a dataset:"),
                            
                            choices = c("Microarray", "Single-Cell-RNAseq","RNAseq", "MeRIPseq"),
                            selected = "Microarray",
                            multiple = F
                          ),
                          
                          tags$p(HTML("PCA is a demensionality reduction technique, which can be used for visualizing high dimensional data. Three coordinates corresponds to three principal components, which capture most of the variance in the data.")),
                          tags$p(HTML("'Microarray' is a tabular data containing expression levels obtained from a set of microarray experiments. There are 22215 genes in 189 samples from 7 human tissues.")),
                          tags$p(HTML("'Single-Cell-RNAseq' comes from single-cell RNA-seq count matrix. There are 511 cells in the rows and 45768 genes in the columns. The expression data were transformed as FPKM.")),
                          tags$p(HTML("'RNAseq' is a summarized result of differential expression analysis with CuffDiff. It contains FPKM values for 26260 genes in three experimental conditions.")),
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
                      
                      
                      tabPanel(

                        "Multi-Dimensional Scaling (MDS)",

                      #   # Sidebar panel for controls.
                        sidebarPanel(
                          pickerInput(
                            "MeasurePicker2", tags$h5("Choose a dataset:"),
                            
                            choices = c("Microarray","Single-Cell-RNAseq"),
                            selected = "Microarray",
                            multiple = F
                          ),
                          tags$p(HTML("MDS is another demensionality reduction technique, which computes distances between all data points and tries to preserve the relative distance between high-dimensional vectors in a low-dimensional space.")),
                          tags$p(HTML("MDS is similar to PCA in that they are both computing eigenvectors and eigenvalues, but the difference is that they are dealing with different matrices, which also results in different computational complexity.")),
                          tags$p(HTML("'Microarray' is a tabular data containing expression levels obtained from a set of microarray experiments. There are 22215 genes in 189 samples from 7 human tissues.")),
                          tags$p(HTML("'Single-Cell-RNAseq' comes from single-cell RNA-seq count matrix. There are 511 cells in the rows and 45768 genes in the columns. The expression data were transformed as FPKM."))
                        ),
                      #   
                      #   # Main panel with table.
                        mainPanel(
                          plotlyOutput("MDS",width = "950px", height = "600px")
                        ),
                        hr(),
                        div(#class = "footer",
                          tags$br(),
                          hr(),
                          p("App created by Di Zhen in April 2021", HTML("&bull;"),
                            "Find the code on Github:", tags$a(href = "https://github.com/JoKerDii/", "JoKerDii", tags$i(class = 'fa fa-github', style = 'color:#5000a5'), target = '_blank'), style = "font-size: 80%"),
                          p(tags$em("Last updated: April 2021"), style = 'font-size:75%'))
                      #   
                      ),
                      
                      tabPanel(
                        
                        "t-distributed Stochastic Neighbor Embedding(t-SNE)",
                        
                        #   # Sidebar panel for controls.
                        sidebarPanel(
                          pickerInput(
                            "MeasurePicker3", tags$h5("Choose a dataset:"),
                            
                            choices = c("Microarray","Single-Cell-RNAseq"),
                            selected = "Microarray",
                            multiple = F
                          ),
                          sliderInput("Perplexity",
                                      tags$h5("Choose the number of neighbors:"),
                                      min = 1,
                                      max = 50,
                                      value = 10,
                                      step = 1
                          ),
                          sliderInput("Max_Iteration",
                                      tags$h5("Choose the maximum number of iterations:"),
                                      min = 10,
                                      max = 2000,
                                      value = 1000,
                                      step = 10
                          ),
                          tags$p(HTML("t-SNE is another demensionality reduction technique for visualizing high-dimensional data.")),
                          tags$p(HTML("t-SNE converts similarities between data points to joint probabilities and tries to minimize the Kullback-Leibler divergence between the joint probabilities of the low-dimensional embedding and the high-dimensional data.")),
                          tags$p(HTML("'Microarray' is a tabular data containing expression levels obtained from a set of microarray experiments. There are 22215 genes in 189 samples from 7 human tissues.")),
                          tags$p(HTML("'Single-Cell-RNAseq' comes from single-cell RNA-seq count matrix. There are 511 cells in the rows and 45768 genes in the columns. The expression data were transformed as FPKM."))
                        ),
                        #   
                        #   # Main panel witht able.
                        mainPanel(
                          plotlyOutput("tsne",width = "950px", height = "600px")
                        ),
                        hr(),
                        div(#class = "footer",
                          tags$br(),
                          hr(),
                          p("App created by Di Zhen in April 2021", HTML("&bull;"),
                            "Find the code on Github:", tags$a(href = "https://github.com/JoKerDii/", "JoKerDii", tags$i(class = 'fa fa-github', style = 'color:#5000a5'), target = '_blank'), style = "font-size: 80%"),
                          p(tags$em("Last updated: April 2021"), style = 'font-size:75%'))
                        #  
                      )
                      # ),
                      # 
                      # # Table for Gene specific functions.
                      # tabPanel(
                      #   
                      #   "What's More",
                      #   
                      #   #   # Sidebar panel for controls.
                      #   sidebarPanel(
                      #     pickerInput(
                      #       "MeasurePicker", tags$h5("Choose a dataset:"),
                      #       
                      #       choices = c("Microarray"),
                      #       selected = "Microarray",
                      #       multiple = F
                      #     ),
                      #     tags$p(HTML("More algorithms are coming soon, such as t-SNE if possible. And maybe use other interesting data such as scRNAseq data.")),
                      #     tags$p(HTML("jajaja..."))
                      #   ),
                      #   #   
                      #   #   # Main panel with table.
                      #   mainPanel(
                      #     # dataTableOutput("genefunctiontable")
                      #   ),
                      #   hr(),
                      #   div(#class = "footer",
                      #     tags$br(),
                      #     hr(),
                      #     p("App created by Di Zhen in April 2021", HTML("&bull;"),
                      #       "Find the code on Github:", tags$a(href = "https://github.com/JoKerDii/", "JoKerDii", tags$i(class = 'fa fa-github', style = 'color:#5000a5'), target = '_blank'), style = "font-size: 80%"),
                      #     p(tags$em("Last updated: April 2021"), style = 'font-size:75%'))
                      # ) # end of the last tabPanel
                    ) # end of the tabsetPanel
        ) # end of the navbarPage
)# end of the fluidPage
) #end of shiny UI

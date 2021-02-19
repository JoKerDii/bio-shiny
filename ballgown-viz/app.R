
# rm(list= ls())
# setwd("/data/zhendi/protocol/ballgown/ballgown-viz")

# initializer <- function(gown_file){
# Load R packages

# library(ballgown)
library(genefilter)
library(ggplot2)
library(tidyr)
library(shiny)
library(shinythemes)
library(dplyr)
library(gplots)


# load("homo_bg.rda")
# # bg <- get(ls()[ls() != gown_file])
# group <- c(rep("experimental", 3), rep("control", 3))
# pData(bg) = data.frame(id=sampleNames(bg), group=group)
# 
# sampleNames_bg <- sampleNames(bg)
# # Analysis with Ballgown
# # Extract FPKM values from the 'bg' object
# fpkm = texpr(bg,meas="FPKM")
# # Transform the FPKM values by adding 1 and convert to a log2 scale
# fpkm = log2(fpkm+1)
# # Load all attributes and gene names
# bg_table = texpr(bg, 'all')
# bg_gene_names = unique(bg_table[, 9:10]) # not required
# # Pull the gene_expression data frame from the ballgown object
# gene_expression = as.data.frame(gexpr(bg))
# # Load the transcript to gene index from the ballgown object
# transcript_gene_table = indexes(bg)$t2g
# #Each row of data represents a transcript. Many of these transcripts represent the same gene. Determine the numbers of transcripts and unique genes
# # Differential expression results
# results_genes = stattest(bg, feature="gene", covariate="group", getFC=TRUE, meas="FPKM") #not required
# results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))
# # save variables
# variables <- c("sampleNames_bg", "fpkm", "bg_table", "gene_expression", "transcript_gene_table", "results_genes")
# ls()
# save(list=variables, file="myvariables.RData")
# # rm(list=ls())
# # ls()
load("myvariables.RData")
sampleNames_bg <- pData_bg[,1]
group <- pData_bg[,2]


# Set some variable values
# colours()
data_colors=c("tomato1","tomato2","tomato3","royalblue1","royalblue2","royalblue3")

# Load helper functions
source('./figure_helper.R', echo=TRUE)
source('./table_helper.R', echo=TRUE)


# shiny
# Define UI for application 
ui <- fluidPage(
    theme = shinytheme('flatly'),

    # App title
    titlePanel(tags$b("Ballgown Visualizer", style = "font-size: 110%, font-family:Helvetica; color:#010151"), windowTitle = "Ballgown Visualizer"),
    hr(),

    # App Description
    p("This is a visualization tool designed specifically for differential expression analysis pipline: HISAT2 - StringTie - Ballbown, aiming to explore the high-throughput sequencing data and the analysis results.", style = "font-size: 100%"),
    hr(),

    # 1. Plot01
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot 1"),
        hr(),
        tags$h5("A boxplot to display summary statistics for the FPKM values for each sample.")
      ),
      # Output of plots
      mainPanel(
        plotOutput("boxPlot01", height = 500)
      )
    ),
    br(),
    hr(),

    # 2. Plot02
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot 2"),
        hr(),
        sliderInput("Gene1",
                    tags$h4("Index of genes:"),
                    min = 1,
                    max = nrow(fpkm),
                    value = 30,
                    step = 1)
      ),
      # Output of plots
      mainPanel(
        plotOutput("boxPlot02", height = 500)
      )
    ),
    br(),
    hr(),


    # 4. Plot04
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot 3"),
        hr(),
        tags$h5("View the distribution of differential expression values as a histogram."),
        tags$h6("Display only those that are significant according to Ballgown.")
      ),
      # Output of plots
      mainPanel(
        plotOutput("diff_exp_hist", height = 500)
      )
    ),
    
    br(),
    hr(),
    # 5. Plot05
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot 4"),
        hr(),
        tags$h5("Display the expression values from experimental and control groups, and mark those that are significantly differentially expressed.")
      ),
      # Output of plots
      mainPanel(
        plotOutput("sig_diff", height = 500)
      )
    ),
    br(),
    hr(),
    
    # plot 06
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot 5"),
        hr(),
        tags$h5("A heatmap vizualizes the expression differences between all samples.")
      ),
      # Output of plots
      mainPanel(
        plotOutput("heatmap", height = 500)
      )
    ),
    br(),
    hr(),
    
    ### can be visualized if gown object is available
    # # 3. Plot03: transcript structure
    # sidebarLayout(
    #   sidebarPanel(
    #     tags$h3("Plot 6"),
    #     hr(),
    #     sliderInput("Gene2",
    #                 tags$h4("Index of genes:"),
    #                 min = 1,
    #                 max = nrow(fpkm),
    #                 value = 30,
    #                 step = 1),
    #     hr(),
    #     selectizeInput(inputId = "Sample",
    #                    label = tags$h4("Select a Sample:"),
    #                    choices = sampleNames_bg,
    #                    # selected = NULL,
    #                    # multiple = T,
    #                    options = list(
    #                      placeholder = 'Please select an option below'
    #                    )),
    #   ),
    #   # Output of plots
    #   mainPanel(
    #     plotOutput("trans_str", height = 500)
    #   )
    # ),
    # br(),
    # hr(),
    
    # makeTable2
    sidebarLayout(
      sidebarPanel(
        tags$h3("Table 1"),
        hr(),
        tags$h5("A table summarized information about differential expressed genes.")
      ),
      mainPanel(
        DT::dataTableOutput("table2")
      )),
    

    br(),
    hr(),
    br(),
    tags$h3("Supplementary"),
    hr(),
    br(),

    # plot #1
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot #1"),
        hr(),
        tags$h5("This plot displays the number of transcripts per gene.")
      ),
      mainPanel(
        plotOutput("transcriptCountPerGene", height = 500)
      )),
    br(),
    hr(),

    # plot #2
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot #2"),
        hr(),
        tags$h5("This plot displays the number of transcripts per gene."),
        tags$h6("If we supplied StringTie with transcript models, the lengths will be those of known transcripts.
    However, if we had used a de novo transcript discovery mode, this step would give us some idea of how well transcripts were being assembled.
    If we had a low coverage library, or other problems, we might get short 'transcripts' that are actually only pieces of real transcripts")
      ),
      mainPanel(
        plotOutput("transcriptSizeLength", height = 500)
      )),
    br(),
    hr(),

    # plot #3
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot #3"),
        hr(),
        tags$h5("View the range of values and general distribution of FPKM values for all samples.")
      ),
      mainPanel(
        plotOutput("GeneFPKMdist", height = 500)
      )),
    br(),
    hr(),
    
    # Plot #4: pair of replicates
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot #4"),
        hr(),
        tags$h5("Plot a pair of replicates to assess reproducibility of technical replicates, or compare control and experimental groups from the same library. "),
        tags$h6("The data are transformed to log2 scale."),

        hr(),
        selectizeInput(inputId = "Sample2",
                       label = tags$h4("Select a Sample:"),
                       choices = sampleNames_bg,
                       selected = c(sampleNames_bg[1], sampleNames_bg[2]),
                       multiple = TRUE, 
                       options = list(
                         maxItems = 2
                       )),
      ),
      # Output of plots
      mainPanel(
        plotOutput("pairOfReplicates", height = 500)
      )
    ),
    br(),
    hr(),
    
    # Plot #5 MDS distance
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot #5"),
        hr(),
        tags$h5("Display the relative differences between libraries."),
        tags$h6("The correlations are converted to 'distance', and multi-dimensional scaled. This step calculates 2-dimensional coordinates to plot points for each library. Libraries with similar expression patterns (highly correlated to each other) should group together.")
      ),
      # Output of plots
      mainPanel(
        plotOutput("MDSdistance", height = 500)
      )),
    br(),
    hr(),
    
    # Table 1
    sidebarLayout(
      sidebarPanel(
        tags$h3("Table #2"),
        hr(),
        tags$h5("Compare the correlation between all replicates.")
      ),
      mainPanel(
        DT::dataTableOutput("table1")
      )),
    br(),
    hr(),
    
    
    # Footer
    tags$br(),
    hr(),
    p("App created by Di Zhen in Feb 2021", HTML("&bull;"), "Find the code on Github:", tags$a(href = "https://github.com/JoKerDii/m6a-seq-analysis-visualizer/", "JoKerDii", tags$i(class = 'fa fa-github', style = 'color:#5000a5'), target = '_blank'), style = "font-size: 80%"),
    p(tags$em("Last updated: Feb 2021"), style = 'font-size:75%')

  )

# Define server logic
server <- function(input, output) {
  # plot01
  output$boxPlot01 <- renderPlot({
    boxplot01(fpkm, pData_bg)
  })
  # plot02
  output$boxPlot02 <- renderPlot({
    boxplot02(fpkm, input$Gene1, pData_bg)
  })
  # plot03 - transcript structure
  output$trans_str <- renderPlot({
    trans_structure(input$Gene2, input$Sample, geneIDs)
  })
  # plot04 - differential expression
  output$diff_exp_hist <- renderPlot({
    diff_exp_hist(results_genes)
  })
  # plot05 - differential expression + significant
  output$sig_diff <- renderPlot({
    sig_diff(gene_expression,results_genes,pData_bg)
  })
  # plot06 - heatmap
  output$heatmap <- renderPlot({
    heatmap(gene_expression,results_genes)
  })
  

  # Table 2 correlation 'distance'
  output$table2 <- DT::renderDataTable(DT::datatable({
    makeTable2(results_genes)
    
  }))
  

  ##### Supplementary 
  # Plot #1 - the number of transcripts per gene.
  output$transcriptCountPerGene <- renderPlot({
    transcript_per_gene(transcript_gene_table)
  })

  # Plot #2 - the distribution of transcript lengths.
  output$transcriptSizeLength <- renderPlot({
    transcript_size_dist(bg_table)
  })

  # Plot #3 - the distribution of FPKM value for gene expression.
  output$GeneFPKMdist <- renderPlot({
    boxplot01(gene_expression, pData_bg,transcript = FALSE)
  })
  
  # Plot #4 - pair of replicates
  output$pairOfReplicates <- renderPlot({
    pairOfReplicates(gene_expression, input$Sample2[1], input$Sample2[2])
  })
  
  # Plot #5 - MDS distance
  output$MDSdistance <- renderPlot({
    MDSdistance(gene_expression,pData_bg)
  })
  
  # Table 1 correlation 'distance'
  output$table1 <- DT::renderDataTable(DT::datatable({
    makeTable1(gene_expression)
  }))


}

# Run the application
shinyApp(ui = ui, server = server)


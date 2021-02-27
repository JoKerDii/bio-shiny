shinyUI(fluidPage(theme = "style1.css",
                  titlePanel(HTML("<h3>DAVID Gene Ontology - visualizer</h3>"),windowTitle = "DAVID Gene Ontology - visualizer"),
                  navbarPage(

    # One tab for each plot/table.
    tabsetPanel(

      type = "tabs",

      # Bar plot of enriched GO terms
      tabPanel(

        "Enriched GO Terms",

        # Sidebar panel for controls.
        sidebarPanel(
          pickerInput(
            "MeasurePicker", "Choose one or more cagetories:",
            
            choices = c("UP_KEYWORDS", "GOTERM_MF_DIRECT", "GOTERM_CC_DIRECT", "UP_SEQ_FEATURE",
                        "GOTERM_BP_DIRECT", "INTERPRO", "SMART", "KEGG_PATHWAY", "BIOCARTA", 
                        "PIR_SUPERFAMILY", "OMIM_DISEASE", "BBID", "COG_ONTOLOGY"),
            selected = "GOTERM_BP_DIRECT",
            multiple = F
          ),
          
          tags$p(HTML("<b>Hover</b> to see the p-value and the ratio a term takes.")),
          tags$p(HTML("Use the <b>category</b> picker to choose the categories of GO terms of interest.")),
          tags$p(HTML("The order of terms is sorted by p-values in ascending order. The bar plot displays the number of genes."))
        ),

        # Main panel with plot.
        mainPanel(
          highchartOutput("GOtermbarplot",
                          width = "800px", height = "600px")
        ),
        hr(),
        div(#class = "footer",
          tags$br(),
          hr(),
          p("App created by Di Zhen in Feb 2021", HTML("&bull;"), 
            "Find the code on Github:", tags$a(href = "https://github.com/JoKerDii/m6a-seq-analysis-visualizer", "JoKerDii", tags$i(class = 'fa fa-github', style = 'color:#5000a5'), target = '_blank'), style = "font-size: 80%"),
          p(tags$em("Last updated: Feb 2021"), style = 'font-size:75%'))

      ),
      
      # Table for Gene specific functions.
      tabPanel(

        "Gene Specific Functions",

        # Sidebar panel for controls.
        sidebarPanel(
          pickerInput(
            "SetThemePicker", "Choose one or more categories:",
            choices = c(), multiple = T
          ),
          pickerInput(
            "SetGenePicker", "Choose one or more enriched genes:",
            choices = c(), multiple = T
          )
        ),

        # Main panel with table.
        mainPanel(
          dataTableOutput("genefunctiontable")
        ),
        hr(),
        div(#class = "footer",
            tags$br(),
            hr(),
            p("App created by Di Zhen in Feb 2021", HTML("&bull;"), 
              "Find the code on Github:", tags$a(href = "https://github.com/JoKerDii/m6a-seq-analysis-visualizer", "JoKerDii", tags$i(class = 'fa fa-github', style = 'color:#5000a5'), target = '_blank'), style = "font-size: 80%"),
            p(tags$em("Last updated: Feb 2021"), style = 'font-size:75%'))

      )

    )

),

)
)

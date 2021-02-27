library(dplyr)
library(ggplot2)
library(shiny)
library(shinyWidgets)
library(shinycssloaders)
library(highcharter)
library(stringr)
library(httr)
library(rdrop2)
library(tidyr)
library(DT)
shinyServer(function(input, output, session) {
  
  source("observer.R", local = T)
  source("load_data.R", local = T)
  # source("global.R", local = T)
  
  

  output$GOtermbarplot = renderHighchart({
    
    # Data manipulation
    selectPvalue <- function(data){
      pvals <- seq(from = 1, to = 6000, by = 10)/100000
      for (i in pvals[order(pvals,decreasing = TRUE)]){
        nRows <- data %>%
          filter( PValue < i) %>%
          nrow()
        if (nRows < 20){
          return(i)
          break
        }
      }
      return(i)
    }
  
    temp.chart <- chart %>% select(c("Category","Count", "Term", "%", "PValue")) %>%
      rename("Ratio"=`%`) %>%
      mutate(Term = as.factor(gsub("^.*?~", "",Term)), 
             Ratio = Ratio / 100) %>%
      drop_na() %>%
      as.data.frame() 
    
    if(input$MeasurePicker == "UP_KEYWORDS") {
      temp.chart = temp.chart %>%
        filter(Category == "UP_KEYWORDS")
      p <- selectPvalue(temp.chart)
      temp.chart = temp.chart %>%
        filter(PValue < p)
    } else if(input$MeasurePicker == "GOTERM_MF_DIRECT") {
      temp.chart = temp.chart %>%
        filter(Category == "GOTERM_MF_DIRECT")
      p <- selectPvalue(temp.chart)
      temp.chart = temp.chart %>%
        filter(PValue < p)
    } else if(input$MeasurePicker == "GOTERM_CC_DIRECT") {
      temp.chart = temp.chart %>%
        filter(Category == "GOTERM_CC_DIRECT")
      p <- selectPvalue(temp.chart)
      temp.chart = temp.chart %>%
        filter(PValue < p)
    } else if(input$MeasurePicker == "UP_SEQ_FEATURE") {
      temp.chart = temp.chart %>%
        filter(Category == "UP_SEQ_FEATURE")
      p <- selectPvalue(temp.chart)
      temp.chart = temp.chart %>%
        filter(PValue < p)
    } else if(input$MeasurePicker == "GOTERM_BP_DIRECT") {
      temp.chart = temp.chart %>%
        filter(Category == "GOTERM_BP_DIRECT")
      p <- selectPvalue(temp.chart)
      temp.chart = temp.chart %>%
        filter(PValue < p)
    } else if(input$MeasurePicker == "INTERPRO") {
      temp.chart = temp.chart %>%
        filter(Category == "INTERPRO")
      p <- selectPvalue(temp.chart)
      temp.chart = temp.chart %>%
        filter(PValue < p)
    } else if(input$MeasurePicker == "SMART") {
      temp.chart = temp.chart %>%
        filter(Category == "SMART")
      p <- selectPvalue(temp.chart)
      temp.chart = temp.chart %>%
        filter(PValue < p)
    } else if(input$MeasurePicker == "KEGG_PATHWAY") {
      temp.chart = temp.chart %>%
        filter(Category == "KEGG_PATHWAY")
      p <- selectPvalue(temp.chart)
      temp.chart = temp.chart %>%
        filter(PValue < p)
    } else if(input$MeasurePicker == "BIOCARTA") {
      temp.chart = temp.chart %>%
        filter(Category == "BIOCARTA")
      p <- selectPvalue(temp.chart)
      temp.chart = temp.chart %>%
        filter(PValue < p)
    } else if(input$MeasurePicker == "PIR_SUPERFAMILY") {
      temp.chart = temp.chart %>%
        filter(Category == "PIR_SUPERFAMILY")
      p <- selectPvalue(temp.chart)
      temp.chart = temp.chart %>%
        filter(PValue < p)
    } else if(input$MeasurePicker == "OMIM_DISEASE") {
      temp.chart = temp.chart %>%
        filter(Category == "OMIM_DISEASE")
      p <- selectPvalue(temp.chart)
      temp.chart = temp.chart %>%
        filter(PValue < p)
    } else if(input$MeasurePicker == "BBID") {
      temp.chart = temp.chart %>%
        filter(Category == "BBID")
      p <- selectPvalue(temp.chart)
      temp.chart = temp.chart %>%
        filter(PValue < p)
    } else if(input$MeasurePicker == "COG_ONTOLOGY") {
      temp.chart = temp.chart %>%
        filter(Category == "COG_ONTOLOGY")
      p <- selectPvalue(temp.chart)
      temp.chart = temp.chart %>%
        filter(PValue < p)
    }
    
    temp.chart = temp.chart %>%
      mutate(measure = Ratio) %>%
      arrange(PValue) %>%
      mutate(num.parts.col = log(Count) / max(log(Count)))
    
    temp.chart$Term = factor(temp.chart$Term,levels = temp.chart$Term)
    
    
    # Set the format for the tooltip
    point.format = paste(
      " (P-value: {point.num_pieces})</span><br/><span>",
      'Ratio: ',
      ":\u00A0{point.y}</span>",
      sep = ""
    )
    
    # Make the plot.
    hc = highchart() %>%
      hc_chart(type = "bar") %>%
      hc_xAxis(categories = temp.chart$Term) %>%
      hc_add_series(pointPadding = 0,
                    data = temp.chart %>%
                      mutate(y = measure,
                             num_pieces = PValue),
                    colorByPoint = T,
                    colors = rgb(colorRamp(c("white", "black"))(temp.chart$num.parts.col),
                                 maxColorValue = 255),
                    borderColor = "#000000") %>%
      hc_tooltip(headerFormat = "<span><b>{point.key}</b>",
                 pointFormat = point.format,
                 valueDecimals = 8) %>%
      hc_legend(enabled = F)
    hc
  })
  

  output$genefunctiontable = renderDataTable({
    
    selected.themes = input$SetThemePicker # groups
    selected.ethnicities = input$SetGenePicker # genes
    
    part.table = function(data.df, column.names, columns.to.hide) {
      temp.df = data.df
      # Create the DataTable.
      datatable(temp.df,
                options = list(pageLength = 10,
                               columnDefs = list(list(targets = columns.to.hide,
                                                      visible = F))),
                rownames = F,
                colnames = column.names)
    }
    
    # Get the dataset to display
    temp.mybindeddf = mybindeddf
    
    if(length(selected.themes) > 0) {
      temp.mybindeddf = temp.mybindeddf %>%
        # filter(group %in% gsub(" \\([0-9]+\\)$", "", selected.themes))
        filter(Group %in% selected.themes)
    }
    if(length(selected.ethnicities) > 0) {
      temp.mybindeddf = temp.mybindeddf %>%
        filter(Gene_Name %in% selected.ethnicities)
    }
    
    part.table(temp.mybindeddf,
               columns.to.hide = c(),
               column.names = c("ID", "Gene_Name", "Terms", "Group"))
  })

  
})

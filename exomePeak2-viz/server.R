library(shiny)
library(ggvis)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
# setwd("/home/zhen.di/exomePeak2-viz")
peaks <- read_csv("data/DiffMod.csv")
source("ggcircos_helpers.R")



chroms <- c(1:22, "X", "Y")
lengths <- c(
  249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
  159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
  115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
  59128983, 63025520, 48129895, 51304566, 155270560, 59373566
)

radians_f <- create_radians(chroms[1:23], lengths[1:23])
radians_m <- create_radians(chroms, lengths)


server <- function(input, output, session) {
  re_values <- reactiveValues(
    scaling_factors = rep(1, 24), previous_radians = radians_f,
    chrom_clicked = "1", chroms_selected = NULL, sig = 0.4
  )

  radians <- reactive({
    
    rads <- create_radians(chroms, lengths * re_values$scaling_factors)
    isolate(mid <- mean(rads[names(rads) == re_values$chrom_clicked]))
    isolate(prev_mid <- mean(re_values$previous_radians[names(rads) == re_values$chrom_clicked]))
    offset <- mid - prev_mid
    rads - offset
  })


  track_radians <- reactive({
    create_track_radians(radians(), points_per_track = rep(20, length(radians())))
  })


  output$top_clinvar_genes <- renderTable(
    {

      top_peaks <- peaks %>%
        drop_na() %>%
        filter(padj < 0.05) %>%
        select(c(colnames(peaks)[c(1,2,3,4,6,22)])) %>%
        arrange(padj) %>%
        head(n = 10) %>%
        select(c(colnames(peaks)[c(1,2,3,4,6)]))
      setNames(top_peaks, c("Chr", "Start", "End", "Peak ID", "Strand"))
    },
    include.rownames = FALSE
  )

  seq_df <- reactive({
    
    scale <- re_values$scaling_factors
    create_seq_df(radians(), scale = scale)
  })

  
  peaks_filt <- reactive({#################################################################################################### peaks_filt
    peaks %>%
      drop_na() %>%
      filter(padj <= input$float) %>% # re_values$sig
      mutate(pos = (chromStart + chromEnd) / 2) %>%
      mutate(chr = as.factor(str_remove(chr, "chr"))) %>%
      mutate(copy_number = as.integer(case_when(strand == "+" ~ 1,
                                                strand != "+"  ~ -1)))
  })
  
  peak_pval_plot_data <- reactive({ ######################################################################################## peak_pval_plot_data

    peaks_filt <- peaks_filt()
    
    peaks_filt$is_sig <- ifelse(peaks_filt$padj < 0.05, "Yes", "No")
    points <- fit_points(peaks_filt$chr, peaks_filt$pos, peaks_filt$padj,
                         0.8, 0.6, seq_df(),
                         metadata = peaks_filt[, c("name", "strand", "blockCount", "chromStart", "chromEnd", "is_sig", "padj","copy_number")],
                         min_value = -0.1, max_value = max(peaks_filt$padj) + 0.1
    )
    
    points$id <- paste0("padj", 1:nrow(points))
    points
    
  })
  
  
  # Histogram
  peaks_plot_data <- reactive({ ################################################################################################## peaks_plot_data
    peaks_filt <- peaks_filt()
    peaks_plot_data <- fit_to_seq(peaks_filt$chr, peaks_filt$pos, seq_df(),
                                  metadata = peaks_filt[, c("name", "strand", "copy_number", "blockCount", "chromStart", "chromEnd", "padj")]
    )
    
    peak_inner <- rep(0.8 + 2/ 90, nrow(peaks_plot_data))
    peak_outer <- 0.8 + (as.numeric(peaks_plot_data$copy_number)+2) / 90
    
    copy_neutral <- peak_inner == peak_outer
    peak_inner[copy_neutral] <- 0.8 + 1.5 / 90
    peak_outer[copy_neutral] <- 0.8 + 2.5 / 90
    
    peaks_plot_data$copy_number <- paste0("peak", peaks_plot_data$copy_number)
    
    peaks_plot_data <- data.frame(rbind(peaks_plot_data, peaks_plot_data),
                                  r = c(peak_inner, peak_outer)
    )
    
    peaks_plot_data$id <- paste0("peak", 1:(nrow(peaks_plot_data) / 2)) # useless, same as name

    
    peaks_plot_data %>% filter(seq %in% re_values$chroms_selected) # when chr is clicked, show up
  })
  
  # Lines
  peaks_line_data <- reactive({ ################################################################################################ peaks_line_data
    
    peaks_filt <- peaks_filt() %>%
      group_by(chr) %>%
      arrange(pos)
    
    peaks_line_data <- fit_points(peaks_filt$chr, peaks_filt$pos, peaks_filt$copy_number, 0.9, 0.8, seq_df(), max_value = 6, min_value = -2)
    
    peaks_line_data$opac <- ifelse(peaks_line_data$seq %in% re_values$chroms_selected, 0, 0.7) # when chr is clicked, disappear
    
    peaks_line_data
  })

  
  ### TEXT  
  text_df <- reactive({
    seq_df() %>% mutate(theta = (seq_start + seq_end) / 2, r = 1.05)
  })

  tooltip_fun <- function(data, session) {
    if ("is_sig" %in% names(data)) {
      
      tt_data <- peak_pval_plot_data()

      row <- tt_data[tt_data$id == data$id, ]

      # record_colour <- paste0("rgb(", sample(170:255, 1), ",", sample(150:255, 1), ",", sample(170:255, 1), ")")

      paste0( # red pts, # p value, # add peak name, number of block count info
        "Start: ", unique(row$chromStart), "<br>",
        "End: ", unique(row$chromEnd), "<br>",
        "Strand: ",row$strand, "<br>",
        "Adjusted p-value: ",round(row$padj,8), "<br>"
      )
    } else if ("copy_number" %in% names(data)) {
      
      
      tt_data <- peaks_plot_data()

      rows <- tt_data[tt_data$theta > data$theta - 0.0001 & tt_data$theta < data$theta + 0.0001, ]

      paste0(
        "Start: ", unique(rows$chromStart), "<br>",
        "End: ", unique(rows$chromEnd), "<br>"
      )
    } 
  }


  tooltip_click_fun <- function(data) {
    str(data)
  }

  click_handle <- function(data, location, session) {
    if (is.null(data)) {
      return(NULL)
    }

    isolate(re_values$chrom_clicked <- data$group)

    isolate(re_values$previous_radians <- radians())

    isolate(re_values$scaling_factors[which(chroms == data$group)] <- ifelse(re_values$scaling_factors[which(chroms == data$group)] == 1,
      input$integer, 1
    ))

    isolate(re_values$chroms_selected <- chroms[which(re_values$scaling_factors > 1)])

    print(data)
  }

  fill_domain <- c(
    c(1:22, "X", "Y"), # chromosomes
    "Yes", "No" # significance
  ) 

  fill_range <- c(
    # chromosome colours from Circos
    "#d50000", "#ff1744", "#dd2c00", "#ff6d00", "#ff9100", "#ffab00", "#ffd600", "#ffea00", "#aeea00",
    "#64dd17", "#76ff03", "#00c853", "#00bfa5", "#1de9b6", "#00b8d4", "#0091ea", "#00b0ff", "#2962ff",
    "#304ffe", "#3d5afe", "#6200ea", "#aa00ff", "#d500f9", "#c51162",

    # colours for padj
    "#c41659", "#83aee6"
  )


  stroke_domain <- c(
    c(1:22, "X", "Y"), # chromosomes
    paste0("peak", c(-1,1)), # copy numbers
    "Yes", "No"
  ) # interchromosomal or not

  stroke_range <- c(
    # chromosome colours from Circos
    "#d50000", "#ff1744", "#dd2c00", "#ff6d00", "#ff9100", "#ffab00", "#ffd600", "#ffea00", "#aeea00",
    "#64dd17", "#76ff03", "#00c853", "#00bfa5", "#1de9b6", "#00b8d4", "#0091ea", "#00b0ff", "#2962ff",
    "#304ffe", "#3d5afe", "#6200ea", "#aa00ff", "#d500f9", "#c51162",

    # colours for copy numbers
    "#2f8bfa", rep("#C11F29", 6),"#0C2F60", "green", 

    # colours for links
    "blue", "grey"
  )


  add_tooltip <- function(vis, html, on = c("hover", "click")) {
    on <- match.arg(on)

    show_tooltip2 <- function(data, location, session, ...) {
      if (is.null(data)) {
        hide_tooltip(session)
        return()
      }

      html <- html(data)
      if (is.null(html)) {
        hide_tooltip(session)
      } else {
        show_tooltip(session, location$x + 5, location$y + 5, html)
      }
    }
    hide_tooltip2 <- function(session) {
      hide_tooltip(session)
    }

    switch(on,
      click = handle_click(vis, show_tooltip2),
      hover = handle_hover(vis, show_tooltip2)
    )
  }

  ggvis() %>%
    ## Add outside colored chromosome blocks
    add_track(track_radians, 1, 0.9, fill = ~group, stroke = ~group, fillOpacity := 0.7, fillOpacity.hover := 1) %>%
    ## Line plots
    layer_paths( 
      data = peaks_line_data %>% group_by(seq), ~ sin(theta) * r, ~ cos(theta) * r, interpolate := "basis",
      strokeWidth := 0.8, opacity := ~opac
    ) %>%
    ## histogram
    layer_paths(
      data = peaks_plot_data %>% group_by(theta), ~ sin(theta) * r, ~ cos(theta) * r, stroke = ~copy_number,
      strokeWidth := 2, strokeWidth.hover := 3, strokeOpacity:=0.4
    ) %>%
    ## Background lines for red points
    add_circles(track_radians, r = 0.8 + 2 / 90, opacity := 0.2) %>%
    add_track(track_radians, 0.8, 0.6, strokeOpacity := 0.5, stroke := "black", strokeWidth := 0.5) %>%
    add_circles(track_radians, seq(0.6, 0.8, length.out = 7)[-c(1, 7)], strokeOpacity := 0.3, strokeWidth := 0.5) %>%
    ## Red points from 0.6 - 0.8
    layer_points(
      data = peak_pval_plot_data, ~ sin(theta) * r, ~ cos(theta) * r, fill = ~is_sig, 
      key := ~id, size = ~is_sig, size.hover := 30, strokeOpacity = ~is_sig, strokeOpacity.hover := 1, fillOpacity:=0.3,
      stroke := "grey", strokeWidth := 0.5
    ) %>%

    ## Links
    # layer_paths(
    #   data = struct_plot_data, ~ sin(theta) * r, ~ cos(theta) * r, stroke = ~name_from,
    #   strokeWidth := ~opac, strokeWidth.hover := 2, strokeOpacity := 1, interpolate := "basis"
    # ) %>%
    ## Chromosome names
    layer_text(
      data = text_df, ~ sin(theta) * r, ~ cos(theta) * r, text := ~seq, align := "center", baseline := "middle",
      angle := ~ 180 * (theta - pi * (cos(theta) < 0)) / pi
    ) %>%
    ## Basic circle functions
    add_tooltip(tooltip_fun, "hover") %>%
    add_tooltip(tooltip_click_fun, "click") %>%
    handle_click(click_handle) %>%
        scale_numeric("x", domain = c(-1, 1), nice = FALSE, clamp = TRUE) %>%
        scale_numeric("y", domain = c(-1, 1), nice = FALSE, clamp = TRUE) %>%
    ## Basic circle structure
    scale_ordinal("fill", domain = fill_domain, range = fill_range) %>%
    scale_ordinal("stroke", domain = stroke_domain, range = stroke_range) %>%
    scale_nominal("shape", domain = c("Yes", "No"), range = c("cross", "circle")) %>%
    scale_nominal("size", domain = c("Yes", "No"), range = c(15, 10)) %>%
    scale_nominal("opacity", domain = c("Yes", "No"), range = c(0, 0)) %>%
    # scale_nominal("fillOpacity", domain = c("Yes", "No"), range = c(0.5, 0.5)) %>%
    hide_axis("x") %>%
    hide_axis("y") %>%
    hide_legend(c("fill", "stroke", "shape", "size")) %>%
    set_options(hover_duration = 0, width = 775, height = 775, keep_aspect = TRUE, duration = 1000) %>%
    bind_shiny("plot")
}

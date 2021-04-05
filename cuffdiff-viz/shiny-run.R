# This script is used to run the application defined in app.R in the background
setwd("/home/zhen.di/cw1/DGE/cuffdiff-viz")
options(shiny.autoreload = TRUE)
shiny::runApp()

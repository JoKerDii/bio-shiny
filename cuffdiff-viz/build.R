install.packages('rsconnect')
library(rsconnect)
setwd("/home/zhen.di/cw1/DGE/cuffdiff-viz/")
options(repos = BiocManager::repositories())
options("repos")
rsconnect::deployApp(account='jokerdii')



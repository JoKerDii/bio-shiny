library(dplyr)
library(SummarizedExperiment)
setwd("D:/git/m6a-seq-analysis-visualizer/tissue-viz")
getwd()
tissue <- read.csv("data/tissue.csv") # 189
e <- read.csv("data/expression.csv", row.names = "probeID")
SE <- SummarizedExperiment(assays = I(e), colData = DataFrame(tissue))
tissue <- tissue[[1]]
mydata <- data.frame(t(assay(SE)))
mydata <- mydata %>% mutate(tissue = tissue)
# mydata$tissue 
pc <- prcomp(mydata[,1:22215], center = TRUE,scale. = TRUE)

myPC <- data.frame(pc$x[,1:3])
myPC <- myPC %>% mutate(tissue = tissue)
write.csv(myPC, "data/myPC.csv")
mynewpc <- read.csv("data/myPC.csv")

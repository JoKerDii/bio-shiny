library(dplyr)
library(SummarizedExperiment)
library(data.table)
library(caret)

setwd("D:/git/bio-shiny/tissue-viz")
getwd()

### tissue data
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
write.csv(myPC, "data/tissue_pc.csv")
mypc_tissue <- read.csv("data/tissue_pc.csv")


### Cuffdiff results from RNAseq data

cuffdiff <- read.csv("data/cuffdiff.csv") 
cuffdiff <- tidyr::spread(cuffdiff, rep_name, fpkm) # %>% as.data.frame()
rownames(cuffdiff) <- cuffdiff$gene_id
sample_names <- colnames(cuffdiff)[2:7]
cuffdiff_t <- transpose(cuffdiff[,2:7])[1:6,] %>% as.data.frame()
rownames(cuffdiff_t) <- sample_names

NZV<- nearZeroVar(cuffdiff_t, saveMetrics = TRUE)
cuffdiff_t_filter <- cuffdiff_t[, which(NZV[, "nzv"] == FALSE)]
pc <- prcomp(cuffdiff_t_filter, center = TRUE, scale. = TRUE)
myPC <- data.frame(pc$x[,1:3])
myPC <- myPC %>% mutate(group = sample_names)
write.csv(myPC, "data/cuffdiff_pc.csv")
mypc_cuffdiff <- read.csv("data/cuffdiff_pc.csv")


### ballgown results from MeRIPseq data

MeRIPseq <- read.csv("data/MeRIPseq.csv") # 189
MeRIPseq_t <- transpose(MeRIPseq)
rownames(MeRIPseq_t) <- c("EX_0", "EX_2","EX_3","CTR_0", "CTR_2","CTR_3")
# for (i in 1:length(colnames(MeRIPseq))) {
#   print(strsplit(colnames(MeRIPseq), "\\.", perl=TRUE)[[i]][2])
# }
NZV<- nearZeroVar(MeRIPseq_t, saveMetrics = TRUE)
MeRIPseq_t_filter <- MeRIPseq_t[, which(NZV[, "nzv"] == FALSE)]
pc <- prcomp(MeRIPseq_t_filter, center = TRUE, scale. = TRUE)
myPC <- data.frame(pc$x[,1:3])
myPC <- myPC %>% mutate(group = rownames(MeRIPseq_t))
write.csv(myPC, "data/MeRIPseq_pc.csv")
mypc_MeRIPseq <- read.csv("data/MeRIPseq_pc.csv")


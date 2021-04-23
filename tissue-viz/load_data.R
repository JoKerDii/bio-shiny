library(dplyr)
library(SummarizedExperiment)
library(data.table)
library(caret)
library(Rtsne)

setwd("D:/git/bio-shiny/tissue-viz")
getwd()

# ### tissue data
# tissue <- read.csv("data/tissue.csv") # 189
# e <- read.csv("data/expression.csv", row.names = "probeID")
# SE <- SummarizedExperiment(assays = I(e), colData = DataFrame(tissue))
# tissue <- tissue[[1]]
# mydata <- data.frame(t(assay(SE)))
# mydata <- mydata %>% mutate(tissue = tissue)
# # mydata$tissue
# # NZV<- nearZeroVar(mydata, saveMetrics = TRUE)
# # mydata_filter <- mydata[, which(NZV[, "nzv"] == FALSE)]
# 
# pc <- prcomp(mydata[,1:22215], center = T, scale. = T)
# myPC <- data.frame(pc$x[,1:3])
# myPC <- myPC %>% mutate(tissue = tissue)
# # write.csv(pc$x[,1:50], "data/allpcs_tissue.csv")
# mypc_tissue <- read.csv("data/tissue_pc.csv")
# allpcs <- read.csv("data/allpcs_tissue.csv")
# 
# ### tsne
# # pc <- prcomp(mydata[,1:22215], center = F, scale. = F)
# allpcs <- read.csv("data/allpcs_tissue.csv")
# tissue <- read.csv("data/tissue.csv")
# tsne.results <- Rtsne(allpcs, check_duplicates = FALSE, dims=3, theta=0.0, pca=FALSE, verbose=FALSE, max_iter=1000, perplexity=10)
# colnames(tsne.results$Y) <- c("V1", "V2", "V3")
# mytsnedata <- data.frame(tsne.results$Y)
# mytsnedata <- mytsnedata %>% mutate(tissue = tissue)
# # write.csv(mytsnedata, "data/tsne_three.csv", row.names = FALSE)
# # mypc_tsne <- read.csv("data/tsne_three.csv")
# 
# 
# ### Cuffdiff results from RNAseq data
# 
# cuffdiff <- read.csv("data/cuffdiff.csv")
# cuffdiff <- read.csv("data/cuffdiff_three.csv")
# cuffdiff <- tidyr::spread(cuffdiff, rep_name, fpkm) # %>% as.data.frame()
# # cuffdiff <- cuffdiff[, c(1:5,10:17)] # three
# rownames(cuffdiff) <- cuffdiff$gene_id
# sample_names <- colnames(cuffdiff)[2:length(colnames(cuffdiff))]
# cuffdiff_t <- transpose(cuffdiff[,2:length(colnames(cuffdiff))])[1:length(colnames(cuffdiff))-1,] %>% as.data.frame()
# rownames(cuffdiff_t) <- sample_names
# 
# NZV<- nearZeroVar(cuffdiff_t, saveMetrics = TRUE)
# cuffdiff_t_filter <- cuffdiff_t[, which(NZV[, "nzv"] == FALSE)]
# pc <- prcomp(cuffdiff_t_filter, center = TRUE, scale. = TRUE)
# myPC <- data.frame(pc$x[,1:3])
# myPC <- myPC %>% mutate(group = sample_names)
# # write.csv(myPC, "data/cuffdiff_three_pc.csv")
# mypc_cuffdiff <- read.csv("data/cuffdiff_three_pc.csv")
# 
# 
# ### ballgown results from MeRIPseq data
# 
# MeRIPseq <- read.csv("data/MeRIPseq.csv") # 189
# MeRIPseq_t <- transpose(MeRIPseq)
# rownames(MeRIPseq_t) <- c("EX_0", "EX_2","EX_3","CTR_0", "CTR_2","CTR_3")
# # for (i in 1:length(colnames(MeRIPseq))) {
# #   print(strsplit(colnames(MeRIPseq), "\\.", perl=TRUE)[[i]][2])
# # }
# NZV<- nearZeroVar(MeRIPseq_t, saveMetrics = TRUE)
# MeRIPseq_t_filter <- MeRIPseq_t[, which(NZV[, "nzv"] == FALSE)]
# pc <- prcomp(MeRIPseq_t_filter, center = TRUE, scale. = TRUE)
# myPC <- data.frame(pc$x[,1:3])
# myPC <- myPC %>% mutate(group = rownames(MeRIPseq_t))
# # write.csv(myPC, "data/MeRIPseq_pc.csv")
# mypc_MeRIPseq <- read.csv("data/MeRIPseq_pc.csv")
# 
# ###############################################################
# # MDS for tissue
# tissue <- read.csv("data/tissue.csv") # 189
# e <- read.csv("data/expression.csv", row.names = "probeID")
# SE <- SummarizedExperiment(assays = I(e), colData = DataFrame(tissue))
# tissue <- tissue[[1]]
# mydata <- data.frame(t(assay(SE)))
# 
# mds_mydata <- mydata %>%
#   dist() %>%
#   cmdscale(k=3) %>%
#   as.data.frame()
# 
# mds_mydata <- mds_mydata %>% mutate(tissue = tissue)
# # write.csv(mds_mydata, "data/tissue_mds.csv")
# mds_tissue <- read.csv("data/tissue_mds.csv")
# 
# 
# ####################################################################
# 
# singlecellRNAseq<-read.csv("data/single-cell-rnaseq.csv")
# singlecellRNAseq <- singlecellRNAseq[,-(which(colSums(singlecellRNAseq) == 0))]
# singlecellRNAseq <- singlecellRNAseq %>% mutate(label=as.factor(label))
# # singlecellRNAseq %>% head()
# # singlecellRNAseq[1:10, 1:10]
# 
# ### PCA
# 
# # ncol(singlecellRNAseq)
# pc <- prcomp(singlecellRNAseq[,1:ncol(singlecellRNAseq)-1], center = T, scale. = T)
# myPC <- data.frame(pc$x[,1:3])
# myPC <- myPC %>% mutate(label = singlecellRNAseq$label)
# # write.csv(myPC, "data/scRNAseq_pc.csv")
# mypc_scRNAseq <- read.csv("data/scRNAseq_pc.csv")
# 
# ### MDS
# mds_singlecellRNAseq <- singlecellRNAseq[,1:ncol(singlecellRNAseq)-1] %>%
#   dist() %>%
#   cmdscale(k=3) %>%
#   as.data.frame()
# 
# mds_singlecellRNAseq <- mds_singlecellRNAseq %>% mutate(label = singlecellRNAseq$label)
# # write.csv(mds_singlecellRNAseq, "data/scRNAseq_mds.csv")
# mds_scRNAseq <- read.csv("data/scRNAseq_mds.csv")
# 
# ### tsne
# # singlecellRNAseq<-read.csv("data/single-cell-rnaseq.csv")
# # singlecellRNAseq <- singlecellRNAseq[,-(which(colSums(singlecellRNAseq) == 0))]
# # singlecellRNAseq <- singlecellRNAseq %>% mutate(label=as.factor(label))
# # pc <- prcomp(mydata[,1:22215], center = F, scale. = F)
# allpcs <- pc$x[,1:50]
# allpcs <- as.data.frame(allpcs)
# allpcs$label <- as.factor(singlecellRNAseq$label)
# # write.csv(allpcs, "data/scRNAseq_allpc.csv")
# # allpcs <- read.csv("data/scRNAseq_allpc.csv")
# 
# tsne.results <- Rtsne(pc$x[,1:50], check_duplicates = FALSE, dims=3, theta=0.0, pca=FALSE, verbose=FALSE, max_iter=1000, perplexity=10)
# colnames(tsne.results$Y) <- c("V1", "V2", "V3")
# mytsnedata <- data.frame(tsne.results$Y)
# mytsnedata$label <- as.factor(singlecellRNAseq$label)
# # write.csv(mytsnedata, "data/tsne_three_scRNAseq.csv", row.names = FALSE)
# mytsne_scRNAseq <- read.csv("data/tsne_three_scRNAseq.csv")
# 
# 
# ####################################################################
# 
# # tissue <- read.csv("data/tissue.csv") # 189
# # e <- read.csv("data/expression.csv", row.names = "probeID")
# # SE <- SummarizedExperiment(assays = I(e), colData = DataFrame(tissue))
# # tissue <- tissue[[1]]
# # mydata <- data.frame(t(assay(SE)))
# #
# # NZV<- nearZeroVar(mydata, saveMetrics = TRUE)
# # mydata_filter <- mydata[, which(NZV[, "nzv"] == FALSE)]


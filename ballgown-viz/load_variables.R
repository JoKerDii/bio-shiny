rm(list= ls())

setwd("/data/zhendi/protocol/ballgown/ballgown-viz")

# initializer <- function(gown_file){
# Load R packages

library(ballgown)
# library(genefilter)
library(ggplot2)
library(tidyr)
# library(shiny)
# library(shinythemes)
library(dplyr)
# library(gplots)


load("homo_bg.rda")

pData(bg)
# bg <- get(ls()[ls() != gown_file])
group <- c(rep("experimental", 3), rep("control", 3))
pData(bg) = data.frame(id=sampleNames(bg), group=group)
pData_bg<-pData(bg)


sampleNames_bg <- pData(bg)[,1]
group <- pData(bg)[,2]

geneIDs <- ballgown::geneIDs(bg)

sampleNames_bg <- sampleNames(bg) # not required
# Analysis with Ballgown
# Extract FPKM values from the 'bg' object
fpkm = texpr(bg,meas="FPKM")
# Transform the FPKM values by adding 1 and convert to a log2 scale
fpkm = log2(fpkm+1)
# Load all attributes and gene names
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10]) # not required
# Pull the gene_expression data frame from the ballgown object
gene_expression = as.data.frame(gexpr(bg))
# Load the transcript to gene index from the ballgown object
transcript_gene_table = indexes(bg)$t2g
#Each row of data represents a transcript. Many of these transcripts represent the same gene. Determine the numbers of transcripts and unique genes
# Differential expression results
results_genes = stattest(bg, feature="gene", covariate="group", getFC=TRUE, meas="FPKM") #not required
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))


# save variables
variables <- c("geneIDs", "pData_bg", "fpkm", "bg_table", "gene_expression", "transcript_gene_table", "results_genes")
ls()
save(list=variables, file="myvariables.RData")
# rm(list=ls())
# ls()
load("myvariables.RData")
# library(dplyr)
# 
# getSamplesFromColnames<-function(fpkmDF){
#   samples<-gsub("_FPKM$","",colnames(fpkmDF)[grepl("_FPKM$",colnames(fpkmDF))])
# }
# 
# ####################################### handle "read_groups.info"
# allfiles = list(sep="\t", header=TRUE, row.names = NULL, quote="", na.string="-")
# filename = "/home/zhen.di/cw1/DGE/read_groups.info"
# full = as.data.frame(read.delim(filename))
# loadRepTable = full %>% mutate(condition=as.character(condition)) %>% mutate(rep_name = paste(condition, replicate_num, sep="_"))
# 
# ####################################### handle "var_model.info"
# allfiles = list(sep="\t", header=TRUE, row.names = NULL, quote="", na.string="-")
# filename = "/home/zhen.di/cw1/DGE/var_model.info"
# full = as.data.frame(read.delim(filename))
# loadVarModelTable = full %>% mutate(condition=as.character(condition)) 
# 
# ####################################### handle genes: five files
## geneFPKM="genes.fpkm_tracking",
# allfiles = list(sep="\t", header=TRUE, row.names = NULL, quote="", na.string="-")
# filename = "/home/zhen.di/cw1/DGE/genes.fpkm_tracking"
# allfiles$file = filename
# full = as.data.frame(do.call(read.table,allfiles))
# samples<-getSamplesFromColnames(full)
# genesTable<-full[,c(1:3,5,7:9)]
# genemelt<-melt(full,id.vars=c("tracking_id"),measure.vars=-c(1:9),variable_name="sample_name")
# colnames(genemelt)[colnames(genemelt)=='variable']<-'sample_name'
# # clean and normalize
# genemelt$measurement = ""
# 
# genemelt$measurement[grepl("_FPKM$",genemelt$sample_name)] = "fpkm"
# genemelt$measurement[grepl("_conf_lo$",genemelt$sample_name)] = "conf_lo"
# genemelt$measurement[grepl("_conf_hi$",genemelt$sample_name)] = "conf_hi"
# genemelt$measurement[grepl("_status$",genemelt$sample_name)] = "status"
# 
# genemelt$sample_name<-gsub("_FPKM$","",genemelt$sample_name)
# genemelt$sample_name<-gsub("_conf_lo$","",genemelt$sample_name)
# genemelt$sample_name<-gsub("_conf_hi$","",genemelt$sample_name)
# genemelt$sample_name<-gsub("_status$","",genemelt$sample_name)
# 
# genemelt = genemelt %>% mutate(sample_name = as.vector(sample_name))
# genemelt<-as.data.frame(dcast(genemelt,...~measurement))
# genemelt =genemelt[,c(1:2,5,3,4,6)]
# 
# # geneDiff="gene_exp.diff",  ############################################## usable
# allfiles = list(sep="\t", header=TRUE, row.names = NULL, quote="", na.string="-")
# filename = "/home/zhen.di/cw1/DGE/gene_exp.diff"
# allfiles$file = filename
# diff = as.data.frame(do.call(read.table,allfiles))
# if(dim(diff)[1]>0){
#   diff = diff %>% mutate(sample_1 = as.vector(sample_1), sample_2 = as.vector(sample_2))
#   diff = diff[, c(1,5:14)]
# }
# 
# # promoterFile="promoters.diff",
# allfiles = list(sep="\t", header=TRUE, row.names = NULL, quote="", na.string="-")
# filename = "/home/zhen.di/cw1/DGE/promoters.diff"
# allfiles$file = filename
# promoter = as.data.frame(do.call(read.table,allfiles))
# if(dim(promoter)[1]>0){
#   promoter = promoter[,c(2,5:14)]
# }
# 
# # geneCount="genes.count_tracking",
# allfiles = list(sep="\t", header=TRUE, row.names = NULL, quote="", na.string="-")
# filename = "/home/zhen.di/cw1/DGE/genes.count_tracking"
# allfiles$file = filename
# counts = as.data.frame(do.call(read.table,allfiles))
# if(dim(counts)[1]>0){
#   countmelt<-melt(counts,id.vars=c("tracking_id"),measure.vars=-c(1))
#   colnames(countmelt)[colnames(countmelt)=='variable']<-'sample_name'
#   # clear and normalize
#   countmelt$measurement = ""
#   
#   countmelt$measurement[grepl("_count$",countmelt$sample_name)] = "count"
#   countmelt$measurement[grepl("_count_variance$",countmelt$sample_name)] = "variance"
#   countmelt$measurement[grepl("_count_uncertainty_var$",countmelt$sample_name)] = "uncertainty"
#   countmelt$measurement[grepl("_count_dispersion_var$",countmelt$sample_name)] = "dispersion"
#   countmelt$measurement[grepl("_status$",countmelt$sample_name)] = "status"
#   
#   countmelt$sample_name<-gsub("_count$","",countmelt$sample_name)
#   countmelt$sample_name<-gsub("_count_variance$","",countmelt$sample_name)
#   countmelt$sample_name<-gsub("_count_uncertainty_var$","",countmelt$sample_name)
#   countmelt$sample_name<-gsub("_count_dispersion_var$","",countmelt$sample_name)
#   countmelt$sample_name<-gsub("_status$","",countmelt$sample_name)
#   
#   countmelt = countmelt %>% mutate(sample_name = as.vector(sample_name))
#   #Adjust sample names with make.db.names
#   countmelt<-as.data.frame(dcast(countmelt,...~measurement))
# }
#   
# # geneRep="genes.read_group_tracking",
# allfiles = list(sep="\t", header=TRUE, row.names = NULL, quote="", na.string="-")
# filename = "/home/zhen.di/cw1/DGE/genes.read_group_tracking"
# allfiles$file = filename
# reps = as.data.frame(do.call(read.table,allfiles))
# if(dim(reps)[1]>0){
#   reps = reps %>% mutate(condition = as.character(condition))
#   reps$rep_name<-paste(reps$condition,reps$replicate,sep="_")
#   colnames(reps)[colnames(reps)=="condition"]<-"sample_name"
# }
# 
# ####################################### handle isoforms: four files
# # isoformFPKM="isoforms.fpkm_tracking",
# allfiles = list(sep="\t", header=TRUE, row.names = NULL, quote="", na.string="-")
# filename = "/home/zhen.di/cw1/DGE/isoforms.fpkm_tracking"
# allfiles$file = filename
# full = as.data.frame(do.call(read.table,allfiles))
# samples<-getSamplesFromColnames(full)
# isoformsTable<-full[,c(1,4,5,6,2,3,7:9)]
# isoformsTable<-cbind(isoformsTable[,1:2],data.frame(CDS_id=rep("NA",dim(isoformsTable)[1])),isoformsTable[,-c(1:2)])
# 
# isoformmelt<-melt(full,id.vars=c("tracking_id"),measure.vars=-idCols,variable_name="sample_name")
# colnames(isoformmelt)[colnames(isoformmelt)=='variable']<-'sample_name'
# 
# #Clean up and normalize data
# isoformmelt$measurement = ""
# 
# isoformmelt$measurement[grepl("_FPKM$",isoformmelt$sample_name)] = "fpkm"
# isoformmelt$measurement[grepl("_conf_lo$",isoformmelt$sample_name)] = "conf_lo"
# isoformmelt$measurement[grepl("_conf_hi$",isoformmelt$sample_name)] = "conf_hi"
# isoformmelt$measurement[grepl("_status$",isoformmelt$sample_name)] = "status"
# 
# isoformmelt$sample_name<-gsub("_FPKM$","",isoformmelt$sample_name)
# isoformmelt$sample_name<-gsub("_conf_lo$","",isoformmelt$sample_name)
# isoformmelt$sample_name<-gsub("_conf_hi$","",isoformmelt$sample_name)
# isoformmelt$sample_name<-gsub("_status$","",isoformmelt$sample_name)
# 
# isoformmelt = isoformmelt %>% mutate(sample_name = as.vector(sample_name))
# isoformmelt<-as.data.frame(dcast(isoformmelt,...~measurement))
# isoformmelt = isoformmelt[,c(1:2,5,3,4,6)]
# 
# # isoformDiff="isoform_exp.diff",
# allfiles = list(sep="\t", header=TRUE, row.names = NULL, quote="", na.string="-")
# filename = "/home/zhen.di/cw1/DGE/isoform_exp.diff"
# allfiles$file = filename
# diff = as.data.frame(do.call(read.table,allfiles))
# if(dim(diff)[1]>0){
#   diff = diff %>% mutate(sample_1 = as.vector(sample_1), sample_2 = as.vector(sample_2))
#   diff = diff[,c(1,5:14)]
# }
# # isoformCount="isoforms.count_tracking",
# allfiles = list(sep="\t", header=TRUE, row.names = NULL, quote="", na.string="-")
# filename = "/home/zhen.di/cw1/DGE/isoforms.count_tracking"
# allfiles$file = filename
# counts = as.data.frame(do.call(read.table,allfiles))
# if(dim(counts)[1]>0){
#   
#   countmelt<-melt(counts,id.vars=c("tracking_id"),measure.vars=-c(1))
#   colnames(countmelt)[colnames(countmelt)=='variable']<-'sample_name'
#   
#   countmelt$measurement = ""
#   
#   countmelt$measurement[grepl("_count$",countmelt$sample_name)] = "count"
#   countmelt$measurement[grepl("_count_variance$",countmelt$sample_name)] = "variance"
#   countmelt$measurement[grepl("_count_uncertainty_var$",countmelt$sample_name)] = "uncertainty"
#   countmelt$measurement[grepl("_count_dispersion_var$",countmelt$sample_name)] = "dispersion"
#   countmelt$measurement[grepl("_status$",countmelt$sample_name)] = "status"
#   
#   countmelt$sample_name<-gsub("_count$","",countmelt$sample_name)
#   countmelt$sample_name<-gsub("_count_variance$","",countmelt$sample_name)
#   countmelt$sample_name<-gsub("_count_uncertainty_var$","",countmelt$sample_name)
#   countmelt$sample_name<-gsub("_count_dispersion_var$","",countmelt$sample_name)
#   countmelt$sample_name<-gsub("_status$","",countmelt$sample_name)
#   
#   countmelt = countmelt %>% mutate(sample_name = as.vector(sample_name))
#   countmelt = as.data.frame(dcast(countmelt,...~measurement))
#   
# }
# 
# # isoformRep="isoforms.read_group_tracking",
# allfiles = list(sep="\t", header=TRUE, row.names = NULL, quote="", na.string="-")
# filename = "/home/zhen.di/cw1/DGE/isoforms.read_group_tracking"
# allfiles$file = filename
# reps = as.data.frame(do.call(read.table,allfiles))
# if(dim(reps)[1]>0){
#   reps = reps %>% mutate(condition = as.character(condition))
#   reps$rep_name<-paste(reps$condition,reps$replicate,sep="_")
#   colnames(reps)[colnames(reps)=="condition"]<-"sample_name"
# }

#######################################################################

library(dplyr)
library(ggplot2)
library(plotly)

handle_diff <- function(file1){
  allfiles = list(sep="\t", header=TRUE, row.names = NULL, quote="", na.string="-")
  allfiles$file = file1
  diff = as.data.frame(do.call(read.table,allfiles))
  if(dim(diff)[1]>0){
    # diff = diff %>% mutate(sample_1 = as.vector(sample_1), sample_2 = as.vector(sample_2))
    diff = diff[, c(1,5:14)]
  }
  return(diff)
}

# file1 <- "./data/gene_exp.diff"
# diff <- handle_diff(file1)

# write.csv(diff, "/home/zhen.di/cw1/DGE/cuffdiff-viz/data/diff.csv", row.names = F)
# diff=read.csv("/home/zhen.di/cw1/DGE/cuffdiff-viz/data/diff.csv")

# geneRep="genes.read_group_tracking",
handle_tracking <- function(file2) {
  allfiles = list(sep="\t", header=TRUE, row.names = NULL, quote="", na.string="-")
  allfiles$file = file2
  reps = as.data.frame(do.call(read.table,allfiles))
  if(dim(reps)[1]>0){
    reps = reps %>% mutate(condition = as.character(condition))
    reps$rep_name<-paste(reps$condition,reps$replicate,sep="_")
    # colnames(reps)[colnames(reps)=="condition"]<-"sample_name"
  }
  
  return(reps)
}

file2 = "./data/genes.read_group_tracking"
reps <- handle_tracking(file2)

# write.csv(reps, "/home/zhen.di/cw1/DGE/cuffdiff-viz/data/reps.csv", row.names = F)
# reps=read.csv("/home/zhen.di/cw1/DGE/cuffdiff-viz/data/reps.csv")

handle_info <- function(file3) {
  allfiles = list(sep="\t", header=TRUE, row.names = NULL, quote="", na.string="-")
  # allfiles$file = file3
  
  info = as.data.frame(read.delim(file3))[,1:3]
  # info$rep_name<-paste(info$condition,info$replicate,sep="_")
  # if(dim(reps)[1]>0){
  #   reps = reps %>% mutate(condition = as.character(condition))
  #   reps$rep_name<-paste(reps$condition,reps$replicate,sep="_")
  #   # colnames(reps)[colnames(reps)=="condition"]<-"sample_name"
  # }
  
  return(info)
}
# file3 = "./data/read_groups.info"
# info <- handle_info(file3)

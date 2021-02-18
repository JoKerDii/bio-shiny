
# customized plotter
source('./plotTranscripts_cust.R', echo=TRUE)

# plot01 and plot #3
boxplot01 <- function(fpkm_data, transcript = TRUE) {
  if (transcript == FALSE){
    # fpkm_data <- gexpr(bg) [gene level expression]
    min_nonzero <- 1
    fpkm_long <- gather(data.frame(log2(fpkm_data+min_nonzero)), "Sample", "log2(FPKM+1)")
  } else {
    # fpkm_data <- log(texpr(bg,meas="FPKM") + 1) [transcript level expression]
    fpkm_long <- gather(data.frame(fpkm_data), "Sample", "log2(FPKM+1)")
  }
  
  index_table <- pData(bg)
  index_table$id <- paste("FPKM.", index_table$id, sep='')
  colnames(index_table) <- c("Sample", "Group")
  fpkm_long <- left_join(fpkm_long, index_table, by = "Sample")
  p <- ggplot(fpkm_long, aes(x=Sample, y=`log2(FPKM+1)`, fill = Group)) + 
    geom_boxplot() + coord_flip() +
    labs(title = "log2(FPKM+1) for Each Sample") +
    theme_bw() + 
    theme(plot.title = element_text(face="bold", size=18, vjust=0.5, hjust=0.5),
          axis.title.x = element_text(face="bold", size=15, vjust=0.5), 
          axis.title.y = element_text(face="bold", size=15, vjust=0.5),
          axis.text.x = element_text(size=12, colour="black", face="bold"),
          axis.text.y = element_text(size=12, colour="black", face="bold"),
          legend.text = element_text(size=15, colour="black", face="bold"),
          legend.title = element_text(size=15, colour="black", face="bold"),
          # legend.background = element_rect(fill = "transparent"),
          # legend.key = element_rect(fill = "transparent"),
          # legend.justification = c(1, 0), legend.position = c(0.99, 0.01),
          plot.margin = margin(10, 100, 10, 100),
          plot.caption = element_text(size=12, colour="gray75", face="italic", hjust = 1, vjust = 1))
  return(p)
}

# plot02
boxplot02 <- function(fpkm_data, Gene) {
  colnames(fpkm_data) <- group
  fpkm_tibble <- tibble(FPKM = fpkm_data[round(Gene),], Group = group)
  p <- ggplot(data.frame(fpkm_tibble), aes(x=Group, y=FPKM, color = Group)) + 
    geom_boxplot() +
    geom_jitter(width=0.15) + labs(title = "log2(FPKM+1) of a Single Gene for Each Sample") +
    theme_bw() + 
    theme(plot.title = element_text(face="bold", size=18, vjust=0.5, hjust=0.5),
          axis.title.x = element_text(face="bold", size=15, vjust=0.5), 
          axis.title.y = element_text(face="bold", size=15, vjust=0.5),
          axis.text.x = element_text(size=12, colour="black", face="bold"),
          axis.text.y = element_text(size=12, colour="black", face="bold"),
          legend.text = element_text(size=15, colour="black", face="bold"),
          legend.title = element_text(size=15, colour="black", face="bold"),
          # legend.background = element_rect(fill = "transparent"),
          # legend.key = element_rect(fill = "transparent"),
          # legend.justification = c(1, 0), legend.position = c(0.99, 0.01),
          plot.margin = margin(10, 100, 10, 100),
          plot.caption = element_text(size=12, colour="gray75", face="italic", hjust = 1, vjust = 1))
  return(p)
}

# plot03
trans_structure <- function(Gene, Sample){
  p = plotTranscripts_cust(ballgown::geneIDs(bg)[Gene], bg, main=paste('Transcript Structure in Sample',Sample), sample=Sample, labelTranscripts=TRUE)
  return(p)
}

# plot 04
diff_exp_hist <- function(results_genes){
  sig=which(results_genes$pval<0.05)
  results_genes[,"de"] = log2(results_genes[,"fc"])
  frame <- data.frame(value = results_genes[sig,"de"])
  p <- ggplot(frame) +
    geom_histogram(aes(x = value), bins = 100, fill='#5d6fe3') + 
    labs(title = "Distribution of Differential Expression Values") + 
    xlab("log2(fold change) experimental vs control")+ ylab("Frequency")+
    theme_bw() + 
    theme(plot.title = element_text(face="bold", size=18, vjust=0.5, hjust=0.5),
          axis.title.x = element_text(face="bold", size=15, vjust=0.5), 
          axis.title.y = element_text(face="bold", size=15, vjust=0.5),
          axis.text.x = element_text(size=12, colour="black", face="bold"),
          axis.text.y = element_text(size=12, colour="black", face="bold"),
          legend.text = element_text(size=15, colour="black", face="bold"),
          legend.title = element_text(size=15, colour="black", face="bold"),
          # legend.background = element_rect(fill = "transparent"),
          # legend.key = element_rect(fill = "transparent"),
          # legend.justification = c(1, 0), legend.position = c(0.99, 0.01),
          plot.margin = margin(10, 100, 10, 100),
          plot.caption = element_text(size=12, colour="gray75", face="italic", hjust = 1, vjust = 1))
  return(p)
}


# plot 05
sig_diff <- function(gene_expression, results_genes){
  gene_expression = as.data.frame(gexpr(bg))
  sig=which(results_genes$pval<0.05)
  # generate label names
  index_table <- pData(bg)
  index_table$id <- paste("FPKM", index_table$id, sep='.')
  colnames(index_table) <- c("Sample", "Group")
  shortNames <- as.character(index_table$Group[match(colnames(gene_expression), index_table$Sample)])
  # colnames(gene_expression) <- shortNames
  # unique(shortNames)[1]
  index_table$Group <- as.character(index_table$Group)
  group1 <- unique(index_table$Group)[1]
  group2 <- unique(index_table$Group)[2]
  gene_expression[,group1]=apply(gene_expression[,index_table$Group == group1], 1, mean)
  gene_expression[,group2]=apply(gene_expression[,index_table$Group == group2], 1, mean)
  
  # data frame
  min_nonzero = 1
  x=log2(gene_expression[,group1]+min_nonzero) 
  y=log2(gene_expression[,group2]+min_nonzero)
  xsig=x[sig]
  ysig=y[sig]
  frame = data.frame(x=x, y=y)
  frame2= data.frame(x = xsig, y = ysig)
  # plot
  p <- ggplot(data = frame) + 
    geom_point(aes(x = x, y = y), pch=20, cex=0.25) + 
    geom_point(frame2, mapping = aes(x = x, y = y), col="#c9388d", pch=20, cex=0.5) +
    xlab(paste(group1," FPKM (log2)")) + ylab(paste(group2," FPKM (log2)")) + 
    labs(title="Experimental vs Control: FPKMs")+
    annotate(
      geom = "text", x = 0.1, y = max(y)-1, 
      label = "Significant", hjust = 0, vjust = 1, size = 6, col = "#c9388d"
    ) +
    theme_bw() + 
    theme(plot.title = element_text(face="bold", size=18, vjust=0.5, hjust=0.5),
          axis.title.x = element_text(face="bold", size=15, vjust=0.5), 
          axis.title.y = element_text(face="bold", size=15, vjust=0.5),
          axis.text.x = element_text(size=12, colour="black", face="bold"),
          axis.text.y = element_text(size=12, colour="black", face="bold"),
          legend.text = element_text(size=15, colour="black", face="bold"),
          legend.title = element_text(size=15, colour="black", face="bold"),
          # legend.background = element_rect(fill = "transparent"),
          # legend.key = element_rect(fill = "transparent"),
          # legend.justification = c(1, 0), legend.position = c(0.99, 0.01),
          plot.margin = margin(10, 100, 10, 100),
          plot.caption = element_text(size=12, colour="gray75", face="italic", hjust = 1, vjust = 1))
  return(p)
  
}
# plot 06
heatmap <- function(gene_expression, results_genes){
  results_genes[,"de"] = log2(results_genes[,"fc"])
  #Get the gene symbols for the top N (according to corrected p-value) and display them on the plot
  topn = order(abs(results_genes[sig,"fc"]), decreasing=TRUE)[1:25]
  topn = order(results_genes[sig,"qval"])[1:25]
  
  #Each should be significant with a log2 fold-change >= 2
  sigpi = which(results_genes[,"pval"]<0.05)
  sigp = results_genes[sigpi,]
  sigde = which(abs(sigp[,"de"]) >= 2)
  sig_tn_de = sigp[sigde,]
  
  mydist=function(c) {dist(c,method="euclidian")}
  myclust=function(c) {hclust(c,method="average")}
  main_title="Significance of Differential Expression"
  par(cex.main=1.5)
  sig_genes_de=sig_tn_de[,"id"]
  sig_gene_names_de=sig_tn_de[,"gene_name"]
  data=log2(as.matrix(gene_expression[as.vector(sig_genes_de),])+1)
  p <- heatmap.2(data, 
                 hclustfun=myclust, 
                 distfun=mydist, 
                 na.rm = TRUE, 
                 scale="none", 
                 dendrogram="both", 
                 margins=c(10,4), 
                 Rowv=TRUE, Colv=TRUE, 
                 symbreaks=FALSE, key=TRUE, symkey=FALSE, 
                 density.info="none", trace="none", 
                 main=main_title, 
                 cexRow=1, cexCol=1.5, 
                 labRow=sig_gene_names_de,
                 col=rev(heat.colors(75)))
  return(p)
  
}


# plot #1
transcript_per_gene <- function(transcript_gene_table) {
  
  counts=table(transcript_gene_table[,"g_id"])
  c_one = length(which(counts == 1))
  c_more_than_one = length(which(counts > 1))
  c_max = max(counts)
  legend_text = c(paste("Genes with one transcript =", c_one, "\nGenes with more than one transcript =", c_more_than_one, "\nMax transcripts for single gene = ", c_max))
  p <- ggplot(data.frame(counts), aes(x=Freq)) + 
    geom_histogram(color="bisque4", fill="bisque4", bins=100) +
    labs(title = "Distribution of Transcript Count per Gene", subtitle = legend_text) +
    xlab("Count") + ylab("Frequency") +
    theme_bw() + 
    theme(plot.title = element_text(face="bold", size=18, vjust=0.5, hjust=0.5),
          plot.subtitle = element_text(size=13),
          axis.title.x = element_text(face="bold", size=15, vjust=0.5), 
          axis.title.y = element_text(face="bold", size=15, vjust=0.5),
          axis.text.x = element_text(size=12, colour="black", face="bold"),
          axis.text.y = element_text(size=12, colour="black", face="bold"),
          legend.text = element_text(size=15, colour="black", face="bold"),
          # legend.title = element_text(size=15, colour="black", face="bold"),
          # legend.background = element_rect(fill = "transparent"),
          # legend.key = element_rect(fill = "transparent"),
          # legend.justification = c(1, 0), legend.position = c(0.99, 0.01),
          plot.margin = margin(10, 100, 10, 100),
          plot.caption = element_text(size=12, colour="gray75", face="italic", hjust = 1, vjust = 1))
  return(p)
}

# Plot #2
transcript_size_dist <- function(bg_table){
  full_table <- bg_table
  p <- ggplot(full_table, aes(x=length)) + 
    geom_histogram(color="steelblue", fill="steelblue", bins=150) +
    labs(title = "Distribution of Transcript Lengths") + xlab("Transcript length (bp)") + ylab("Frequency") +
    theme_bw() + 
    theme(plot.title = element_text(face="bold", size=18, vjust=0.5, hjust=0.5),
          axis.title.x = element_text(face="bold", size=15, vjust=0.5), 
          axis.title.y = element_text(face="bold", size=15, vjust=0.5),
          axis.text.x = element_text(size=12, colour="black", face="bold"),
          axis.text.y = element_text(size=12, colour="black", face="bold"),
          legend.text = element_text(size=15, colour="black", face="bold"),
          # legend.title = element_text(size=15, colour="black", face="bold"),
          # legend.background = element_rect(fill = "transparent"),
          # legend.key = element_rect(fill = "transparent"),
          # legend.justification = c(1, 0), legend.position = c(0.99, 0.01),
          plot.margin = margin(10, 100, 10, 100),
          plot.caption = element_text(size=12, colour="gray75", face="italic", hjust = 1, vjust = 1))
  return(p)
}

# plot #4
pairOfReplicates <- function(gene_expression, choice1 = sampleNames(bg)[1], choice2 = sampleNames(bg)[2]) {
  choice1 <- paste("FPKM",choice1, sep='.')
  choice2 <- paste("FPKM",choice2, sep='.')
  x = gene_expression[,choice1]
  y = gene_expression[,choice2]
  min_nonzero = 1
  frame = data.frame(x = log2(x+min_nonzero), y = log2(y+min_nonzero))
  p <- ggplot(frame, aes(x = x, y = y)) +
    geom_point(pch=16, col="#2ed9d3", cex=0.25) +
    xlab(choice1) + ylab(choice2) + labs(title = "Comparison of Expression Values") +
    theme_bw() + 
    theme(plot.title = element_text(face="bold", size=18, vjust=0.5, hjust=0.5),
          axis.title.x = element_text(face="bold", size=15, vjust=0.5), 
          axis.title.y = element_text(face="bold", size=15, vjust=0.5),
          axis.text.x = element_text(size=12, colour="black", face="bold"),
          axis.text.y = element_text(size=12, colour="black", face="bold"),
          legend.text = element_text(size=15, colour="black", face="bold"),
          legend.title = element_text(size=15, colour="black", face="bold"),
          # legend.background = element_rect(fill = "transparent"),
          # legend.key = element_rect(fill = "transparent"),
          # legend.justification = c(1, 0), legend.position = c(0.99, 0.01),
          plot.margin = margin(10, 100, 10, 100),
          plot.caption = element_text(size=12, colour="gray75", face="italic", hjust = 1, vjust = 1))
  return(p)
}

# plot #5
MDSdistance <- function(gene_expression){
  data_columns <- colnames(gene_expression)
  gene_expression_sum=apply(gene_expression[,data_columns], 1, sum)
  #Identify the genes with a grand sum FPKM of at least 5 -  filter out the genes with very low expression across the board
  i = which(gene_expression_sum > 5)
  #Calculate the correlation between all pairs of data
  r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
  #distance
  d=1-r
  mds=cmdscale(d, k=2, eig=TRUE)
  data_colors=c("tomato1","tomato2","tomato3","royalblue1","royalblue2","royalblue3")
  # generate label names
  index_table <- pData(bg)
  index_table$id <- paste("FPKM", index_table$id, sep='.')
  colnames(index_table) <- c("Sample", "Group")
  shortNames <- as.character(index_table$Group[match(colnames(gene_expression), index_table$Sample)])
  
  frame = data.frame(mds$points)
  p <- ggplot(data = frame, aes(x = frame[,1], y = frame[,2])) +
    geom_point(col="grey", cex=2, pch=16) +
    xlim(-0.08,0.08) + ylim(-0.05,0.05) + xlab("") + ylab("")+
    geom_text(label = shortNames, col = data_colors) + 
    labs(title = "MDS Distance Plot (all non-zero genes)") +
    theme_bw() + 
    theme(plot.title = element_text(face="bold", size=18, vjust=0.5, hjust=0.5),
          axis.title.x = element_text(face="bold", size=15, vjust=0.5), 
          axis.title.y = element_text(face="bold", size=15, vjust=0.5),
          axis.text.x = element_text(size=12, colour="black", face="bold"),
          axis.text.y = element_text(size=12, colour="black", face="bold"),
          legend.text = element_text(size=15, colour="black", face="bold"),
          legend.title = element_text(size=15, colour="black", face="bold"),
          # legend.background = element_rect(fill = "transparent"),
          # legend.key = element_rect(fill = "transparent"),
          # legend.justification = c(1, 0), legend.position = c(0.99, 0.01),
          plot.margin = margin(10, 100, 10, 100),
          plot.caption = element_text(size=12, colour="gray75", face="italic", hjust = 1, vjust = 1))
  return(p)
}

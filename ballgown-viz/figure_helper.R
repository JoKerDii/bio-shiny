
# customized plotter
source('./plotTranscripts_cust.R', echo=TRUE)

# plot01
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

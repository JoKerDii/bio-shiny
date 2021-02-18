makeTable1 <- function(gene_expression){
  data_columns <- colnames(gene_expression)
  gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)
  
  #Identify the genes with a grand sum FPKM of at least 5 - filter out the genes with very low expression across the board
  i = which(gene_expression[,"sum"] > 5)
  
  #Calculate the correlation between all pairs of data
  r=round(cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson"),6)
  
  #Print out these correlation values
  return(r)
}

makeTable2 <- function(results_genes){
  sig = which(results_genes$pval<0.05)
  results_genes[,"de"] = log2(results_genes[,"fc"])
  #Get the gene symbols for the top N (according to corrected p-value) and display them on the plot
  topn = order(abs(results_genes[sig,"fc"]), decreasing=TRUE)[1:25]
  topn = order(results_genes[sig,"qval"])[1:25]

  #Each should be significant with a log2 fold-change >= 2
  sigpi = which(results_genes[,"pval"]<0.05)
  sigp = results_genes[sigpi,]
  sigde = which(abs(sigp[,"de"]) >= 2)
  sig_tn_de = sigp[sigde,]

  #Order the output by or p-value and then break ties using fold-change
  o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"de"]), decreasing=FALSE)

  output = sig_tn_de[o,c("gene_name","id","fc","pval","qval","de")]
  output[,c("fc","pval","qval","de")] = round(output[,c("fc","pval","qval","de")],5)
  return(output)
}
options(shiny.maxRequestSize = 30*1024^2)

library(ggplot2)
library(tidyr)
library(shinythemes)
# unloadNamespace("shiny")
require(shiny)
library(RColorBrewer)
library(dplyr)
library(data.table)
library(gplots)


# load data
source("./handle_cuffdiff.R")

# ##### new data
# if (is.null(input$file2) == FALSE & is.null(input$file3) == FALSE) {
#   diff <- handle_diff(input$file2)
#   reps <- handle_tracking(input$file3)
#   index_table = data.frame(Sample = unique(reps$rep_name))
#   index_table$Group = index_table %>% separate(Sample, "Group")
#   newreps <- reps %>%
#     select(c('tracking_id','rep_name' ,'FPKM')) %>%
#     spread('rep_name', 'FPKM')
#   rownames(newreps) <- newreps$tracking_id
#   newreps <- newreps %>% select(-c('tracking_id'))
# }


# Define server logic
server <- function(input, output) {
  # load data
  diff <- read.csv("./data/diff.csv")
  reps <- read.csv("./data/reps.csv")
  index_table = data.frame(Sample = unique(reps$rep_name))
  index_table$Group = index_table %>% separate(Sample, "Group")
  newreps <- reps %>%
    select(c('tracking_id','rep_name' ,'FPKM')) %>%
    spread('rep_name', 'FPKM')
  rownames(newreps) <- newreps$tracking_id
  newreps <- newreps %>% select(-c('tracking_id'))
  
  
  ################################################### information about uploaded file ended
  # output$contents <- renderTable({
  #   req(input$file1)
  #   tryCatch(
  #     {
  #       info <- handle_info(input$file1$datapath)
  #     },
  #     error = function(e) {
  #       # return a safeError if a parsing error occurs
  #       stop(safeError(e))
  #     }
  #   )
  #   head(info)
  # })
  output$contents <- DT::renderDataTable(DT::datatable({
    req(input$file1)
    
    # print(input$file1$datapath)
    tryCatch(
      {
        # print(input$file1$datapath)
        info <- handle_info(input$file1$datapath)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    head(info)
  }))
  
  # # load new data
  # if (is.null(input$file2$datapath) == FALSE & is.null(input$file3$datapath) == FALSE) {
  #   diff <- handle_diff(input$file2$datapath)
  #   reps <- handle_tracking(input$file3$datapath)
  #   index_table = data.frame(Sample = unique(reps$rep_name))
  #   index_table$Group = index_table %>% separate(Sample, "Group")
  #   newreps <- reps %>%
  #     select(c('tracking_id','rep_name' ,'FPKM')) %>%
  #     spread('rep_name', 'FPKM')
  #   rownames(newreps) <- newreps$tracking_id
  #   newreps <- newreps %>% select(-c('tracking_id'))
  # }
  ################################################### information about uploaded file ended
  
  output$logFPKM <- renderPlotly({
    if (is.null(input$file2$datapath) == FALSE & is.null(input$file3$datapath) == FALSE) {
      # load new data
      # DIFF <- handle_diff(input$file2$datapath)
      REPS <- handle_tracking(input$file3$datapath)
      # INDEX_TABLE = data.frame(Sample = unique(REPS$rep_name))
      # INDEX_TABLE$Group = INDEX_TABLE %>% separate(Sample, "Group")
      # NEWREPS <- REPS %>%
      #   select(c('tracking_id','rep_name' ,'FPKM')) %>%
      #   spread('rep_name', 'FPKM')
      # rownames(NEWREPS) <- NEWREPS$tracking_id
      # NEWREPS <- NEWREPS %>% select(-c('tracking_id'))
      
      # plot
      p <- ggplot(REPS, aes(x=rep_name, y=log2(FPKM+1),fill = condition)) +
        geom_boxplot()  +
        xlab("Sample") +
        labs(title = "log2(FPKM+1) for Each Sample") +
        # labs(title = "log10(FPKM) for Each Sample") +
        scale_fill_brewer(palette="Set2")
      
      fig <- ggplotly(p)
      fig
    }else{
      p <- ggplot(reps, aes(x=rep_name, y=log2(FPKM+1),fill = condition)) +
        geom_boxplot()  +
        xlab("Sample") +
        labs(title = "log2(FPKM+1) for Each Sample") +
        # labs(title = "log10(FPKM) for Each Sample") +
        scale_fill_brewer(palette="Set2")
      
      fig <- ggplotly(p)
      fig
    }
    
  })
  
  output$logFPKM_gene <- renderPlot({
    SINGLE_ROW = input$Gene1

    
    if (is.null(input$file2$datapath) == FALSE & is.null(input$file3$datapath) == FALSE) {
      # load new data
      # DIFF <- handle_diff(input$file2$datapath)
      REPS <- handle_tracking(input$file3$datapath)
      INDEX_TABLE = data.frame(Sample = unique(REPS$rep_name))
      INDEX_TABLE$Group = INDEX_TABLE %>% separate(Sample, "Group")
      NEWREPS <- REPS %>%
        select(c('tracking_id','rep_name' ,'FPKM')) %>%
        spread('rep_name', 'FPKM')
      rownames(NEWREPS) <- NEWREPS$tracking_id
      NEWREPS <- NEWREPS %>% select(-c('tracking_id'))
      
      # plot
      newreps_t <- transpose(NEWREPS[SINGLE_ROW,])
      newreps_t$group <- colnames(NEWREPS[SINGLE_ROW,])
      colnames(newreps_t) <- c('FPKM', 'group')
      newreps_t_g <- newreps_t %>% 
        separate(group, "group")
      
      # fpkm_tibble <- tibble(FPKM = newreps[SINGLE_ROW,], group = colnames(newreps))
      p <- ggplot(newreps_t_g, aes(x=group, y=log2(FPKM+1), color = group)) +
        geom_boxplot() +
        geom_jitter(width=0.15) +
        labs(title = paste("log2(FPKM+1) of ",rownames(newreps[SINGLE_ROW,]), " Gene for Each Condition")) +
        scale_fill_brewer(palette="Set2")+
        theme(plot.title = element_text(size=18),
              axis.title.x = element_text(size=15, vjust=0.5),
              axis.title.y = element_text(size=15, vjust=0.5),
              axis.text.x = element_text(size=12, colour="black"),
              axis.text.y = element_text(size=12, colour="black"),
              legend.text = element_text(size=15, colour="black"),
              legend.title = element_text(size=15, colour="black"),
              # legend.background = element_rect(fill = "transparent"),
              # legend.key = element_rect(fill = "transparent"),
              # legend.justification = c(1, 0), legend.position = c(0.99, 0.01),
              # plot.margin = margin(10, 100, 10, 100),
              plot.caption = element_text(size=12, colour="gray75", face="italic", hjust = 1, vjust = 1))
      p

    }else{
      newreps_t <- transpose(newreps[SINGLE_ROW,])
      newreps_t$group <- colnames(newreps[SINGLE_ROW,])
      colnames(newreps_t) <- c('FPKM', 'group')
      newreps_t_g <- newreps_t %>% 
        separate(group, "group")
      
      # fpkm_tibble <- tibble(FPKM = newreps[SINGLE_ROW,], group = colnames(newreps))
      p <- ggplot(newreps_t_g, aes(x=group, y=log2(FPKM+1), color = group)) +
        geom_boxplot() +
        geom_jitter(width=0.15) +
        labs(title = paste("log2(FPKM+1) of ",rownames(newreps[SINGLE_ROW,]), " Gene for Each Condition")) +
        scale_fill_brewer(palette="Set2")+
        theme(plot.title = element_text(size=18),
              axis.title.x = element_text(size=15, vjust=0.5),
              axis.title.y = element_text(size=15, vjust=0.5),
              axis.text.x = element_text(size=12, colour="black"),
              axis.text.y = element_text(size=12, colour="black"),
              legend.text = element_text(size=15, colour="black"),
              legend.title = element_text(size=15, colour="black"),
              # legend.background = element_rect(fill = "transparent"),
              # legend.key = element_rect(fill = "transparent"),
              # legend.justification = c(1, 0), legend.position = c(0.99, 0.01),
              # plot.margin = margin(10, 100, 10, 100),
              plot.caption = element_text(size=12, colour="gray75", face="italic", hjust = 1, vjust = 1))
      p
    }
    
    
  })
  # 
  output$csDensity <- renderPlotly({
    if (is.null(input$file2$datapath) == FALSE & is.null(input$file3$datapath) == FALSE) {
      # load new data
      # DIFF <- handle_diff(input$file2$datapath)
      REPS <- handle_tracking(input$file3$datapath)
      # INDEX_TABLE = data.frame(Sample = unique(REPS$rep_name))
      # INDEX_TABLE$Group = INDEX_TABLE %>% separate(Sample, "Group")
      # NEWREPS <- REPS %>%
      #   select(c('tracking_id','rep_name' ,'FPKM')) %>%
      #   spread('rep_name', 'FPKM')
      # rownames(NEWREPS) <- NEWREPS$tracking_id
      # NEWREPS <- NEWREPS %>% select(-c('tracking_id'))
      
      # plot
      p<-ggplot(REPS)+
        geom_density(aes(x= log10(FPKM),group=rep_name,color=condition,fill=condition),alpha=I(1/3)) +
        labs(title = "Density Plot")
      fig <- ggplotly(p)
      fig
    }else{
      p<-ggplot(reps)+
        geom_density(aes(x= log10(FPKM),group=rep_name,color=condition,fill=condition),alpha=I(1/3)) +
        labs(title = "Density Plot")
      fig <- ggplotly(p)
      fig
    }
  })
  
  output$scatter <- renderPlotly({
    if (is.null(input$file2$datapath) == FALSE & is.null(input$file3$datapath) == FALSE) {
      # load new data
      DIFF <- handle_diff(input$file2$datapath)
      REPS <- handle_tracking(input$file3$datapath)
      INDEX_TABLE = data.frame(Sample = unique(REPS$rep_name))
      INDEX_TABLE$Group = INDEX_TABLE %>% separate(Sample, "Group")
      NEWREPS <- REPS %>%
        select(c('tracking_id','rep_name' ,'FPKM')) %>%
        spread('rep_name', 'FPKM')
      rownames(NEWREPS) <- NEWREPS$tracking_id
      NEWREPS <- NEWREPS %>% select(-c('tracking_id'))
      
      # plot
      sig=which(DIFF$significant == "yes")
      group1 <- unique(INDEX_TABLE$Group)[1,]
      group2 <- unique(INDEX_TABLE$Group)[2,]
      
      NEWREPS[,group1]=apply(NEWREPS[,index_table$Group == group1], 1, mean)
      NEWREPS[,group2]=apply(NEWREPS[,index_table$Group == group2], 1, mean)
      
      # data frame
      x=log2(NEWREPS[,group1]+1)
      y=log2(NEWREPS[,group2]+1)
      xsig=x[sig]
      ysig=y[sig]
      frame = data.frame(x=x, y=y)
      frame2= data.frame(x = xsig, y = ysig)
      # plot
      p <- ggplot(data = frame) +
        geom_point(aes(x = x, y = y), pch=20, cex=0.25, alpha = I(1/3)) +
        geom_point(frame2, mapping = aes(x = x, y = y), col="#c9388d", pch=20, cex=0.5, alpha = I(1/3)) +
        xlab(paste(group1," FPKM (log2)")) + ylab(paste(group2," FPKM (log2)")) +
        labs(title="Experimental vs Control: FPKMs")+
        annotate(
          geom = "text", x = 1, y = max(y)-1,
          label = "Significant", hjust = 0, vjust = 1, size = 5, col = "#c9388d"
        ) 
      ggplotly(p)
    }else{
      
      sig=which(diff$significant == "yes")
      group1 <- unique(index_table$Group)[1,]
      group2 <- unique(index_table$Group)[2,]
  
      newreps[,group1]=apply(newreps[,index_table$Group == group1], 1, mean)
      newreps[,group2]=apply(newreps[,index_table$Group == group2], 1, mean)
      
      # data frame
      x=log2(newreps[,group1]+1)
      y=log2(newreps[,group2]+1)
      xsig=x[sig]
      ysig=y[sig]
      frame = data.frame(x=x, y=y)
      frame2= data.frame(x = xsig, y = ysig)
      # plot
      p <- ggplot(data = frame) +
        geom_point(aes(x = x, y = y), pch=20, cex=0.25, alpha = I(1/3)) +
        geom_point(frame2, mapping = aes(x = x, y = y), col="#c9388d", pch=20, cex=0.5, alpha = I(1/3)) +
        xlab(paste(group1," FPKM (log2)")) + ylab(paste(group2," FPKM (log2)")) +
        labs(title="Experimental vs Control: FPKMs")+
        annotate(
          geom = "text", x = 1, y = max(y)-1,
          label = "Significant", hjust = 0, vjust = 1, size = 5, col = "#c9388d"
        ) 
      ggplotly(p)
    }  
  })

  output$Heatmap <- renderPlot({
    
    Colors =c("#27408B","#3A5FCD","#7995EA","#B1C3F8","#D9DADF","#F4F4A0","#FFD365","#FF9C5B","#FF704D","#FF3D0D","#A02422")

    
    if (is.null(input$file2$datapath) == FALSE & is.null(input$file3$datapath) == FALSE) {
      # load new data
      DIFF <- handle_diff(input$file2$datapath)
      REPS <- handle_tracking(input$file3$datapath)
      INDEX_TABLE = data.frame(Sample = unique(REPS$rep_name))
      INDEX_TABLE$Group = INDEX_TABLE %>% separate(Sample, "Group")
      NEWREPS <- REPS %>%
        select(c('tracking_id','rep_name' ,'FPKM')) %>%
        spread('rep_name', 'FPKM')
      rownames(NEWREPS) <- NEWREPS$tracking_id
      NEWREPS <- NEWREPS %>% select(-c('tracking_id'))
      
      # plot
      sigpi = which(DIFF[,"p_value"]<0.05)
      topn = order(abs(DIFF[sigpi,"log2.fold_change."]), decreasing=TRUE)[1:25]
      # topn = order(results_genes[sigpi,"qval"])[1:25]
      sigp = DIFF[sigpi,]
      sigde = which(abs(sigp[,"log2.fold_change."]) >= 2)
      sig_tn_de = sigp[sigde,]
      mydist=function(c) {dist(c,method="euclidian")}
      myclust=function(c) {hclust(c,method="average")}
      main_title="Heatmap of Differential Expression"
      # par(cex.main=1.3)
      sig_genes_de=sig_tn_de[,"test_id"]
      sig_gene_names_de=sig_tn_de[,"test_id"]
      newreps_id <- REPS %>%
        select(c('tracking_id','rep_name' ,'FPKM')) %>%
        spread('rep_name', 'FPKM') 
      newreps_sig <- newreps_id[which(newreps_id$tracking_id %in% sig_genes_de),]
      rownames(newreps_sig) <- newreps_sig$tracking_id
      newreps_sig <- newreps_sig %>%
        select(-c("tracking_id"))
      heatdata=as.matrix(log2(newreps_sig+1))
      
      p = heatmap.2(heatdata,
                    hclustfun=myclust,
                    distfun=mydist,
                    na.rm = TRUE,
                    scale="none",
                    dendrogram="both",
                    # margins=c(12,9),
                    Rowv=TRUE, Colv=TRUE,
                    symbreaks=FALSE, key=TRUE, symkey=FALSE,
                    density.info="none", trace="none",
                    main=main_title,
                    cexRow=1, cexCol=1.5,
                    key.title = "log(FPKMs+1)",
                    labRow=sig_gene_names_de,
                    # col=rev(heat.colors(75)),
                    col=Colors)
      p
    }else{
      sigpi = which(diff[,"p_value"]<0.05)
      topn = order(abs(diff[sigpi,"log2.fold_change."]), decreasing=TRUE)[1:25]
      # topn = order(results_genes[sigpi,"qval"])[1:25]
      sigp = diff[sigpi,]
      sigde = which(abs(sigp[,"log2.fold_change."]) >= 2)
      sig_tn_de = sigp[sigde,]
      mydist=function(c) {dist(c,method="euclidian")}
      myclust=function(c) {hclust(c,method="average")}
      main_title="Heatmap of Differential Expression"
      # par(cex.main=1.3)
      sig_genes_de=sig_tn_de[,"test_id"]
      sig_gene_names_de=sig_tn_de[,"test_id"]
      newreps_id <- reps %>%
        select(c('tracking_id','rep_name' ,'FPKM')) %>%
        spread('rep_name', 'FPKM') 
      newreps_sig <- newreps_id[which(newreps_id$tracking_id %in% sig_genes_de),]
      rownames(newreps_sig) <- newreps_sig$tracking_id
      newreps_sig <- newreps_sig %>%
        select(-c("tracking_id"))
      heatdata=as.matrix(log2(newreps_sig+1))
      
      p = heatmap.2(heatdata,
                    hclustfun=myclust,
                    distfun=mydist,
                    na.rm = TRUE,
                    scale="none",
                    dendrogram="both",
                    # margins=c(12,9),
                    Rowv=TRUE, Colv=TRUE,
                    symbreaks=FALSE, key=TRUE, symkey=FALSE,
                    density.info="none", trace="none",
                    main=main_title,
                    cexRow=1, cexCol=1.5,
                    key.title = "log(FPKMs+1)",
                    labRow=sig_gene_names_de,
                    # col=rev(heat.colors(75)),
                    col=Colors)
      p
    }
  })

  output$volcano <- renderPlotly({
    if (is.null(input$file2$datapath) == FALSE & is.null(input$file3$datapath) == FALSE) {
      # load new data
      DIFF <- handle_diff(input$file2$datapath)
      # REPS <- handle_tracking(input$file3$datapath)
      # INDEX_TABLE = data.frame(Sample = unique(REPS$rep_name))
      # INDEX_TABLE$Group = INDEX_TABLE %>% separate(Sample, "Group")
      # NEWREPS <- REPS %>%
      #   select(c('tracking_id','rep_name' ,'FPKM')) %>%
      #   spread('rep_name', 'FPKM')
      # rownames(NEWREPS) <- NEWREPS$tracking_id
      # NEWREPS <- NEWREPS %>% select(-c('tracking_id'))
      
      # plot
      p <- ggplot(DIFF) +
        geom_point(aes(x=`log2.fold_change.`,y=-log10(p_value),color=significant),size=1.2, alpha = I(1/3)) +
        labs(title="Volcano Plot") + 
        scale_colour_manual(values = c("black","red"))
      ggplotly(p)
    }else{
      p <- ggplot(diff) +
        geom_point(aes(x=`log2.fold_change.`,y=-log10(p_value),color=significant),size=1.2, alpha = I(1/3)) +
        labs(title="Volcano Plot") + 
        scale_colour_manual(values = c("black","red"))
      ggplotly(p)
    }
  })

  output$CorrelationMatrix <- renderPlot({
    ### Plot correlation matrix
    if (is.null(input$file2$datapath) == FALSE & is.null(input$file3$datapath) == FALSE) {
      # load new data
      DIFF <- handle_diff(input$file2$datapath)
      REPS <- handle_tracking(input$file3$datapath)
      INDEX_TABLE = data.frame(Sample = unique(REPS$rep_name))
      INDEX_TABLE$Group = INDEX_TABLE %>% separate(Sample, "Group")
      NEWREPS <- REPS %>%
        select(c('tracking_id','rep_name' ,'FPKM')) %>%
        spread('rep_name', 'FPKM')
      rownames(NEWREPS) <- NEWREPS$tracking_id
      NEWREPS <- NEWREPS %>% select(-c('tracking_id'))
      
      # plot
      sig=which(DIFF$significant == "yes")
      data_columns <- colnames(NEWREPS)
      NEWREPS[,"sum"]=apply(NEWREPS[,data_columns], 1, sum)
      
      #Identify the genes with a grand sum FPKM of at least 5 - we will filter out the genes with very low expression across the board
      i = which(NEWREPS[,"sum"] > 5)
      #Calculate the correlation between all pairs of data
      r=round(cor(NEWREPS[i,data_columns], use="pairwise.complete.obs", method="pearson"),6)
      
      dist.df = reshape2::melt(as.matrix(r),varnames=c("X1","X2"))
      dist.df$value<-as.numeric(format(dist.df$value,digits=3))
      # labels = labels(obj.dists)
      g = ggplot(dist.df, aes(x=X1, y=X2, fill=value)) +
        geom_tile(color="black") +
        geom_text(aes(label=value)) +
        xlab("") + ylab("") +
        labs(title="Correlation Matrix") +
        theme(axis.text.x=element_text(angle=-90, hjust=0), axis.text.y=element_text(angle=0, hjust=1))+
        theme(panel.grid.minor=element_line(colour=NA), panel.grid.major=element_line(colour=NA), panel.background=element_rect(fill=NA, colour=NA))+
        # scale_fill_gradient2(low=heatscale[1], mid=heatscale[2], high=heatscale[3])
        scale_fill_gradientn(colours = brewer.pal(8,"YlGnBu"))+
        theme(plot.title = element_text(size=18),
              axis.title.x = element_text(size=15, vjust=0.5),
              axis.title.y = element_text(size=15, vjust=0.5),
              axis.text.x = element_text(size=12, colour="black"),
              axis.text.y = element_text(size=12, colour="black"),
              legend.text = element_text(size=15, colour="black"),
              legend.title = element_text(size=15, colour="black"),
              # legend.background = element_rect(fill = "transparent"),
              # legend.key = element_rect(fill = "transparent"),
              # legend.justification = c(1, 0), legend.position = c(0.99, 0.01),
              # plot.margin = margin(10, 100, 10, 100),
              plot.caption = element_text(size=12, colour="gray75", face="italic", hjust = 1, vjust = 1))
      g
    }else{
      sig=which(diff$significant == "yes")
      
      # again because newreps has beea manipulated above
      newreps <- reps %>%
        select(c('tracking_id','rep_name' ,'FPKM')) %>%
        spread('rep_name', 'FPKM')
      rownames(newreps) <- newreps$tracking_id
      newreps <- newreps %>% select(-c('tracking_id'))
      
      data_columns <- colnames(newreps)
      newreps[,"sum"]=apply(newreps[,data_columns], 1, sum)
      
      #Identify the genes with a grand sum FPKM of at least 5 - we will filter out the genes with very low expression across the board
      i = which(newreps[,"sum"] > 5)
      #Calculate the correlation between all pairs of data
      r=round(cor(newreps[i,data_columns], use="pairwise.complete.obs", method="pearson"),6)
      
      dist.df = reshape2::melt(as.matrix(r),varnames=c("X1","X2"))
      dist.df$value<-as.numeric(format(dist.df$value,digits=3))
      # labels = labels(obj.dists)
      g = ggplot(dist.df, aes(x=X1, y=X2, fill=value)) +
        geom_tile(color="black") +
        geom_text(aes(label=value)) +
        xlab("") + ylab("") +
        labs(title="Correlation Matrix") +
        theme(axis.text.x=element_text(angle=-90, hjust=0), axis.text.y=element_text(angle=0, hjust=1))+
        theme(panel.grid.minor=element_line(colour=NA), panel.grid.major=element_line(colour=NA), panel.background=element_rect(fill=NA, colour=NA))+
        # scale_fill_gradient2(low=heatscale[1], mid=heatscale[2], high=heatscale[3])
        scale_fill_gradientn(colours = brewer.pal(8,"YlGnBu"))+
        theme(plot.title = element_text(size=18),
              axis.title.x = element_text(size=15, vjust=0.5),
              axis.title.y = element_text(size=15, vjust=0.5),
              axis.text.x = element_text(size=12, colour="black"),
              axis.text.y = element_text(size=12, colour="black"),
              legend.text = element_text(size=15, colour="black"),
              legend.title = element_text(size=15, colour="black"),
              # legend.background = element_rect(fill = "transparent"),
              # legend.key = element_rect(fill = "transparent"),
              # legend.justification = c(1, 0), legend.position = c(0.99, 0.01),
              # plot.margin = margin(10, 100, 10, 100),
              plot.caption = element_text(size=12, colour="gray75", face="italic", hjust = 1, vjust = 1))
      g
    }
  })

  output$PCAplot <- renderPlot({
    if (is.null(input$file2$datapath) == FALSE & is.null(input$file3$datapath) == FALSE) {
      # load new data
      # DIFF <- handle_diff(input$file2$datapath)
      REPS <- handle_tracking(input$file3$datapath)
      INDEX_TABLE = data.frame(Sample = unique(REPS$rep_name))
      INDEX_TABLE$Group = INDEX_TABLE %>% separate(Sample, "Group")
      NEWREPS <- REPS %>%
        select(c('tracking_id','rep_name' ,'FPKM')) %>%
        spread('rep_name', 'FPKM')
      rownames(NEWREPS) <- NEWREPS$tracking_id
      NEWREPS <- NEWREPS %>% select(-c('tracking_id'))
      
      # plot
      fpkms<-log10(NEWREPS+1)
      pc = prcomp(fpkms, center = TRUE, scale. = TRUE)
      dat <- data.frame(obsnames=row.names(pc$x), pc$x)
      datapc <- data.frame(Sample=rownames(pc$rotation), pc$rotation)
      x="PC1"
      y="PC2"
      mult <- min(
        (max(dat[,y]) - min(dat[,y])/(max(datapc[,y])-min(datapc[,y]))),
        (max(dat[,x]) - min(dat[,x])/(max(datapc[,x])-min(datapc[,x])))
      )
      datapc <- transform(datapc,
                          v1 = .7 * mult * (get(x)),
                          v2 = .7 * mult * (get(y))
      )
      p <- ggplot(dat, aes_string(x=x, y=y)) +
        geom_point(alpha=I(1/3), size=1.2, aes(label=obsnames)) +
        geom_hline(aes(yintercept=0), size=.2) +
        geom_vline(aes(xintercept=0), size=.2) +
        geom_text(data=datapc, aes(x=v1, y=v2, label=Sample,color=Sample), vjust=1)+
        geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2,color=Sample), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75) +
        labs(title="PCA Plot") +
        theme(plot.title = element_text(size=18),
              axis.title.x = element_text(size=15, vjust=0.5),
              axis.title.y = element_text(size=15, vjust=0.5),
              axis.text.x = element_text(size=12, colour="black"),
              axis.text.y = element_text(size=12, colour="black"),
              legend.text = element_text(size=15, colour="black"),
              legend.title = element_text(size=15, colour="black"),
              # legend.background = element_rect(fill = "transparent"),
              # legend.key = element_rect(fill = "transparent"),
              # legend.justification = c(1, 0), legend.position = c(0.99, 0.01),
              # plot.margin = margin(10, 100, 10, 100),
              plot.caption = element_text(size=12, colour="gray75", face="italic", hjust = 1, vjust = 1))
      p
    }else{
      
      fpkms<-log10(newreps+1)
      pc = prcomp(fpkms, center = TRUE, scale. = TRUE)
      dat <- data.frame(obsnames=row.names(pc$x), pc$x)
      datapc <- data.frame(Sample=rownames(pc$rotation), pc$rotation)
      x="PC1"
      y="PC2"
      mult <- min(
        (max(dat[,y]) - min(dat[,y])/(max(datapc[,y])-min(datapc[,y]))),
        (max(dat[,x]) - min(dat[,x])/(max(datapc[,x])-min(datapc[,x])))
      )
      datapc <- transform(datapc,
                          v1 = .7 * mult * (get(x)),
                          v2 = .7 * mult * (get(y))
      )
      p <- ggplot(dat, aes_string(x=x, y=y)) +
        geom_point(alpha=I(1/3), size=1.2, aes(label=obsnames)) +
        geom_hline(aes(yintercept=0), size=.2) +
        geom_vline(aes(xintercept=0), size=.2) +
        geom_text(data=datapc, aes(x=v1, y=v2, label=Sample,color=Sample), vjust=1)+
        geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2,color=Sample), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75) +
        labs(title="PCA Plot") +
        theme(plot.title = element_text(size=18),
              axis.title.x = element_text(size=15, vjust=0.5),
              axis.title.y = element_text(size=15, vjust=0.5),
              axis.text.x = element_text(size=12, colour="black"),
              axis.text.y = element_text(size=12, colour="black"),
              legend.text = element_text(size=15, colour="black"),
              legend.title = element_text(size=15, colour="black"),
              # legend.background = element_rect(fill = "transparent"),
              # legend.key = element_rect(fill = "transparent"),
              # legend.justification = c(1, 0), legend.position = c(0.99, 0.01),
              # plot.margin = margin(10, 100, 10, 100),
              plot.caption = element_text(size=12, colour="gray75", face="italic", hjust = 1, vjust = 1))
      p
    }
  })


  # Table 1 diffgenes
  output$diffgenes <- DT::renderDataTable(DT::datatable({
    if (is.null(input$file2$datapath) == FALSE & is.null(input$file3$datapath) == FALSE) {
      # load new data
      DIFF <- handle_diff(input$file2$datapath)
      
      # table
      DIFF %>%
        filter(significant == "yes") %>%
        rename("Gene" = "test_id") %>%
        select(c("Gene","log2.fold_change.", "test_stat", "p_value", "q_value"))
    }else{
      diff %>%
        filter(significant == "yes") %>%
        rename("Gene" = "test_id") %>%
        select(c("Gene","log2.fold_change.", "test_stat", "p_value", "q_value"))
    }
  }))

}


# shiny
# Define UI for application 
ui <- fluidPage(
    theme = shinytheme('flatly'),
    
    # App title
    titlePanel(HTML("<h3>CuffDiff - Visualizer</h3>"), windowTitle = "CuffDiff Visualizer"),
    tags$head(
      includeCSS("style1.css")
    ),
    hr(),
    
    # # App Description
    # div(id="app_info",tags$h4("This app is designed for visualizing differential gene expression (DGE) analysis results from CuffDiff. The examplar RNAseq datasets were downloaded from GEO database, and then went through quality control (fastqc), adapter trimming (trim_galore), read alignment (Hisat2), and DGE (Cuffdiff). The required inputs to this App is `.diff` and `.tracking`.", style = "font-size: 80%")),
    # hr(),
    
    
    ########################################################## upload files
    sidebarLayout(
      # Sidebar panel for inputs ----
      sidebarPanel(
        # Input: Select a file ----
        fileInput("file1", "Upload read_groups.info File",
                  multiple = FALSE,
                  accept = c("read_groups.info", ".info")),
        tags$h5("Your sample info will be displayed on the right given `read_groups.info`."),
        # Horizontal line ----
        tags$hr(),
        # Input: Select a file ----
        fileInput("file2", "Upload gene_exp.diff File",
                  multiple = FALSE,
                  accept = c("gene_exp.diff", ".diff")),
        # # Horizontal line ----
        # tags$hr(),
        # Input: Select a file ----
        fileInput("file3", "Upload genes.read_group_tracking File",
                  multiple = FALSE,
                  accept = c("genes.read_group_tracking", ".read_group_tracking")),
        tags$h5("When gene_exp.diff and genes.read_group_tracking have not been uploaded, the figures and tables below are default examples."),
        # Horizontal line ----
        tags$hr()
      ),
      # Main panel for displaying outputs ----
      mainPanel(
        # Output: Data file ----
        # tableOutput("contents")
        # App Description
        div(id="app_info",tags$h4("This app is designed for visualizing differential gene expression (DGE) analysis results from CuffDiff. The examplar RNAseq datasets were downloaded from GEO database, and then went through quality control (fastqc), adapter trimming (trim_galore), read alignment (Hisat2), and DGE (Cuffdiff). The required inputs to this App are `gene_exp.diff` and `genes.read_group_tracking`.", style = "font-size: 100%")),
        hr(),
        DT::dataTableOutput("contents")
      )
    ),
    ########################################################## upload ended
    
    # 1. Plot01 logFPKM
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot 1"),
        hr(),
        tags$h5("A boxplot showing the distribution of log2(FPKM+1) for each sample.")
      ),
      # Output of plots
      mainPanel(
        plotlyOutput("logFPKM", height = 500)
      )
    ),
    br(),
    hr(),
    
    # 2. Plot02 logFPKM_gene
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot 2"),
        hr(),
        tags$h5("A boxplot showing the distribution of log2(FPKM+1) of a single gene for each condition. The maximum number of genes that can be displayed is 100 due to app limitation."),
        hr(),

        sliderInput("Gene1",
                    tags$h4("Index of genes:"),
                    min = 1,
                    # max = nrow(newreps),
                    max = 100,
                    value = 30,
                    step = 1)
      ),
      # Output of plots
      mainPanel(
        plotOutput("logFPKM_gene", height = 500)
      )
    ),
    br(),
    hr(),
    
    # 3. Plot03 csDensity
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot 3"),
        hr(),
        tags$h5("A density plot showing the distributions of log10(FPKM) values across conditions.")
      ),
      # Output of plots
      mainPanel(
        plotlyOutput("csDensity", height = 500)
      )
    ),
    br(),
    hr(),
    
    # 4. Plot04 scatter
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot 4"),
        hr(),
        tags$h5("A scatter plot displaying the log2(FPKM+1) values from two conditions, with significant values marked. Only first two conditions of samples can be compared due to app limitation.")
        # selectizeInput(
        #   "ConditionPicker", "Choose two conditions:",
        #   choices =unique(index_table$Group)[[1]],
        #   selected = unique(index_table$Group)[[1]][1:2],
        #   options = list(placeholder = 'Please choose two conditions', maxItems = 2),
        #   multiple = TRUE
        # )
      ),
      
      # Output of plots
      mainPanel(
        plotlyOutput("scatter", height = 500)
      )
    ),
    br(),
    hr(),
    
    # 5. Plot05 Heatmap
    
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot 5"),
        hr(),
        tags$h5("A heatmap vizualizes the expression differences among all samples.")
      ),
      # Output of plots
      mainPanel(
        plotOutput("Heatmap", height = 500)
      )
    ),
    br(),
    hr(),
    
    # 6. Plot06 volcano
    
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot 6"),
        hr(),
        tags$h5("A volcano plot comparing fold change between any two conditions and significance (-log P-values)")
      ),
      # Output of plots
      mainPanel(
        plotlyOutput("volcano", height = 500)
      )
    ),
    br(),
    hr(),
    
    
    # 7. Plot7 CorrelationMatrix
    
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot 7"),
        hr(),
        tags$h5("A correlation matrix showing the relationship / similarities among all samples and conditions.")
      ),
      # Output of plots
      mainPanel(
        plotOutput("CorrelationMatrix", height = 500)
      )
    ),
    br(),
    hr(),
    
    
    # 8. Plot8 PCAplot
    sidebarLayout(
      sidebarPanel(
        tags$h3("Plot 8"),
        hr(),
        tags$h5("A PCA plot with first two components as x and y coordinates, showing the differences between conditions by reducing demensionalities to the most important features.")
      ),
      # Output of plots
      mainPanel(
        plotOutput("PCAplot", height = 500)
      )
    ),
    br(),
    hr(),
    
    # makeTable1 - diffgenes
    sidebarLayout(
      sidebarPanel(
        tags$h3("Table 1"),
        hr(),
        tags$h5("A table summarized differential expression results. Non-significant genes are filtered out.")
      ),
      mainPanel(
        DT::dataTableOutput("diffgenes")
      )),
    br(),
    hr(),
    
    # Footer
    tags$br(),
    hr(),
    p("App created by Di Zhen in Apr 2021", HTML("&bull;"), "Find the code on Github:", tags$a(href = "https://github.com/JoKerDii/", "JoKerDii", tags$i(class = 'fa fa-github', style = 'color:#5000a5'), target = '_blank'), style = "font-size: 80%"),
    p(tags$em("Last updated: Apr 2021"), style = 'font-size:75%')
    
  )




# Run the application
shinyApp(ui = ui, server = server)


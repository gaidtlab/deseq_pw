#library(tidyverse)
#library(ggplot2)
#library(patchwork)
#library(viridis)
#library(ggpointdensity)
#library(data.table)

#library(DESeq2)
#library(apeglm)
#library(ShortRead)
#library(RUVSeq)
#library(GGally)


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ShortRead")
#BiocManager::install("DESeq2")
#BiocManager::install("apeglm")
#BiocManager::install("RUVSeq")
#install.packages("tidyverse")
#install.packages("patchwork")
#install.packages("GGally")
#install.packages("viridis")
#install.packages("ggpointdensity")

#
# ---- functions ----
#
coords_to_ucsc_links <- function(df=NULL,chr=NULL,start=NULL,end=NULL,session=NULL) {
  
  url <- paste0(session,"&position=",df[,chr],":",df[,start],"-",df[,end])
  html <- df %>% 
    mutate(link=paste0("<a href=",url,">",df[,chr],":",df[,start],"-",df[,end],"</a>")) %>% 
    relocate(link)
  return(html)
  
}



plot_counts <- function(dge,sections,region="gene") {
  
  # plot raw count sums per gene
  print("plotting bargraph ...")
  p1 <- as_tibble(rowSums(dge)) %>%
    mutate(bins=cut(value,breaks=sections)) %>%
    ggplot(aes(bins)) +
    geom_bar() +
    theme_light() +
    theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) +
    ggtitle(paste0("Binned count sum distribution per ",region,
    "\n(",region,"-wise count sums)"))
  
  print("plotting overall count sum ditribution ...")
  p2 <- as_tibble(rowSums(dge)) %>%
    ggplot(aes(x=value)) +
    scale_x_log10() +
    geom_density()+
    ggtitle(paste0("Count sum distribution per ",
    region,"\n(",region,"-wise count sums)")) +
    geom_vline(xintercept = sections,color="red") +
    theme_light()
  
  print("plotting overall count sum ditribution, by time and replicate ...")
  p3 <- as_tibble(dge) %>% 
    pivot_longer(cols=everything(),names_to = "condition",values_to = "counts") %>%
    separate(condition,into=c("label1","label2","label3"),sep="_",remove=F) %>%
    ggplot(aes(x=counts,col=label1)) +
    facet_grid(rows=vars(label3),cols=vars(label2)) +
    scale_x_log10() +
    geom_density() +
    ggtitle(paste0("Count distribution, per ",region,
    " and condition\nexcluding zero counts")) +
    theme_light()
  
  print("final plot ...")
  (p1 / p2) | p3
}

my_get_density <- function(data, mapping, N=25) {
  
  get_density <- function(x, y, n ) {
    dens <- MASS::kde2d(x = x, y = y, n = n)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  
  X <- eval_data_col(data, mapping$x)
  Y <- eval_data_col(data, mapping$y)
  
  return(get_density(x=X, y=Y, n=N))
  
}

# rm NA (results in nicer density colors)
my_filter_na <- function(df=NULL) {
  
  #print("my_filter_na ...")
  #summary(df)
  #table(is.na(df))
  
  # filter out any NA values
  df.nona <- df[!is.na(df$padj),]
  #table(is.na(df.nona))
  #summary(df.nona)
  
  # original size
  n <- dim(df)
  
  # filtered size
  n.appl <- dim(df.nona)
  
  # clean up
  res <- as.tibble(df.nona) %>% mutate(genes=rownames(df.nona))
  
  return(list(res=res,n=n,n.appl=n.appl))
  
}

my_fn <- function(data, mapping, N=25, ...){
  
  #get_density <- function(x, y, n ) {
  #  dens <- MASS::kde2d(x = x, y = y, n = n)
  #  ix <- findInterval(x, dens$x)
  #  iy <- findInterval(y, dens$y)
  #  ii <- cbind(ix, iy)
  #  return(dens$z[ii])
  #}
  
  #X <- eval_data_col(data, mapping$x)
  #Y <- eval_data_col(data, mapping$y)
  
  #data$density <- get_density(x=X, y=Y, n=N)
  data$density <- my_get_density(data,mapping,N=N)
  
  p <- ggplot(data, mapping) +
    geom_point(aes(colour=density), ...) +
    scale_color_viridis(option="H") +
    geom_abline(slope=1)
  p
}

my_volcano <- function(df=NULL,baseMean_co=0,color=F,filter_df_na=F,name="",outfile=NULL) {
 
  if(filter_df_na==T) {
    
    resl <- my_filter_na(df)
    
  } else {
    
    resl <- list(res=as_tibble(df) %>% mutate(genes=rownames(df)),
                 n=dim(df),
                 n.appl=dim(df))
  } 
  
  # prepare final dataset
  resl$res.filtered <- resl$res %>% filter(baseMean>=baseMean_co)
  baseMean_filtered <- dim(resl$res)[1] - dim(resl$res.filtered)[1]
  
  title <- paste0(name,"\nVolcano-plot, n=",dim(resl$res.filtered)[1],
                 " (minus ",resl$n-resl$n.appl," NA, minus ",
                 baseMean_filtered," baseMean),\nbaseMean=",baseMean_co)
  
  if(color==T) {
    
    p <- ggplot(data = resl$res.filtered, 
           mapping = aes(x=log2FoldChange,y=-log10(padj))) +
      geom_pointdensity(adjust=4) +
      scale_color_viridis() +
      geom_vline(xintercept =c(-1,1),color="gray") +
      geom_hline(yintercept =-log10(0.01),color="gray") + 
      geom_text_repel(data=resl$res.filtered %>% 
                        filter( ( log2FoldChange< -1 | log2FoldChange > 1) & padj<0.01),
                      aes(x=log2FoldChange,y=-log10(padj),label=genes),
                      max.overlaps=15) +
      theme_light() +
      ggtitle(title)
    
  } else {
    
    p <- ggplot(data = resl$res.filtered, 
           mapping = aes(x=log2FoldChange,y=-log10(padj))) +
      geom_point() +
      geom_vline(xintercept =c(-1,1),color="gray") +
      geom_hline(yintercept =-log10(0.01),color="gray") + 
      geom_text_repel(data=resl$res.filtered %>% 
                        filter( ( log2FoldChange< -1 | log2FoldChange > 1) & padj<0.01),
                      aes(x=log2FoldChange,y=-log10(padj),label=genes),
                      max.overlaps=15) +
      theme_light() +
      ggtitle(title)
  }
  #print(p)
  if(!is.null(outfile)) {
    ggsave(p,filename=outfile)
  } 
  return(p)
}

# df = data frame
# baseMean_co = cutoff
# padj_co = cutoff
# color = color or black and white?
# filter_df_na = NA remove?
my_ma <- function(padj_co=0.1,deobj=NULL,df=NULL,baseMean_co=0,color=F,filter_df_na=F,fc_max=0,name="",outfile=NULL) {
  
  #
  # DESeq2 MA plot
  #
  if(!is.null(deobj)) {
    #print("preparing deobj ma plot ...")
    #print("baseMean_co not applicable")
    #print("color not applicable")
    #print("filter_df_na not applicable")
    p1 <- DESeq2::plotMA(deobj,alpha=padj_co)
    if(!is.null(outfile)) {
      ggsave(p1,filename=outfile)
    }
    #p1
    return(p1)  # no other ggplot can be shown in same function
  }
  
  #
  # manual MA plot
  #
  
  if(!is.null(df)) {
  
    # filter out any NA values
    if(filter_df_na==T) {
      resl <- my_filter_na(df)
    } else {
      resl <- list(res=as.tibble(df) %>% mutate(genes=rownames(df)),
                   n=dim(df),
                   n.appl=dim(df))
    } 
    
    # (possibly) NA filtered dataset, now filter for bas mean
    resl$res.filtered <- resl$res %>% filter(baseMean>=baseMean_co)
    
    # how many were filtered out by base mean?
    baseMean_filtered <- dim(resl$res)[1] - dim(resl$res.filtered)[1]
    
    # further filter by padj
    res.sign <- resl$res.filtered %>% filter(padj<=padj_co)
    
    title <- paste0(name,"\nMA-plot, n=",dim(resl$res.filtered)[1],
                   " (minus ",resl$n-resl$n.appl," NA, minus ",
                   baseMean_filtered," baseMean),\nbaseMean=",baseMean_co,
                   ", padj=",padj_co)
    
    if(fc_max==0) {
      ylim=c(min(resl$res.filtered %>% dplyr::select(log2FoldChange) %>% pull()),
             max(resl$res.filtered %>% dplyr::select(log2FoldChange) %>% pull()))
    } else {
      ylim=c(-fc_max,fc_max)
    }
    
    if(color==T) {
    
      #print("preparing color ma plot ...")
      #print(dim(resl$res.filtered))
      p2 <- ggplot(data = resl$res.filtered, 
               mapping = aes(x=baseMean,y=log2FoldChange)) +
        geom_pointdensity(adjust=4) +
        scale_color_viridis() +
        scale_x_log10() +
        geom_smooth(se=F,color="red") +
        geom_hline(yintercept=0) +
        geom_hline(yintercept=-1,lwd=0.5,color="grey") +
        geom_hline(yintercept=1,lwd=0.5,color="grey") +
        geom_point(data = res.sign,
                   aes(x=baseMean,y=log2FoldChange),
                   color="red") +
        geom_text_repel(data=resl$res.filtered %>% 
                          filter(padj<padj_co & 
                                   baseMean>100 & 
                                   (log2FoldChange< -1 | log2FoldChange > 1)),
                        aes(x=baseMean,y=log2FoldChange,label=genes),
                        max.overlaps=15) +
        ylim(ylim) +
        theme_light() +
        ggtitle(title)
      #p2
      
    } else {
      
      #print("preparing black and white MA plot ...")
      #print(dim(resl$res.filtered))
      p2 <- ggplot(data = resl$res.filtered, 
                   mapping = aes(x=baseMean,y=log2FoldChange)) +
        geom_point() +
        scale_x_log10() +
        geom_smooth(se=F,color="red") +
        geom_hline(yintercept=0) +
        geom_hline(yintercept=-1,lwd=0.5,color="grey") +
        geom_hline(yintercept=1,lwd=0.5,color="grey") +
        geom_point(data = res.sign,
                   aes(x=baseMean,y=log2FoldChange),
                   color="red") +
        geom_text_repel(data=resl$res.filtered %>% 
                          filter(padj<padj_co & 
                                   baseMean>100 & 
                                   (log2FoldChange< -1 | log2FoldChange > 1)),
                        aes(x=baseMean,y=log2FoldChange,label=genes),
                        max.overlaps=15) +
        ylim(ylim) +
        theme_light() +
        ggtitle(title)
      
    }
    
    #print("preparing fc plot ...")
    #print(dim(resl$res.filtered))
    log2FoldChange.mean <- round(mean(resl$res.filtered %>% pull(log2FoldChange),
                                      na.rm=T), digits=2)
    p3 <- ggplot(data = resl$res.filtered, 
                 mapping = aes(x=log2FoldChange)) +
      geom_density() +
      geom_rug(data=res.sign,aes(x=log2FoldChange),color="red") +
      geom_vline(xintercept=0) +
      ylim(ylim) +
      coord_flip() +
      theme_light() +
      ggtitle(paste(name,"\nlog2FoldChange\nmean=",log2FoldChange.mean))
    #p3
  }
  
  # plt
  if(!is.null(df)) {
    #print("df")
    layout <- "
    AAAB
    "
    pf <- p2 + p3 + plot_layout(design=layout)
    if(!is.null(outfile)) {
      ggsave(pf,filename=outfile)
    }
    
    #p2 + p3 + plot_layout(design=layout)
  } else {
    print("nothing to plot ...")
    pf <- NULL
  }
  
  return(pf)
}

deseq2_RUV <- function(dds=NULL,covar=NULL,pval_co=0.1,top_co=5000,ruv_type="empirical_deseq2") {
  
  print(paste("deseq2_RUV ...",ruv_type))
  
  # make new expression set
  print("create new expression set ...")
  x <- as.factor(covar)
  set <- newSeqExpressionSet(counts(dds),
                             phenoData = data.frame(x, 
                                                    row.names=colnames(counts(dds))))
  print(pData(set))
  
  #plotRLE(set)
  #plotPCA(set)
  
  # first normalization
  print("betweenLaneNormalization (upper quantile) ...")
  set <- betweenLaneNormalization(set, which="upper")
  
  plotRLE(set)
  plotPCA(set)
  
  if(ruv_type=="empirical_deseq2") {
    
    # --- empirical_deseq2 ----
  
    # empirical control genes - DESeq2
    print("defining empirical control genes ...")
    not.sig <- rownames(results(dds))[which(results(dds)$pvalue > pval_co)]
    #length(not.sig)
    empirical.deseq2 <- rownames(set)[ rownames(set) %in% not.sig ]
    print(length(empirical.deseq2))
  
    # RUGg
    print("RUVg ...")
    set2 <- RUVg(set, empirical.deseq2, k=1)
  
    print(pData(set2))
  
    plotRLE(set2)
    plotPCA(set2)
  
    # DESeq2
    print("make DESeqDataSet ...")
    ddsruv <- DESeqDataSetFromMatrix(countData = counts(set2),
                                  colData = pData(set2),
                                  design = ~ W_1 + x)
  
    return(ddsruv)
    
  } else if(ruv_type=="empirical_edger") {
    
    # ---- empirical_edger ----
    print("running EdgeR normalizations ...")
    design <- model.matrix(~x, data=pData(set))
    y <- DGEList(counts=counts(set), group=x)
    y <- calcNormFactors(y, method="upperquartile")
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    
    print("fitting EdgeR model ...")
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef=2)
    
    print(paste("define empirical genes ... not in top",top_co))
    # Here, we consider all but the top 5,000 genes as ranked by edgeR p-values.
    top <- topTags(lrt, n=nrow(set))$table
    empirical.edger <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:top_co]))]
    print(length(empirical.edger))
    
    print("RUVg ...")
    set2 <- RUVg(set, empirical.edger, k=1)
    
    print(pData(set2))
    
    plotRLE(set2)
    plotPCA(set2)
    
    # DESeq2
    print("make DESeqDataSet ...")
    ddsruv <- DESeqDataSetFromMatrix(countData = counts(set2),
                                  colData = pData(set2),
                                  design = ~ W_1 + x)
    
    return(ddsruv)
    
  } else if(ruv_type=="factor_estimation") {
    
    # ---- factor_estimation ----
    differences <- makeGroups(x)
    set3 <- RUVs(set, rownames(counts(set)), k=1, differences)
    print(pData(set3))
    
    plotRLE(set3)
    plotPCA(set3)
    
    # DESeq2
    print("make DESeqDataSet ...")
    ddsruv <- DESeqDataSetFromMatrix(countData = counts(set3),
                                     colData = pData(set3),
                                     design = ~ W_1 + x)
    
    return(ddsruv)
    
  } else if(ruv_type=="residuals") {
    
    # ---- residuals ----
    
    print("running EdgeR normalizations ...")
    design <- model.matrix(~x, data=pData(set))
    y <- DGEList(counts=counts(set), group=x)
    y <- calcNormFactors(y, method="upperquartile")
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    
    print("fitting EdgeR model ...")
    fit <- glmFit(y, design)
    res <- residuals(fit, type="deviance")
    
    print("RUVr ...")
    set4 <- RUVr(set, rownames(counts(set)), k=1, res)
    
    print(pData(set4))
    
    plotRLE(set4)
    plotPCA(set4)
    
    # DESeq2
    print("make DESeqDataSet ...")
    ddsruv <- DESeqDataSetFromMatrix(countData = counts(set4),
                                     colData = pData(set4),
                                     design = ~ W_1 + x)
    
    return(ddsruv)
    
    
  } else {
    print("no ruv_type given")
    return(NULL)
  }
}
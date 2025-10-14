write_extended_results_table <- function(deobj=NULL,vst=NULL,res=NULL,ofile=NULL,genelength=NULL) {
  
  #print(length(genelength))
  
  raw.counts <- counts(deobj)
  raw.counts.means <- add_condition_mean_sds(counts=raw.counts,count.type="raw")
  #print(dim(raw.counts))
  
  norm.counts <- counts(deobj,normalized=T)
  norm.counts.means <- add_condition_mean_sds(counts=norm.counts,count.type="norm")
  #print(dim(norm.counts))
  
  vst.counts <- assay(vst)
  colnames(vst.counts) <- colnames(vst.counts) %>% str_remove("_[ATGC].*")
  vst.counts.means <- add_condition_mean_sds(counts=vst.counts,count.type="vst")
  #print(dim(vst.counts))
  
  tpm.counts <- convertCounts(countsMatrix = counts(deobj), unit="TPM",geneLength=genelength)
  tpm.counts.means <- add_condition_mean_sds(counts=tpm.counts,count.type="tpm")
  
  res.tab <- mutate(res,
                    raw.counts.means,
                    tpm.counts.means,
                    dplyr::select(norm.counts.means,c(contains("mean"),contains("sd"))),
                    dplyr::select(vst.counts.means,c(contains("mean"),contains("sd"))))
  
  write_delim(res.tab,file=ofile,delim = "\t")
}

# we assume the colnames and with replicate numbers (_1, _2 etc)
add_condition_mean_sds <- function(counts=NULL,count.type=NULL) {
  
  # remove replicate info
  cn <- colnames(counts) %>% 
    str_remove_all("_[123]$")

  counts.t <- as_tibble(counts)
  #condition <- unique(cn)[1]
  for(condition in unique(cn)) {
    print(condition)
    condition.subset <- counts.t[,which(cn %in% condition)]
    condition.mean <- mutate(counts.t,
                             mean=rowMeans(dplyr::select(counts.t,which(cn %in% condition)),na.rm=F),
                             sd=rowSds(x=as.matrix(condition.subset))) %>%
      dplyr::select(mean,sd)
    colnames(condition.mean) <- c(paste0(condition,"_mean"),paste0(condition,"_sd"))
    #condition.mean
    counts.t <- mutate(counts.t,condition.mean)
  }
  colnames(counts.t) <- paste0(colnames(counts.t),"_",count.type)
  return(counts.t)
}


coords_to_ucsc_links <- function(df=NULL,chr=NULL,start=NULL,end=NULL,session=NULL) {
  
  url <- paste0(session,"position=",df[,chr],":",df[,start],"-",df[,end])
  html <- df %>% 
    mutate(link=paste0("<a href=",url,">",df[,chr],":",df[,start],"-",df[,end],"</a>")) %>% 
    relocate(link)
  return(html)
  
}

coords_to_igv_links <- function(df=NULL,chr=NULL,start=NULL,end=NULL,session=NULL,pad=2000) {
  
  url <- paste0(session,"&locus=",df[,chr],":",df[,start]-pad,"-",df[,end]+pad)
  html <- df %>% 
    mutate(link=paste0("<a href=",url,">",df[,chr],":",df[,start],"-",df[,end],"</a>")) %>% 
    relocate(link)
  return(html)
  
}


prepare_dge <- function(rtype = "exon", ctype = "umicount", 
  AllCounts=NULL, meta=NULL) {
  
  dge <- as.matrix(AllCounts[[ctype]][[rtype]][["all"]])
  
  #
  # edit column names
  #
  
  # remove zUMI sample barcode
  colnames(dge) <- colnames(dge) %>% str_sub(9,100)
  
  # change to annotated column name
  colnames(dge)
  meta[match(colnames(dge),meta$barcode),]
  colnames(dge) <- meta[match(colnames(dge),meta$barcode),]$annot
  
  # reorder
  dge <- dge[,meta$annot]
  
  # remove sample barcode sequences
  colnames(dge) <- colnames(dge) %>% str_remove("_.[ATGC]*$")
  
  print("DGE:")
  print(dim(dge))
  dge[c(1:3),c(1:3)]
  
  return(dge)
}

runplots_umap <- function(rtype = "exon", ctype = "umicount", 
                          AllCounts=NULL, meta=NULL, ncomp=3) {
  
  dge <- prepare_dge(rtype,ctype,AllCounts,meta)
    
  m <- str_split(colnames(dge),"_",simplify = T)
  
  (groups <- g.celltype <- as.factor(paste(m[,1])))
  (groups <- g.time <- as.factor(paste(m[,2])))
  (groups <- g.replicate <- as.factor(paste(m[,3])))
  (groups <- g.barcode <- as.factor(paste(m[,4])))
  
  
  
  #
  # get most varying features
  #
  var <- rowVars(dge)
  counts.sel <- dge[var>quantile(var,0.99),]
  dim(counts.sel)
  print(counts.sel[1:3,1:3])
  
  # make pca
  samples.pca <- pca(t(counts.sel),ncomp=ncomp,scale=TRUE)
  
  # run umap on top pca
  counts.umap <- umap(samples.pca$variates$X,random_state=42)
  
  # reformat
  um <- counts.umap$layout
  colnames(um) <- c("dim1","dim2")
  #head(um)
  dim(um)
  
  dum <- tibble(data.frame(um),celltype=g.celltype)
  p1 <- ggplot(dum,aes(x=dim1,y=dim2,color=celltype)) +
    geom_point(show.legend=T) +
    ggtitle(paste("UMAP based on first",ncomp,"PCs")) +
    theme_classic()
  
  dum <- tibble(data.frame(um),time=g.time)
  p2 <- ggplot(dum,aes(x=dim1,y=dim2,color=time)) +
    geom_point(show.legend=T) +
    ggtitle(paste("UMAP based on first",ncomp,"PCs")) +
    theme_classic()
  
  dum <- tibble(data.frame(um),replicate=g.replicate)
  p3 <- ggplot(dum,aes(x=dim1,y=dim2,color=replicate)) +
    geom_point(show.legend=T) +
    ggtitle(paste("UMAP based on first",ncomp,"PCs")) +
    theme_classic()
  
  dum <- tibble(data.frame(um),barcode=g.barcode)
  p4 <- ggplot(dum,aes(x=dim1,y=dim2,color=barcode)) +
    geom_point(show.legend=F) +
    ggtitle(paste("UMAP based on first",ncomp,"PCs")) +
    theme_classic()
  
  print(p1 + p2 + p3 + p4 + plot_layout(nrow = 2, byrow = FALSE))
  
  
}

runplots_pca <- function(rtype = "exon", ctype = "umicount",
                         AllCounts=NULL, meta=NULL) {
  
  dge <- prepare_dge(rtype,ctype,AllCounts,meta)
  
  m <- str_split(colnames(dge),"_",simplify = T)
  
  (groups <- g.celltype <- as.factor(paste(m[,1])))
  (groups <- g.time <- as.factor(paste(m[,2])))
  (groups <- g.replicate <- as.factor(paste(m[,3])))
  (groups <- g.barcode <- as.factor(paste(m[,4])))
  
  #
  # get most varying features
  #
  var <- rowVars(dge)
  counts.sel <- dge[var>quantile(var,0.99),]
  dim(counts.sel)
  print(counts.sel[1:3,1:3])
  
  # calc pca
  samples.pca <- pca(t(counts.sel),ncomp=10,scale=TRUE)
  
  # plot pca
  p1 <- plotIndiv(samples.pca,ind.names = colnames(counts.sel),legend = TRUE,
                  group = g.celltype)
  
  # plot explained variance
  p2 <- tibble(cum_expl_var=samples.pca$cum.var,pc=c(1:length(samples.pca$cum.var))) %>% 
    ggplot(aes(x=pc,y=cum_expl_var)) + 
    geom_col() + 
    ggtitle("First 10 PCA") + 
    xlab("PC") + 
    ylab("cumulative sum of explained variance") + 
    theme_classic()
  
  # biplot
  p3 <- biplot(samples.pca)
  
  print(p2 + p3 + plot_layout(ncol = 2))

}

add_annot <- function(dge=NULL,gene_annot=NULL) {
  # we merge the annotation to the data tabl  
  # (gene_annot is a subset of dge)
  dge.t <- tibble(data.frame(dge)) %>% 
    mutate(gene_id=rownames(dge))
  dge.t <- left_join(dge.t,gene_annot,by="gene_id") %>%
    relocate(gene_id,gene_symbol,gene_type,gene_chr,gene_start,gene_end,gene_strand)
  print(dge.t[c(1:3),c(1:8)])
  return(dge.t)
}

runplots_hm <- function(rtype = "exon", ctype = "umicount", 
                     AllCounts=NULL, meta=NULL,odir=".",gene_annot=NULL) {
  
  dge <- prepare_dge(rtype,ctype,AllCounts,meta)
  
  dge.t <- add_annot(dge,gene_annot)
  
  #
  # convert to matrix for heatmap
  #
  
  dge.hm <- dge.t %>% 
    dplyr::select(c(2,8:dim(dge.t)[2])) %>%
    relocate("gene_symbol",colnames(dge.t[,c(-1,-2,-3,-4,-5,-6,-7)]))
  dge.m <- as.matrix(dge.hm[,-1],ncol = dim(dge)[2])
  rownames(dge.m) <- dge.hm$gene_symbol
  dge.m[1:3,1:3]
  str(dge.m)
  
  # plot heatmap
  plotHeatmap(dge.m,cutoff = 100000,
              pre=paste0(rtype,".",ctype),
              odir=odir)
  
}

set_HM_colors <- function(xls) {
  
  # barcodes
  col.barcode <- set_colors(xls$barcodes)
  
  # samples (zUMI barcodes)
  col.sample_id <- set_colors(xls$samples)
  col.zUMI <- set_colors(xls$zUMI.barcode)
  
  # name + NA
  col.sample_name <- set_colors(xls$descriptions)
  
  return(list(barcode=col.barcode,
              sample_id=col.sample_id,
              zUMI=col.zUMI,
              sample_name=col.sample_name))
}

set_colors <- function(x) {
  
  y <- unique(x)
  col.y <- terrain.colors(length(y))
  names(col.y) <- y
  #col.y <- as.list(col.y)
  
  return(col.y)
}

plotHeatmap_selection <- function(dge.xls.sel,pre,bin,odir) {
  
  #print("plotHeatmap_selection:")
  
  #print(pre)
  #print(dim(dge.xls.sel))
  #print(max.plot)
  
  # subsample
  sample.n <- min(dim(dge.xls.sel)[1],max.plot)
  
  mat <- dge.xls.sel[sample(c(1:dim(dge.xls.sel)[1]),sample.n),]
  
  #print("subsampled:")
  #print(dim(mat))
  #mat <- dge.xls.sel
  
  # create annotations: HCT116_24h_3_AAGCTTACCGGA
  #print(head(colnames(dge.xls.sel)))
  
  m <- str_split(colnames(dge.xls.sel),"_",simplify = T)
  
  ha <- HeatmapAnnotation(celltype=m[,1],
                          time=m[,2],
                          replicate=m[,3])
                    
  
  fn <- paste0(odir,"/",pre,".hm.",bin,".jpeg")
  if(!file.exists(fn)) {
    
    #print("creating heatmap ...")
    ht <- Heatmap(log10(mat+pseudo.count),
                  name=paste0("log10(",pre,"+",pseudo.count,")"),
                  top_annotation=ha,
                  show_column_dend = F,
                  cluster_columns = T,
                  show_row_names = T,
                  column_title = paste0("Genes with ",bin-10," to ",bin," percent of zero counts\n",
                                        dim(mat)[1]," randomly selected genes (of ",dim(dge.xls.sel)[1],")"),
                  show_column_names= F,
                  row_names_max_width = max_text_width(
                    rownames(mat), 
                    gp = gpar(fontsize = 10)))
    
    #print(fn)
    jpeg(filename=paste0(fn),
         height=60,width=60,units="cm",res=100)
    draw(ht)
    dev.off()
  } else {
    #print(paste(fn,"already exists, skipping ..."))
  }
}


plotHeatmap <- function(dge,cutoff=NULL,pre="unk",odir=".") {
  
  #print("plotHeatmap:")
  
  # summary matrix
  sums <- colSums(t(dge))
  means <- colMeans(t(dge))
  maxs <- apply(dge,1,max)
  print("Top genes with highest mean counts:")
  print(head(tibble(gene=names(sums),
                    count_sum=sums,
                    count_mean=means,
                    count_max=maxs) %>%
               dplyr::arrange(desc(count_mean))))
  
  # density of sums across all experiments
  fn <- paste0(odir,"/",pre,".count_sums.jpeg")
  if(!file.exists(fn)) {
    p1 <- tibble(sums=colSums(t(dge))) %>% 
      ggplot(aes(x=sums)) + 
      geom_density() + scale_x_log10() +
      geom_vline(xintercept=cutoff,col="gray") +
      theme_classic() +
      ggtitle(paste0(pre," gene count sums accross all samples\npreset sum cutoff=",cutoff))
    ggsave(fn)
  }
  
  # histogram of sums across all experiments
  fn <- paste0(odir,"/",pre,".count_sums.hist.jpeg")
  if(!file.exists(fn)) {
    p1 <- tibble(sums=colSums(t(dge))) %>% 
      ggplot(aes(x=sums)) + 
      geom_histogram() + scale_x_log10() +
      theme_classic() +
      ggtitle(paste0(pre," gene count sums accross all samples\npreset sum cutoff=",cutoff))
    ggsave(fn)
  }
  
  # check zero sums
  dge.zero_sum <- apply(dge,1,function(x) { return(sum(x==0)) })
  dge.zero_sum_percentage <- 100 * (dge.zero_sum / dim(dge)[2])
  
  # cumsum
  h <- hist(dge.zero_sum_percentage,breaks = 100)
  fn <- paste0(odir,"/",pre,".zero_sums.cumsum.jpeg")
  if(!file.exists(fn)) {
    p1 <- tibble(cs=cumsum(h$counts)) %>% 
      ggplot(aes(x=rank(cs),y=cs)) + 
      geom_line() + 
      theme_classic() +
      ggtitle(paste0(pre," percentage zero counts, cumulative, accross all samples\npreset sum cutoff=",cutoff))
    ggsave(fn)
  }
  
  fn <- paste0(odir,"/",pre,".percentage_zeros.jpeg")
  if(!file.exists(fn)) {
    p2 <- ggplot(as.data.frame(dge.zero_sum_percentage),aes(x=dge.zero_sum_percentage)) + 
      geom_histogram() + 
      xlab("percentage of zeros") + 
      ylab("count (transcripts)") + 
      ggtitle("Percentage of 0 counts per transcript") + 
      theme_light()
    ggsave(fn)
  }
  
  fn <- paste0(odir,"/",pre,".percentage_zeros.yzoom.jpeg")
  if(!file.exists(fn)) {
    p3 <- ggplot(as.data.frame(dge.zero_sum_percentage),aes(x=dge.zero_sum_percentage)) + 
      geom_histogram(bins=100) + 
      coord_cartesian(ylim=c(0,3000)) +
      xlab("percentage of zeros") + 
      ylab("count (transcripts)") + 
      ggtitle("Percentage of 0 counts per transcript") + 
      theme_light()
    ggsave(fn)
  }
  
  # select genes with certain % of zero counts
  for (bin in seq(from=10,to=100,by=10)) {
    print(bin)
    dge.sel <- dge.zero_sum_percentage > bin-10 & dge.zero_sum_percentage < bin
    dge.xls.sel <- dge[dge.sel==TRUE,]
    #print(str(dge.xls.sel))
    if(length(dim(dge.xls.sel)[1])>0) {
      plotHeatmap_selection(dge.xls.sel,pre,bin,odir)
    } else {
      print(paste("too few genes to plot heatmap ... bin=",bin))
    }
  }
  #draw(ht)
  #print(paste0(odir,"/",pre,"umicounts.hm.jpeg written."))
}


# rtype = exon | intron | exin
# ctype = umicount | readcount
# scols = single.samples.colnames | pooled.samples.colnames
# stype = single | pooled
# txt_fname = gn
plotSelected <- function(AllCounts,meta,rtype=NULL,ctype=NULL,scols=NULL,stype=NULL,txt_fname,
                         cutoff=NULL,odir=".") {
  
  #print("plotSelected:")
  #print(length(scols))
  dge <- as.matrix(AllCounts[[ctype]][[rtype]][["all"]])
  #print("DGE:")
  #print(dim(dge))
  dge <- dge[,scols]
  #print(dim(dge))
  #print(dge[c(1:3),c(1:3)])
  
  ofile <- paste0(odir,"/",rtype,".",ctype,".",stype,".txt")
  if(!file.exists(ofile)) {
    write.table(dge,ofile,quote=F,sep="\t")
  }
  
  dge.gn <- add_gene_symbols(dge,txt_fname)
  #print("DGE.gn:")
  #print(dim(dge.gn))
  #print(dge.gn[1:3,1:3])
  
  rl <- prepare_data(dge.gn,meta=meta)
  
  print("RL.DGE:")
  print(dim(rl$dge))
  #print(rl$dge[c(1:3),1:3])
  
  ##
  ## NEW Barcodes!
  ##
  
  # change colnames
  dge.new <- rl$dge.new
  
  if(dim(dge.new)[2]!=0) {
    colnames(dge.new) <- paste0("000000_UNK/UNK.B0.T0_XXX_UNK_",colnames(rl$dge.new))
    
    # plot heatmaps
    print("New barcodes:")
    plotHeatmap(dge.new,cutoff = cutoff,
                pre=paste0(rtype,".",ctype,".",stype,".newBarcodes"),
                odir=odir)
  }
  
  ##
  ## KNOWN barcodes!
  ##
  
  #
  # HEATMAP
  #
  print("Known barcodes:")
  plotHeatmap(rl$dge,cutoff = cutoff,
              pre=paste0(rtype,".",ctype,".",stype),
              odir=odir)
  
  
  # export table with xls barcodes only
  ofile <- paste0(odir,"/",rtype,".",ctype,".",stype,".common.txt")
  if(!file.exists(ofile)) {
    write.table(rl$dge,ofile,quote=F,sep="\t")
  }
  
  # checking for missing counts
  dge.gn.sum <- apply(dge.gn,1,sum)
  hist(log10(dge.gn.sum+1),breaks=50,
       main=paste0(rtype,".",ctype,".",stype," count sums"),
       xlab="log10(count+1)")
  
  
  dge.gn.mean <- apply(dge.gn,1,mean)
  dge.gn.median <- apply(dge.gn,1,median)
  plot(log10(dge.gn.median+1)~log10(dge.gn.mean+1),
       main=paste0(rtype,".",ctype,".",stype," median vs mean count"),
       xlab="log10(mean+1)",
       ylab="log10(median+1)")
  
  return(rl$missing)
}


# rtype = exon | intron | exin
# ctype = umicount | readcount
# scols = single.samples.colnames | pooled.samples.colnames
# stype = single | pooled
# txt_fname = gn
plotSelected_PCA_UMAP <- function(AllCounts,meta,rtype=NULL,ctype=NULL,scols=NULL,stype=NULL,txt_fname,
                                  cutoff=NULL,odir=".") {
  
  print("plotSelected_PCA_UMAP:")
  #print(length(scols))
  dge <- as.matrix(AllCounts[[ctype]][[rtype]][["all"]])
  #print("DGE:")
  #print(dim(dge))
  dge <- dge[,scols]
  #print(dim(dge))
  #print(dge[c(1:3),])
  
  dge.gn <- add_gene_symbols(dge,txt_fname)
  #print("DGE.gn:")
  #print(dim(dge.gn))
  #print(dge.gn[1:3,1:3])
  
  rl <- prepare_data(dge.gn,meta=meta)
  
  #print("RL.DGE:")
  #print(dim(rl$dge))
  #print(rl$dge[c(1:3),1:3])
  
  ##
  ## KNOWN barcodes!
  ##
  
  #
  # PCA
  #
  
  # get most varying features
  var <- rowVars(rl$dge)
  counts.sel <- rl$dge[var>quantile(var,0.99),]
  dim(counts.sel)
  print(counts.sel[1:3,1:3])
  
  # set colors
  m <- str_split(colnames(counts.sel),"_",simplify = T)
  m2 <- str_split(m[,2],"\\.",simplify=T)
  
  # (groups <- as.factor(paste(m[,2],m[,3],sep = "_"))) # 726 levels!
  # (groups <- as.factor(paste(m[,3])))  # 183 levels, to detailed
  (groups <- g.virus <- as.factor(paste(m2[,1])))  # (2) virus only
  (groups <- g.repl <- as.factor(paste(m2[,2],m2[,3],sep = "_"))) # (4) replicates biol_tech only
  (groups <- g.biol <- as.factor(paste(m2[,2],sep = "_"))) # (2) biological replicates only
  (groups <- g.tech <- as.factor(paste(m2[,3],sep = "_"))) # (2) technical replicates only
  (groups <- g.virrep <- as.factor(paste(m[,2])))  # (8) virus + replicates
  (groups <- g.samples <- as.factor(paste(m[,1])))  # (8) virus + replicates
  
  
  # calc pca
  samples.pca <- pca(t(counts.sel),ncomp=10,scale=TRUE)
  
  # plot pca
  p1 <- plotIndiv(samples.pca,ind.names = colnames(counts.sel),legend = TRUE,
                  group = g.samples)
  
  # plot explained variance
  p2 <- tibble(cum_expl_var=samples.pca$cum.var,pc=c(1:length(samples.pca$cum.var))) %>% 
    ggplot(aes(x=pc,y=cum_expl_var)) + 
    geom_col() + 
    ggtitle("First 10 PCA") + 
    xlab("PC") + 
    ylab("cumulative sum of explained variance") + 
    theme_classic()
  
  # biplot
  p3 <- biplot(samples.pca)
  
  print(p2 + p3 + plot_layout(ncol = 2))
  #print(p2)
  
  
  #
  # UMAP
  #
  
  # make pca
  samples.pca <- pca(t(counts.sel),ncomp=10,scale=TRUE)
  
  # run umap on top 10 pca
  counts.umap <- umap(samples.pca$variates$X,random_state=42)
  
  # reformat
  um <- counts.umap$layout
  colnames(um) <- c("dim1","dim2")
  #head(um)
  dim(um)
  
  dum <- tibble(data.frame(um),samples=g.samples)
  p1 <- ggplot(dum,aes(x=dim1,y=dim2,color=samples)) +
    geom_point(show.legend=T) +
    ggtitle("UMAP based on first 10 PCs") +
    theme_classic()
  print(p1)
  
  dum <- tibble(data.frame(um),samples=g.biol)
  p2 <- ggplot(dum,aes(x=dim1,y=dim2,color=samples)) +
    geom_point(show.legend=F) +
    ggtitle("UMAP based on first 10 PCs") +
    theme_classic()
  
  dum <- tibble(data.frame(um),samples=g.tech)
  p3 <- ggplot(dum,aes(x=dim1,y=dim2,color=samples)) +
    geom_point(show.legend=F) +
    ggtitle("UMAP based on first 10 PCs") +
    theme_classic()
  
  dum <- tibble(data.frame(um),samples=g.virrep)
  p4 <- ggplot(dum,aes(x=dim1,y=dim2,color=samples)) +
    geom_point(show.legend=F) +
    ggtitle("UMAP based on first 10 PCs") +
    theme_classic()
  
  
  print(p1 + p2 + p3 + p4 + plot_layout(nrow = 2, byrow = FALSE))
  
  return(rl$missing)
}



write_image_code <- function(rtype = "exon" , ctype = "umicount", stype= "single") {
  
  rmd=paste0("* ",ctype," max 10, 50 and 100% zero counts\n* ")
  
  pre=paste0(rtype,".",ctype,".",stype)
  
  for (ext in c("count_sums","percentage_zeros","percentage_zeros.yzoom")) {
    rmd=paste0(rmd,
               paste0("<a href=",pre,".",ext,".jpeg><img src=",pre,".",ext,".jpeg width=200/></a>\n"))
  }
  
  rmd=paste0(rmd,"\n* ")
  for (bin in c(10,50,100)) {
    rmd=paste0(rmd,
               paste0("<a href=",pre,".hm.",bin,".jpeg><img src=",pre,".hm.",bin,".jpeg width=200/></a>\n"))
  }
  rmd=paste0(rmd,paste0("* [original count table incl. additional barcodes](",pre,".txt)\n"))
  rmd=paste0(rmd,paste0("* [count table with only barcodes provided in xls](",pre,".common.txt)\n"))
  
  
  write(rmd,file=paste0(pre,".rmd"))
}




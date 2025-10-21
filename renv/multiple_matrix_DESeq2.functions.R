mPCA <- function(m=NULL,annot=NULL,varCutoff=0.9,
	colgroups_factor=NULL,pchgroups=NULL,pchlegend=NULL,collegend=NULL,
	title="PCA",ofile="pca.png",plotLabels=TRUE,show_biplot=T) {
  
  #print("mPCA ...")
  
  #print("make matrix ...")
  counts <- as.matrix(m)
  #print(dim(counts))
  
  # handle pch (convert to number of sign)
  #print("pchgroups set to:")
  if(!is.null(pchgroups)) {
  	#print(pchgroups)
  	pchlevels <- pchgroups
  	#pchs <- c(1:length(sort(unique(pchgroups))))
  	pchs <- c(1:length(unique(pchgroups)))
  	#names(pchs) <- sort(unique(pchgroups))
  	names(pchs) <- unique(pchgroups)
  	pchgroups <- pchs[pchgroups]
  	#print("pchgroups set to:")
  	#print(pchgroups)
  } else {
	  pchlevels=NULL
  }
  
  # variance
  print("calculate variance ...")
  rvars <- rowVars(counts)
  co <- quantile(rvars,varCutoff)
  
  print("plotting variances ...")
  
  #print("violin ...")
  p1 <- tibble(rvars) %>% ggplot(aes(y=log2(rvars),x="var")) +
    geom_violin() +
    geom_hline(yintercept = log2(co),col="red",lty=2) +
    theme_light()
  #ggsave(paste0(ofile,".p1.png"))
  
  #print("density ...")
  p2 <- tibble(rvars) %>% ggplot(aes(x=log2(rvars))) +
    geom_density() +
    coord_flip() +
    scale_y_reverse() +
    geom_vline(xintercept = log2(co),col="red",lty=2) +
    theme_light()
  #ggsave(paste0(ofile,".p2.png"))
  
  # variance versus mean value
  #print("plotting variances versus means ...")
  rmeans <- rowMeans(counts)
  
  #print("density ...")
  p4 <-  tibble(rmeans) %>% ggplot(aes(x=rmeans)) +
    geom_density() +
    geom_vline(xintercept = median(rmeans),col="red",lty=2) +
    theme_light()
  #ggsave(paste0(ofile,".p4.png"))
  
  #print("scatter ...")
  p3 <- tibble(cbind(rvars,rmeans)) %>% ggplot(aes(x=rmeans,y=log2(rvars))) +
    #geom_point() +
    #stat_density_2d_filled() +
    #scale_fill_brewer() +
    geom_hex(bins = 70) +
    scale_fill_continuous(type = "viridis") +
    geom_hline(yintercept = log2(co),col="red",lty=2) +
    geom_vline(xintercept = median(rmeans),col="red",lty=2) +
    theme_light()
  #ggsave(paste0(ofile,".p3.png"))
  
  #print("patchwork ...")
  p1 + p4 + p2 + p3 + plot_annotation(title=paste0("Distribution of variance estimates,\n",
                              "cutoff=",varCutoff,", ",round(co)))
  #p1 + p4 + p2 + p3 
  
  ggsave(paste0(ofile,".var.png"))
  
  #
  # select by variance cutoff only
  #
  counts.sel <- counts[rvars>=co,]
  print("selected features (variances only):")
  #print(dim(counts.sel)[1])
  
  # print out annot
  annot.sel <- annot[rvars>=co,]
  #print("top selected features (variances only):")
  #print(head(annot.sel))
  write_delim(data.frame(annot.sel),file=paste0(ofile,".var.bed"),delim="\t",
              col_names = F)
  
  #####
  ##### PCA
  #####
  labels <- colnames(counts.sel) %>% 
    str_replace("\\.","___") %>% 
    str_extract("^.*___") %>% 
    str_remove("___$") %>% 
    str_remove("\'")
  if(sum(is.na(labels))>0) {
    print("taking default labels ....")
    labels <- colnames(counts.sel)
  }
  #print("labels:")
  #print(labels)
  
  #print("pca ...")
  samples.pca <- pca(t(counts.sel),ncomp=3,scale=TRUE)
  
  #print("plotIndiv ...")
  #print("colgroups_factor:")
  #print(colgroups_factor)
  #print("pchgroup:")
  #print(pchgroups)
  #print("pchlevels:")
  #print(pchlevels)
  plotIndiv(samples.pca,
  	group=colgroups_factor,
	  pch=pchgroups,
	  pch.levels=pchlevels,
	  legend.title=collegend,
	  legend.title.pch=pchlegend,
	  ind.names = labels,legend = TRUE,title = title)
  ggsave(paste0(ofile,".var.pca.png"))
  
  
  if(show_biplot==TRUE) {
	  biplot(samples.pca)
	  ggsave(paste0(ofile,".var.pca.biplot.png"))
  }
  
  #
  # by var + median
  #
  counts.sel <- counts[rvars>co & rmeans>median(rmeans),]
  print("selected features (variance and means):")
  #print(paste(co,median(rmeans)))
  #print(dim(counts.sel)[1])
  
  if(dim(counts.sel)[1]>0) {
  
    # print out annot
    annot.sel <- annot[rvars>co & rmeans>median(rmeans),]
    print("top selected features (variance and means):")
    print(head(annot.sel))
    write_delim(data.frame(annot.sel),file=paste0(ofile,".varmeans.bed"),delim="\t",
                col_names = F)
    
    # pca
    samples.pca <- pca(t(counts.sel),ncomp=3,scale=TRUE)
     plotIndiv(samples.pca,
              group=colgroups_factor,
			  legend.title=collegend,
			  pch=pchgroups,
			  pch.levels=pchlevels,
			  legend.title.pch=pchlegend,
              ind.names = labels,
              legend = TRUE,
              title = title)
    ggsave(paste0(ofile,".varmeans.pca.png"))
    
	if(show_biplot==TRUE) {
    	biplot(samples.pca)
    	ggsave(paste0(ofile,".varmeans.pca.biplot.png"))
	}
  } else {
    print("no features selected ...")
  }
  
  #
  # by no var filtering (all data)
  # (just filter out var==0)
  counts.sel <- counts
  print("selected features (no filtering, only variance == 0 filtered out):")
  #print(dim(counts.sel)[1])
  
  # print out annot
  annot.sel <- annot[rvars>0,]
  print("top selected features (no filtering, only variance == 0 filtered out):")
  #print(head(annot.sel))
  write_delim(data.frame(annot.sel),file=paste0(ofile,".bed"),delim="\t",
              col_names = F)
  
  # pca
  counts.sel <- counts[rvars>0,]
  print("selected features (variance>0):")
  #print(dim(counts.sel)[1])
  
  r <- try(samples.pca <- pca(t(counts.sel),ncomp=3,scale=TRUE))
  
  if(class(r)=="try-error") {
    print("cannot perform PCA!")
  } else {
    plotIndiv(samples.pca,
		group=colgroups_factor,
		legend.title=collegend,
		pch=pchgroups,
		pch.levels=pchlevels,
		legend.title.pch=pchlegend,
		ind.names = labels,legend = TRUE,title = title)
    ggsave(paste0(ofile,".pca.png"))
    
	if(show_biplot==TRUE) {
    	biplot(samples.pca)
    	ggsave(paste0(ofile,".pca.biplot.png"))
	}
  }
}

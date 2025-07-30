## Helper functions for scRNA-seq pseudobulk edgeR DE analysis

fit_model <- function(pb,colDat,formula,geneMap,groupID='group',MDS_groups = c('Genotype','donorID','sex'),coef=NULL){
  
  if(!grepl('%s',formula,fixed=TRUE) || grepl('%s',sub('%s','',formula,fixed=TRUE),fixed=TRUE))
    stop("Invalid formula, must contain one instance of %s")
  #Convert it now in case it is still invalid in more subtle ways
  formula = as.formula(sprintf(formula,groupID))
  
  
  # create an edgeR object with counts and grouping factor
  y <- DGEList(pb, group = colDat$Genotype)
  y$samples = cbind(y$samples,colDat[match(rownames(y$samples),colDat$donorID),!colnames(colDat) %in% colnames(y$samples)])
  
  ## Plot library sizes
  df = y$samples
  df$donorID = rownames(df)
  p = ggplot(df,aes(reorder(donorID,`lib.size`),fill=group,y=`lib.size`))+
    geom_col()+
    scale_fill_manual(values = col25)+
    theme_classic(base_size = 12)+xlab('')+
    theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
          axis.line = element_blank(),axis.text = element_text(color='black'),
          panel.border = element_rect(fill=F,colour = 'black'))
  
  
  
  print(p)
  
  
  
  # filter out genes with low counts
  print("Dimensions before subsetting:")
  print(dim(y))
  print("")
  all(y$samples$lib.size>1e5)
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  print("Dimensions after subsetting:")
  print(dim(y))
  print("")
  
  
  
  ## Normalization 
  y <- calcNormFactors(y)
  # MDS plot
  cluster<-as.factor(colDat$Genotype[match(rownames(y$samples),colDat$donorID)])
  df = plotMDS(y, pch=16,col=col25[cluster], main="MDS") 
  df = data.frame(x = df$x,y=df$y,donorID = names(df$x))
  df = cbind(df,colDat[match(df$donorID,colDat$donorID),!colnames(colDat) %in% colnames(df)])
  for(f in MDS_groups){
    df$group = df[[f]]
    p = ggplot(df,aes(x,y,col=group))  +
      geom_point()+
      theme_classic(base_size = 10)+
      scale_color_manual(values = col25)+
      theme(axis.line = element_blank(),axis.text = element_text(color='black'),
            panel.border = element_rect(fill=F,colour = 'black'))
    print(p)
  }
  
  
  # create a design matrix: here we have multiple donors so also consider that in the design matrix
  design <- model.matrix(formula,data=y$samples)
  
  # estimate dispersion
  y <- estimateDisp(y, design = design,robust = T)
  # fit the model
  fit <- glmQLFit(y, design)
  plotQLDisp(fit)
  
  # Extract coefficients
  if(is.null(coef)){
    contrast = makeContrasts(groupT21-groupdiploid, levels=colnames(design))  
    qlf<-glmQLFTest(fit, contrast=contrast)
  }else{
    qlf<-glmQLFTest(fit, coef = coef)
  }
  
  
  # get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
  tt <- topTags(qlf, n = Inf,p.value = 0.05,)
  tt <- tt$table
  tt = annotateGenes(tt,geneMap = geneMap)
  
  ## Calculate logCPM
  y$logCPM <- cpm(y, log=TRUE, prior.count = 1)
  
  return(list("fit"=fit, "design"=design, "y"=y, 'tt'=tt))
}



plotPCA_inhouse_edgeR = function (count, mdat,intgroup = "condition", ntop = 500, returnData = FALSE,pc.x = 'PC1',pc.y='PC2') 
{
  rv <- rowVars(count)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(count[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% colnames(mdat))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(mdat[, intgroup, 
                                    drop = FALSE])
  
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }else {
    mdat[[intgroup]]
  }
  d <- data.frame(x = pca$x[, as.numeric(gsub('PC','',pc.x))], y = pca$x[, as.numeric(gsub('PC','',pc.y))], group = group, 
                  intgroup.df, name = colnames(count))
  colnames(d)[1:2] = c(pc.x,pc.y)
  
  p = ggplot(data = d, aes_string(x = pc.x, y = pc.y, color = "group")) + 
    geom_point(size = 3) + xlab(paste0(pc.x,": ", round(percentVar[as.numeric(gsub('PC','',pc.x))] * 
                                                          100), "% variance")) + ylab(paste0(pc.y,": ", round(percentVar[as.numeric(gsub('PC','',pc.y))] * 
                                                                                                                100), "% variance")) + coord_fixed()
  print(p)
  if (returnData) {
    #attr(d, "percentVar") <- percentVar[c(as.numeric(gsub('PC','',pc.x)),as.numeric(gsub('PC','',pc.y)))]
    
    pca_data = as.data.frame(pca$x)
    #pca_data$sampleID = rownames(pca_data)
    #m = match(pca_data$sampleID,rownames(mdat))
    pca_data = merge(pca_data, as.data.frame(mdat),by=0)
    colnames(pca_data)[colnames(pca_data) == 'Row.names'] = 'sampleID'
    return(list(p,pca_data,percentVar))
  } else{
    # just return the plot
    return(p)
  }
}


annotateGenes = function(geneTable,geneMap=NULL,tfMap=NULL,membrane_protein=NULL,cosmicGenes=NULL,tsg=NULL){
  if(is.null(geneMap)){
    ## Generic gene map
    geneMap = read.delim('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_45842_SB_Leuk13104278_GRCh38-2020-A/filtered_feature_bc_matrix/features.tsv.gz',header = F)
    colnames(geneMap) = c('ensID','geneSym','GEX')
    geneMap$geneSym[duplicated(geneMap$geneSym)] = paste0(geneMap$geneSym[duplicated(geneMap$geneSym)],'.1')
    geneMap$geneSym = gsub('_','-',geneMap$geneSym)
    
    
    #Define genomic coordinates
    require(GenomicFeatures)
    gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
    txdb = makeTxDbFromGFF(gtf)
    gns = genes(txdb)
    geneMap$chr = as.character(seqnames(gns[match(geneMap$ensID,gns$gene_id)]))
    
  }
  
  if(is.null(tfMap)){
    ## TF list
    # https://humantfs.ccbr.utoronto.ca/download.php
    #tfMap = read.table('~/lustre_mt22/generalResources/TF_list/Homo_sapiens_transcription_factors_gene_list.txt',sep = '\t')
    tfMap = read.csv('~/lustre_mt22/generalResources/TF_list/DatabaseExtract_v_1.01.csv',row.names = 1)
    # colnames(tfMap) = tfMap[1,]
    # tfMap = tfMap[-1,]  
  }
  
  if(is.null(membrane_protein)){
    ## membrane protein list
    membrane_protein = read.csv('~/lustre_mt22/generalResources/proteomic_databases/curated_transmembrane_cellsurface_proteins.csv')
    # Keep high confidence genes
    membrane_protein = membrane_protein[membrane_protein$confidence_score >=2,]
    membrane_protein$ensID = geneMap$ensID[match(membrane_protein$gene_name,geneMap$geneSym)]
    
  }
  
  if(is.null(cosmicGenes)){
    ## COSMIC gene list
    #cosmicGenes = read.csv('~/lustre_mt22/generalResources/cosmic_cancer_gene_list.csv')
    cosmicGenes = read.delim('~/lustre_mt22/generalResources/COSMIC_v100_202408/Cosmic_CancerGeneCensus_v100_GRCh38.tsv.gz',sep = '\t')
    
  }
  
  if(is.null(tsg)){
    # https://bioinfo.uth.edu/TSGene/download.cgi
    tsg = read.delim('~/lustre_mt22/generalResources/Human_TSGs.txt',header = T)
  }
  
  if(!'ensID' %in% colnames(geneTable)){
    geneTable$ensID = rownames(geneTable)  
  }
  
  geneTable$geneSym = geneMap$geneSym[match(geneTable$ensID,geneMap$ensID)]
  geneTable$chr = geneMap$chr[match(geneTable$ensID,geneMap$ensID)]
  
  # Mark TFs
  geneTable$isTF = geneTable$ensID %in% tfMap$`Ensembl ID`
  # Mark CellSurfaceMarkers
  geneTable$isCSM = geneTable$ensID %in% membrane_protein$ensID
  # Mark COSMIC genes
  geneTable$isCosmic = geneTable$geneSym %in% cosmicGenes$Gene.Symbol
  geneTable$cosmicTier = ifelse(geneTable$isCosmic, cosmicGenes$Tier[match(geneTable$geneSym,cosmicGenes$Gene.Symbol)],'')
  geneTable$tumourType = ifelse(geneTable$isCosmic, cosmicGenes$Tumour.Types.Somatic.[match(geneTable$geneSym,cosmicGenes$Gene.Symbol)],'')
  # Mark TSG
  geneTable$isTSG = geneTable$geneSym %in% tsg$GeneSymbol
  
  return(geneTable)
}

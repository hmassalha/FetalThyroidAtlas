library(Seurat)
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes



standard_clustering = function(srat,nPCs=75,clusteringRes=1,skipCCS=FALSE,s.genes=NULL,g2m.genes=NULL,runHarmony=F,harmonyVar = NULL,doPlot=T,...){
  require(Seurat)
  require(cowplot)
  srat = NormalizeData(srat)
  srat = FindVariableFeatures(srat)
  srat = ScaleData(srat,...)
  srat = RunPCA(srat, npcs = nPCs)
  
  if(runHarmony){
    # check the harmony variables provided
    if(is.null(harmonyVar) | !all(harmonyVar %in% colnames(srat@meta.data))){
      stop('No or not all variables for Harmony provided exist...')
    }else{
      require(harmony)
      
      srat = RunHarmony(srat,harmonyVar, plot_convergence = TRUE)
      # To directly access the new Harmony embeddings, use the Embeddings command.
      harmony_embeddings <- Embeddings(srat, 'harmony')
    
      if(doPlot){
        options(repr.plot.height = 5, repr.plot.width = 12)
        p1 <- DimPlot(object = srat, reduction = "harmony", pt.size = .1, group.by = harmonyVar[1])
        p2 <- VlnPlot(object = srat, features = "harmony_1", group.by = harmonyVar[1], pt.size = .1)
        plot_grid(p1,p2)  
      }
    }
  }
  #ElbowPlot(srat, ndims = nPCs)

  # Always cluster without Harmony first
  srat = FindNeighbors(srat, dims=seq(nPCs))
  srat = FindClusters(srat, resolution = clusteringRes)
  srat = RunUMAP(srat, dims=seq(nPCs))
  
  if(!runHarmony | is.null(harmonyVar)){doPlot=F}
  if(doPlot){
    p = list()
    for(var in harmonyVar){
      p[[var]] = DimPlot(object = srat, reduction = "umap", pt.size = .1, group.by = var) + ggtitle('Pre-Harmony')
    }
    plot_grid(plotlist = p)
  }
  
  if(runHarmony){
    srat = FindNeighbors(srat, dims=seq(nPCs), reduction = ifelse(runHarmony,"harmony",'pca'))
    srat = FindClusters(srat, resolution = clusteringRes, reduction = ifelse(runHarmony,"harmony",'pca'))
    srat = RunUMAP(srat, dims=seq(nPCs), reduction = ifelse(runHarmony,"harmony",'pca'))
    
    if(doPlot){
      pHarm = list()
      for(var in harmonyVar){
        pHarm[[var]] = DimPlot(object = srat, reduction = "umap", pt.size = .1, group.by = var) + ggtitle('Post-Harmony')
      }
      plot_grid(plotlist = pHarm)
    }
  }
  


  # CellCycle Scoring
  if(!skipCCS){
    if(is.null(s.genes)){
      s.genes <- cc.genes$s.genes  
    }
    if(is.null(g2m.genes)){
      g2m.genes <- cc.genes$g2m.genes  
    }
    
    srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes)  
  }
  
  return(srat)
}




## Most of the time, Seurat::merge() is fine. But because it keeps the union of genes between the 2 objects
# this might be problematic if the 2 objects were generated and pre-processed differently (eg. published external data or different genome versions etc.)
# also, should always match on geneID instead of gene names

merge_seurat_objects = function(srat1, srat2, keepAllGenes=F, genomeVersions = NULL){
  require(tidyverse)
  ## Check genome version used to map data from each of the 2 objects
  if(is.null(genomeVersions) | length(genomeVersions) != 2){
    stop('Please provide genome versions used for the 2 seurat objects')
  }
  
  ## Map gene names via ensID
  if(genomeVersions[1] == genomeVersions[2]){
    message(sprintf('Same genome version (%s) was used for both objects - No further geneID processing is needed',unique(genomeVersions)))
  }else if(genomeVersions[1] != genomeVersions[2]){
    message(sprintf('Mapping geneID for both objects to v38'))
  }
  
  ## use Seurat function to generate merged metadata
  srat1@meta.data$cellID = rownames(srat1@meta.data)
  srat2@meta.data$cellID = rownames(srat2@meta.data)
  merged.mdat = plyr::rbind.fill(srat1@meta.data,srat2@meta.data)
  
  if(sum(duplicated(merged.mdat$cellID)) > 0){
    rownames(merged.mdat) =  make.names(merged.mdat$cellID, unique=TRUE)  
  }else{
    rownames(merged.mdat) = merged.mdat$cellID
  }
  
  #rownames(merged.mdat) =  gsub('^X','',rownames(merged.mdat))
  #sum(duplicated(colnames(srat2)))
  
  
  # If we wish to check that this is equivalent to what Seurat::merge() does to merge metadata, here's the code for it
  # I've checked though and its good!
  #require(Seurat)
  #merged.srat = merge(srat1,srat2)
  #merged.mdat2 = merged.srat@meta.data
  #merged.mdat2$cellID = rownames(merged.mdat2)
  #indx <- sapply(merged.mdat, is.factor)
  #merged.mdat[indx] <- lapply(merged.mdat[indx], function(x) as.character(x))
  #dplyr::all_equal(merged.mdat,merged.mdat2)
  require(Seurat)
  
  if(keepAllGenes){
    # Use Seurat::merge() method
    merged.srat = merge(srat1,srat2)
    return(merged.srat)
  }
  
  ## If NOT keepAllGenes - ie. only keep genes that are common between the 2 srat objects. 
  
  ## get common genes
  # removing genes with 0 counts everywhere
  gene1 = rownames(srat1@assays$RNA@counts)
  gene2 = rownames(srat2@assays$RNA@counts)
  common_genes = intersect(gene1,gene2)
  
  ## Calculate the fraction of total counts attributed to those genes that got removed
  srat1_unique_genes = rownames(srat1@assays$RNA@counts)[!rownames(srat1@assays$RNA@counts) %in% common_genes]
  srat1_frac_count_lost = sum(srat1@assays$RNA@counts[rownames(srat1@assays$RNA@counts) %in% srat1_unique_genes,])/sum(srat1@assays$RNA@counts)  
  
  srat2_unique_genes = rownames(srat2@assays$RNA@counts)[!rownames(srat2@assays$RNA@counts) %in% common_genes]
  srat2_frac_count_lost = sum(srat2@assays$RNA@counts[rownames(srat2@assays$RNA@counts) %in% srat2_unique_genes,])/sum(srat2@assays$RNA@counts)  
  
  message(sprintf('# genes in srat1: %d \n# genes in srat2: %d \n# genes in common: %d\n\nFrac count lost in srat1: %f \nFrac count lost in srat2: %f \n',
                  n_distinct(rownames(srat1@assays$RNA@counts)),
                  n_distinct(rownames(srat2@assays$RNA@counts)),
                  n_distinct(common_genes),
                  srat1_frac_count_lost,
                  srat2_frac_count_lost))
  
  
  
  cnt.mtx1 = srat1@assays$RNA@counts[rownames(srat1@assays$RNA@counts) %in% common_genes,]
  cnt.mtx2 = srat2@assays$RNA@counts[rownames(cnt.mtx1),]
  if(nrow(cnt.mtx1) != nrow(cnt.mtx2)){
    stop('There is issue with the merge!')
  }
  mtx = cbind(cnt.mtx1,cnt.mtx2)
  srat = CreateSeuratObject(mtx,meta.data = merged.mdat)
  
  return(srat)
}
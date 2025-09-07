##----   Processing leukaemia scRNAseq datasets    -----##

outDir = "~/FetalThyroidAtlas/Data/published_scRNAseq/Zheng_etal_2025"

if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}



##---------------##
#### Libraries ####
##---------------##
library(tidyverse)
#source("R/utils/misc.R")
source("R/utils/sc_basicQC.R")



##----------------------------##
##   Set Global parameters  ####
##----------------------------##

maxMT = 30
minGenes = 300
minUMIs = 500  
maxBadFrac = 0.5
numPCs = 75
clusteringRes = 10
skipScrub = F
skipSoup = T
scrubScoreMax = 0.5
scrubPath='./scrubletScores.tsv'
scPath="./strainedCounts"
doPlot=T
verbose = T
skipIfExists=T
keepMTCells=T



##-----------------------------##
##     Zheng et al., 2025    ####
##-----------------------------##

sub_outDir = file.path(outDir,'Zheng25')
plotDir = file.path(sub_outDir,'Zheng25_')
outPath = file.path(sub_outDir,'Zheng25')
metadata = NULL
matchBy = NULL
cleanCountDir = file.path(sub_outDir,'cleanCount')



cleanSrat_fp = ifelse(keepMTCells,paste0(sub_outDir,'/Zheng25_clean_withMTCells.RDS'),paste0(sub_outDir,'/Zheng25_clean_noMTCells.RDS'))
if(file.exists(cleanSrat_fp) & skipIfExists){
  cleanSrat = readRDS(cleanSrat_fp)
}else{
  if(!dir.exists(sub_outDir)){
    message(sprintf('Creating output directory'))
    dir.create(sub_outDir,recursive = T)
  }
  
  #setwd(sub_outDir)
  
  # List location of cellranger outputs
  dataDirs = list.dirs('/lustre/scratch127/cellgen/cellgeni/tickets/tic-4050/results/HRA003726',full.names = T,recursive = T)
  dataDirs = dataDirs[grepl('/output/Gene/filtered',dataDirs)]
  
  names(dataDirs) = gsub('.*results/HRA003726/|/output.*$','',dirname(dataDirs))
  
  length(names(dataDirs))
  dataDirs=dataDirs[file.exists(dataDirs)]
  dataDirs=dataDirs[sapply(dataDirs, function(x){length(list.files(x))>0})]
  print(n_distinct(dataDirs))
  
  
  # Run basicQC
  message('\nPerforming scRNAseq QC...')
  QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,
                      clusteringRes=clusteringRes,cleanCountDir=cleanCountDir,
                      skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
                      metadata=metadata,matchBy=matchBy,scPath=scPath,outPath=outPath,skipIfExists=skipIfExists,
                      doPlot=doPlot,plotDir=plotDir,verbose=verbose)
  print(n_distinct(QC.output[[1]]$orig.ident))
}

srat = readRDS('~/FetalThyroidAtlas/Data/published_scRNAseq/Zheng_etal_2025/Zheng25/Zheng25_clean_noMTCells.RDS')
FeaturePlot(srat,c('TPO','TSHR','GLIS3','HMGA2','LRRK2','LMO3'))
srat$donorID = ifelse(srat$orig.ident %in% c('HRX617273','HRX617274', 'HRX617275'),'P1 (with LN-met)',
                      ifelse(srat$orig.ident %in% c('HRX617276', 'HRX617277', 'HRX617278'),'P2 (no met)','others'))

srat$tissue_type = ifelse(srat$orig.ident %in% c('HRX617273','HRX617276'),'para-carcinoma',
                          ifelse(srat$orig.ident %in% c('HRX617274','HRX617277'),'tumour',
                                 ifelse(srat$orig.ident %in% c('HRX617275','HRX617278'),'lymph_node','others')))
srat$cellID = rownames(srat@meta.data)
Idents(srat)

DimPlot(srat,cells.highlight = srat$cellID[srat$seurat_clusters == 13])
DimPlot(srat,group.by = 'annot',label = T,repel = T,label.box = T)

srat$annot = ifelse(srat$seurat_clusters ==13,'Tumour',
                    ifelse(srat$seurat_clusters %in% c(11,14,16,18),'Thyrocytes','others'))
mdat = cbind(srat@meta.data[,!colnames(srat@meta.data) %in% c('UMAP_1','UMAP_2')],srat@reductions$umap@cell.embeddings)
write.csv(mdat,'~/FetalThyroidAtlas/Data/published_scRNAseq/Zheng_etal_2025/Zheng25/Zheng25_clean_noMTCells_mdat.csv')                                 



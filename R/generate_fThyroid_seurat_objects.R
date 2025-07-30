## Generate fetal thyroid 2n/2n_T21 Seurat objects for R analysis
## Need to first run FetalThyroidAtlas/R/utils/h5ad_to_SeuratObject_conversion.ipynb to make mtx/barcode/feature files which are compatible with Seurat

setwd('~/FetalThyroidAtlas')

##----------------##
##   Libraries  ####
##----------------##
library(Seurat)
source('R/utils/sc_utils.R')


##-----------------------##
##   Fetal thyroid 2n  ####
##-----------------------##
fThyroid_fp = 'Data/fThyroid_2n_atlas.RDS'
fThy_fp = 'Data/fThyrocytes_2n_atlas.RDS'

if(!file.exists(fThy_fp)){
  # Create a seurat object of just thyrocytes from diploid fetuses
  fThyroid = Read10X('/nfs/team292/Thyroid_hm11_mt22/fThyroid_2n')
  fThyroid = CreateSeuratObject(fThyroid)
  # Add metadata
  mdat = read.csv('/nfs/team292/Thyroid_hm11_mt22/fThyroid_2n/fThyroid_2n_obs.csv.gz',row.names = 1)
  table(rownames(mdat) %in% rownames(fThyroid@meta.data))
  table(colnames(fThyroid@meta.data) %in% colnames(mdat))
  fThyroid@meta.data =  cbind(fThyroid@meta.data,mdat[match(rownames(fThyroid@meta.data),rownames(mdat)),!colnames(mdat) %in% colnames(fThyroid@meta.data)])
  fThyroid$cellID = rownames(fThyroid@meta.data)
  # Add UMAP_coords
  fThyroid = standard_clustering(fThyroid)
  umap = read.csv('/nfs/team292/Thyroid_hm11_mt22/fThyroid_2n/fThyroid_2n_UMAP.csv.gz')
  umap = column_to_rownames(umap,'X')
  fThyroid@reductions$umap@cell.embeddings[,'UMAP_1'] = umap$UMAP_1[match(rownames(fThyroid@meta.data),rownames(umap))]
  fThyroid@reductions$umap@cell.embeddings[,'UMAP_2'] = umap$UMAP_2[match(rownames(fThyroid@meta.data),rownames(umap))]
  DimPlot(fThyroid,group.by = 'cluster',label = T,label.box = T,repel = T) + NoLegend()
  saveRDS(fThyroid,fThyroid_fp)
  
  ## Thyrocytes only
  fThy = subset(fThyroid,subset = cluster == 'Thyrocytes')
  fThy = standard_clustering(fThy)
  saveRDS(fThy,fThy_fp)
}



##------------------------##
##   Fetal thyroid T21  ####
##------------------------##
## Fetal thyroid 2n-T21 age-matched object

fThyroid_2nT21_fp = 'Data/fThyroid_2nT21_agematched_atlas.RDS'

if(!file.exists(fThyroid_2nT21_fp)){
  # Create a seurat object of age-matched 2n-T21 all thyroid cell types
  a = Matrix::readMM('/nfs/team292/Thyroid_hm11_mt22/fThyroid_2nT21/')
  bc = read.delim('/nfs/team292/Thyroid_hm11_mt22/fThyroid_2nT21/barcodes.tsv.gz',sep = '\t',header = F)
  feature = read.delim('/nfs/team292/Thyroid_hm11_mt22/fThyroid_2nT21/features.tsv.gz',sep = '\t')
  colnames(a) = bc$V1
  rownames(a) = feature$V1
  
  fThyroid_2nT21 = CreateSeuratObject(a)
  
  # Add metadata
  mdat = read.csv('/nfs/team292/Thyroid_hm11_mt22/fThyroid_2nT21/fThyroid_2nT21_obs.csv.gz',row.names = 1)
  table(rownames(mdat) %in% rownames(fThyroid_2nT21@meta.data))
  table(colnames(fThyroid_2nT21@meta.data) %in% colnames(mdat))
  fThyroid_2nT21@meta.data =  cbind(fThyroid_2nT21@meta.data,mdat[match(rownames(fThyroid_2nT21@meta.data),rownames(mdat)),!colnames(mdat) %in% colnames(fThyroid_2nT21@meta.data)])
  fThyroid_2nT21$cellID = rownames(fThyroid_2nT21@meta.data)
  
  # Add UMAP_coords
  fThyroid_2nT21 = standard_clustering(fThyroid_2nT21)
  umap = read.csv('/nfs/team292/Thyroid_hm11_mt22/fThyroid_2nT21/fThyroid_2nT21_UMAP.csv.gz')
  umap = column_to_rownames(umap,'X')
  fThyroid_2nT21@reductions$umap@cell.embeddings[,'UMAP_1'] = umap$UMAP_1[match(rownames(fThyroid_2nT21@meta.data),rownames(umap))]
  fThyroid_2nT21@reductions$umap@cell.embeddings[,'UMAP_2'] = umap$UMAP_2[match(rownames(fThyroid_2nT21@meta.data),rownames(umap))]
  DimPlot(fThyroid_2nT21,group.by = 'cluster',label = T,label.box = T,repel = T) + NoLegend()
  saveRDS(fThyroid_2nT21,'Data/fThyroid_2nT21_agematched_atlas.RDS')
}


  
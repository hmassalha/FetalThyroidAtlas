aThy_Mosteiro23 = read.csv('/nfs/team292/Thyroid_hm11_mt22/public_normal_datasets_mtx/Mosteiro_2023/Mosteiro_2023_obs.csv.gz')
table(data_toPlot$cellID[data_toPlot$dataset == 'aThy_Mosteiro23'] %in% aThy_Mosteiro23$X)
wang_obs = read.csv('/nfs/team292/Thyroid_hm11_mt22/public_normal_datasets_mtx/Wang_2022/Wang_2022_obs.csv.gz')
table(gsub('_.*$','',data_toPlot$cellID[data_toPlot$dataset == 'aThy_Pu21']) %in% pu21_scanpy_obs$X)
gsub('_.*$','',data_toPlot$cellID[data_toPlot$dataset == 'aThy_Wang22'])[!gsub('_.*$','',data_toPlot$cellID[data_toPlot$dataset == 'aThy_Wang22']) %in% wang_obs$X]

pu21_scanpy_obs$cellID = paste0(pu21_scanpy_obs$donor_id,'_',pu21_scanpy_obs$cell_id)
data_toPlot$cellID[data_toPlot$dataset == 'aThy_Pu21'][1]

## Generate published adult thyrocytes-only objects, using scanpy processing


##--- Preprocess the published single-cell RNAseq datasets of adult thyroid tissues
##    Datasets included are: (more details in Supplementary Figure 5)
# 1. Mosteiro et al 2023 (adult normal) - processed seperately as only raw sequencing data available
# 2. Wang et al 2022 (adult normal) - processed seperately in this script, using Author's provided code to regenerate the same sc object
# 3. Pu et al 2021 (adult normal + PTC)
# 4. Hong et al 2023 (adult normal)
# 6. Lu et al 2023 (adult normal + PTC)
# 7. Peng et al 2020 (adult normal + PTC)

setwd('~/FetalThyroidAtlas/')
##----------------##
##   Libraries  ####
##----------------##
library(Seurat)
library(tidyverse)
#source("R/utils/misc.R")
source("R/utils/sc_basicQC.R")

preprocessed_with_scanpy=T

dataset = c('Wang22','Pu21','Hong23','Lu23','Mosteiro23','Han20')
output_objects = c('Wang22'='Data/published_scRNAseq/Wang_etal_2022/Wang_etal_2022_thyrocytes_scanpyProcessing.RDS',
                   'Pu21'='Data/published_scRNAseq/Pu_etal_2021/Pu_etal_2021_thyrocytes_scanpyProcessing.RDS',
                   'Hong23'='Data/published_scRNAseq/Hong_etal_2023/Hong_etal_2023_thyrocytes_scanpyProcessing.RDS',
                   'Mosteiro23'='Data/published_scRNAseq/Mosteiro_etal_2023/Mosteiro_etal_2023_thyrocytes_scanpyProcessing.RDS',
                   'Lu23' = 'Data/published_scRNAseq/Lu_etal_2023/Lu_etal_2023_thyrocytes_scanpyProcessing.RDS',
                   'Han20'='Data/published_scRNAseq/Han_etal_2020/Han_etal_2020_thyrocytes_scanpyProcessing.RDS')
dataset_toProcess = names(output_objects[!file.exists(output_objects)])
for(dataset in dataset_toProcess){
  if(dataset=='Mosteiro23'){
    dataDir = '/nfs/team292/Thyroid_hm11_mt22/public_normal_datasets_mtx/Mosteiro_2023/mtx/'
    if(!dir.exists(dataDir)){
      print('Cannot find dataDir, please check!')
    }
    
    srat = Read10X(dataDir)
    
    colnames(srat)
    mdat = read.csv(file.path(dataDir,'../Mosteiro_2023_obs.csv.gz'))
    colnames(mdat)[colnames(mdat) == 'X'] = 'cellID'
    rownames(mdat) = mdat$cellID
    table(mdat$cellID %in% colnames(srat))
    
    mdat = mdat[match(colnames(srat),mdat$cellID),]
    
    srat = CreateSeuratObject(srat,meta.data = mdat)
    srat = standard_clustering(srat)
    
    umap = read.csv(file.path(dataDir,'../Mosteiro_2023_UMAP.csv.gz'),row.names = 1)
    umap = as.matrix(umap[match(colnames(srat),mdat$cellID),])
    rownames(umap) = mdat$cellID
    colnames(umap) = c('UMAP_1','UMAP_2')
    srat@reductions$umap@cell.embeddings = umap
    
    ## Add metadata
    srat$annot = ifelse(srat$leiden_scvi == 0,'aTFC1',
                        ifelse(srat$leiden_scvi == 1,'aTFC2',
                               ifelse(srat$leiden_scvi == 2,'aTFC3',
                                      ifelse(srat$leiden_scvi == 4,'aTFC4',
                                             ifelse(srat$leiden_scvi == 5,'aTFC5',
                                                    ifelse(srat$leiden_scvi == 3,'aEndo/Mes','others'))))))
    View(srat@meta.data)
    DimPlot(srat,group.by = 'annot',label = T,label.box = T,repel = T) + NoLegend()
    
    saveRDS(srat,output_objects[dataset])
  }else if(dataset == 'Wang22'){
    ## Processing code extract from Author's provided script here:
    #  https://github.com/shenglei1988/scRNAseq-bilateral-PTC/blob/main/R%20script%20for%20Figure%201.R
    
    #input h5 data
    T1L.data <- Read10X_h5("Data/published_scRNAseq/Wang_etal_2022/GSM5743021_T1L.h5", use.names = T)
    T1R.data <- Read10X_h5("Data/published_scRNAseq/Wang_etal_2022/GSM5743022_T1R.h5", use.names = T)
    T2L.data <- Read10X_h5("Data/published_scRNAseq/Wang_etal_2022/GSM5743023_T2L.h5", use.names = T)
    T2R.data <- Read10X_h5("Data/published_scRNAseq/Wang_etal_2022/GSM5743024_T2R.h5", use.names = T)
    T3L.data <- Read10X_h5("Data/published_scRNAseq/Wang_etal_2022/GSM5743025_T3L.h5", use.names = T)
    T3R.data <- Read10X_h5("Data/published_scRNAseq/Wang_etal_2022/GSM5743026_T3R.h5", use.names = T)
    NT.data <- Read10X_h5("Data/published_scRNAseq/Wang_etal_2022/GSM5743027_NT.h5", use.names = T)
    
    # Create Seurat object and quality control
    T1L <- CreateSeuratObject(counts = T1L.data, project = "Thyroid_L1", min.cells = 3, min.features = 200)
    T1L$stim <- "LEFT1"
    T1L[["percent.mt"]] <- PercentageFeatureSet(T1L, pattern = "^MT-")
    #VlnPlot(T1L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    T2L <- CreateSeuratObject(counts = T2L.data, project = "Thyroid_L2", min.cells = 3, min.features = 200)
    T2L$stim <- "LEFT2"
    T2L[["percent.mt"]] <- PercentageFeatureSet(T2L, pattern = "^MT-")
    #VlnPlot(T2L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    
    T3L <- CreateSeuratObject(counts = T3L.data, project = "Thyroid_L3", min.cells = 3, min.features = 200)
    T3L$stim <- "LEFT3"
    T3L[["percent.mt"]] <- PercentageFeatureSet(T3L, pattern = "^MT-")
    #VlnPlot(T3L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    
    NT <- CreateSeuratObject(counts = NT.data, project = "Thyroid_L4", min.cells = 3, min.features = 200)
    NT$stim <- "Normal"
    NT[["percent.mt"]] <- PercentageFeatureSet(NT, pattern = "^MT-")
    #VlnPlot(NT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    
    T1R <- CreateSeuratObject(counts = T1R.data, project = "Thyroid_R1", min.cells = 3, min.features = 200)
    T1R$stim <- "RIGHT1"
    T1R[["percent.mt"]] <- PercentageFeatureSet(T1R, pattern = "^MT-")
    #VlnPlot(T1R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    
    T2R <- CreateSeuratObject(counts = T2R.data, project = "Thyroid_R2", min.cells = 3, min.features = 200)
    T2R$stim <- "RIGHT2"
    T2R[["percent.mt"]] <- PercentageFeatureSet(T2R, pattern = "^MT-")
    #VlnPlot(T2R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    
    T3R <- CreateSeuratObject(counts = T3R.data, project = "Thyroid_R3", min.cells = 3, min.features = 200)
    T3R$stim <- "T3R"
    T3R[["percent.mt"]] <- PercentageFeatureSet(T3R, pattern = "^MT-")
    #VlnPlot(T3R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    
    ##Perform integration using CCA method
    
    thyroid.anchors <- FindIntegrationAnchors(object.list = list(T1L, T2L, T3L, NT, T1R, T2R, T3R), dims = 1:20)
    thyroid.combined <- IntegrateData(anchorset = thyroid.anchors, dims = 1:20)
    thyroid.combined$cellID = gsub('_.*$','',rownames(thyroid.combined@meta.data))
    ## Add metadata
    thyroid.combined$mutation = ifelse(thyroid.combined$stim == 'T1L','RET_FARP1_fusion',
                                       ifelse(thyroid.combined$stim == 'NT','Normal','BRAF_V600E'))
    thyroid.combined$PTC_subtype = ifelse(thyroid.combined$stim %in% c('T2L','T2R'),'Follicular',
                                          ifelse(thyroid.combined$stim %in% c('NT'),'Normal','Classical'))
    DefaultAssay(thyroid.combined) = 'RNA'
    
    ## Subset to just normal thyrocytes from scanpy object 
    wang22_scanpy_obs = read.csv('/nfs/team292/Thyroid_hm11_mt22/public_normal_datasets_mtx/Wang_2022/Wang_2022_obs.csv.gz')
    table(thyroid.combined$cellID %in% wang22_scanpy_obs$X)
    table(wang22_scanpy_obs$X %in% thyroid.combined$cellID)
    cells = thyroid.combined$cellID[thyroid.combined$cellID %in% wang22_scanpy_obs$X &
                                      thyroid.combined$orig.ident == 'Thyroid_L4']
    wang22_thyrocytes = subset(thyroid.combined, cellID %in% cells)
    wang22_thyrocytes = standard_clustering(wang22_thyrocytes)
    
    ## Add umap from scanpy object
    umap = read.csv('/nfs/team292/Thyroid_hm11_mt22/public_normal_datasets_mtx/Wang_2022/Wang_2022_UMAP.csv.gz',row.names = 1)
    umap = umap[match(wang22_thyrocytes$cellID,rownames(umap)),]
    rownames(umap) = rownames(wang22_thyrocytes@reductions$umap@cell.embeddings)
    wang22_thyrocytes@reductions$umap@cell.embeddings = as.matrix(umap)
    
    DimPlot(wang22_thyrocytes)
    mdat = cbind(wang22_thyrocytes@meta.data,as.data.frame(wang22_thyrocytes@reductions$umap@cell.embeddings))
    write.csv(mdat,gsub('.RDS','_mdat.csv',output_objects[dataset]))
    saveRDS(thyroid.combined,output_objects[dataset])
    
  }else if(dataset == 'Pu21'){
    # As this dataset did not published annotation, please run R/fetal_thyrocytes_2n/02.0_annotate_Pu21.R
    message('Please run R/fetal_thyrocytes_2n/02.0_annotate_Pu21.R')
    pu21 = readRDS('Data/published_scRNAseq/Pu_etal_2021/Pu_etal_2021.RDS')
    pu21$donorID = gsub('\\..*_.*$','',pu21$cellID)
    pu21$tissue_type_2 = gsub('^.*\\.|_.*$','',pu21$cellID)
    pu21$cellID_2 = paste0(pu21$donorID,'_',gsub('.*_|-\\d$','',pu21$cellID))
    pu21 = subset(pu21,subset = tissue_type == 'para-tumour')
    
    ## Subset to just normal thyrocytes from scanpy object 
    pu21_scanpy_obs = read.csv('/nfs/team292/Thyroid_hm11_mt22/public_normal_datasets_mtx/Pu_2021/Pu_2021_obs.csv.gz')
    pu21_scanpy_obs$cellID = paste0(pu21_scanpy_obs$donor_id,'_',gsub('-\\d$','',pu21_scanpy_obs$cell_id))
    pu21_scanpy_obs$cellID_matched = pu21$cellID[match(pu21_scanpy_obs$cellID,pu21$cellID_2)]
    table(pu21$cellID_2 %in% pu21_scanpy_obs$cellID)
    table(pu21_scanpy_obs$cellID %in% pu21$cellID_2)
    
    pu21_thyrocytes = subset(pu21, cellID_2 %in% pu21_scanpy_obs$cellID)
    
    DimPlot(pu21_thyrocytes)
    
    ## Add umap from scanpy object
    umap = read.csv('/nfs/team292/Thyroid_hm11_mt22/public_normal_datasets_mtx/Pu_2021/Pu_2021_UMAP.csv.gz',row.names = 1)
    rownames(umap) = pu21_scanpy_obs$cellID_matched[match(rownames(umap),pu21_scanpy_obs$cell_id)]
    umap = umap[match(pu21_thyrocytes$cellID,rownames(umap)),]
    all(rownames(umap) == rownames(pu21_thyrocytes@reductions$umap@cell.embeddings))
    pu21_thyrocytes@reductions$umap@cell.embeddings = as.matrix(umap)
    
    DimPlot(pu21_thyrocytes)
    mdat = cbind(pu21_thyrocytes@meta.data,as.data.frame(pu21_thyrocytes@reductions$umap@cell.embeddings))
    write.csv(mdat,gsub('.RDS','_mdat.csv',output_objects[dataset]))
    saveRDS(pu21_thyrocytes,output_objects[dataset])
    
  }else if(dataset == 'Hong23'){
    srat_mtx = read.delim('Data/published_scRNAseq/Hong_etal_2023/GSE182416_Thyroid_normal_7samples_54726cells_raw_count.txt.gz',sep = '\t',header=T)
    rownames(srat_mtx)
    colnames(srat_mtx)[1:3]
    srat_mtx = as.matrix(srat_mtx[,-1])
    srat_mtx = Matrix(srat_mtx, sparse = TRUE)
    
    mdat = read.delim('Data/published_scRNAseq/Hong_etal_2023/GSE182416_Thyroid_normal_7samples_metadata.txt.gz',sep = '\t')
    mdat$Index = gsub('-','.',mdat$Index)
    rownames(mdat) = mdat$Index
    
    srat = CreateSeuratObject(srat_mtx,meta.data = mdat)
    srat = standard_clustering(srat)
    
    ## Add metadata
    srat$cellID = rownames(srat@meta.data)
    sum(!colnames(mdat) %in% colnames(srat@meta.data))
    #srat@meta.data = merge(srat@meta.data,mdat,by=0,all=T)
    
    
    View(srat@meta.data)
    DimPlot(srat,group.by = 'orig.ident',label = T,label.box = T,repel = T,cols = col25[-6]) + NoLegend()
    
    ## Subset to just cells from scanpy object
    hong23_mdat = read.csv('/nfs/team292/Thyroid_hm11_mt22/public_normal_datasets_mtx/Hong_2023/Hong_2023_obs.csv.gz')
    table(gsub('-','.',hong23_mdat$Index) %in% srat$cellID)
    srat = subset(srat,subset = cellID %in% gsub('-','.',hong23_mdat$Index))
    umap = read.csv('/nfs/team292/Thyroid_hm11_mt22/public_normal_datasets_mtx/Hong_2023/Hong_2023_UMAP.csv.gz',row.names = 1)
    rownames(umap) = gsub('-','.',rownames(umap))
    srat@reductions$umap@cell.embeddings = as.matrix(umap[match(srat$cellID,rownames(umap)),])
    DimPlot(srat)
    
    mdat = cbind(srat@meta.data,as.data.frame(srat@reductions$umap@cell.embeddings))
    write.csv(mdat,gsub('.RDS','_mdat.csv',output_objects[dataset]))
    saveRDS(srat,output_objects[dataset])
    
  }else if(dataset == 'Lu23'){
    srat = NULL
    for(sample in list.files('Data/published_scRNAseq/Lu_etal_2023/',pattern = 'UMI\\.txt\\.gz',full.names = T)){
      s = read.delim(sample,sep = '\t')
      s = CreateSeuratObject(s)
      if(is.null(srat)){
        srat = s
      }else{
        srat = merge_seurat_objects(srat,s,keepAllGenes = F,genomeVersions = c('v38','v38'))
      }
    }
    
    mdat = read.delim('Data/published_scRNAseq/Lu_etal_2023/GSE193581_celltype_annotation.txt.gz', sep = '\t')
    mdat = mdat[!grepl('^ATC',rownames(mdat)),]
    mdat$cellID = gsub('-','.',rownames(mdat))
    table(mdat$cellID %in% srat$cellID)
    rownames(mdat) = mdat$cellID
    table(rownames(mdat) %in% srat$cellID)
    
    srat@meta.data = cbind(srat@meta.data,mdat[match(colnames(srat),mdat$cellID),!colnames(mdat) %in% colnames(srat@meta.data)])
    srat = standard_clustering(srat)
    DimPlot(srat,group.by = 'celltype')
    
    ## Subset to just cells from scanpy object
    lu23_mdat = read.csv('/nfs/team292/Thyroid_hm11_mt22/public_normal_datasets_mtx/Lu_2023/Lu_2023_obs.csv.gz')
    table(gsub('-','.',lu23_mdat$X) %in% srat$cellID)
    srat = subset(srat,subset = cellID %in% gsub('-','.',lu23_mdat$X))
    umap = read.csv('/nfs/team292/Thyroid_hm11_mt22/public_normal_datasets_mtx/Lu_2023/Lu_2023_UMAP.csv.gz',row.names = 1)
    rownames(umap) = gsub('-','.',rownames(umap))
    srat@reductions$umap@cell.embeddings = as.matrix(umap[match(srat$cellID,rownames(umap)),])
    DimPlot(srat)
    
    mdat = cbind(srat@meta.data,as.data.frame(srat@reductions$umap@cell.embeddings))
    write.csv(mdat,gsub('.RDS','_mdat.csv',output_objects[dataset]))
    saveRDS(srat,output_objects[dataset])
    
  }else if(dataset == 'Han20'){
    # As this dataset did not published annotation, please run R/fetal_thyrocytes_2n/02.0_annotate_Peng21.R
    message('Please run R/fetal_thyrocytes_2n/02.0_annotate_Peng21.R')
  }
}
























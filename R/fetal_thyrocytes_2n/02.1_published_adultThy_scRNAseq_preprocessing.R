##--- Preprocess the published single-cell RNAseq datasets of adult thyroid tissues
##    Datasets included are: (more details in Supplementary Figure 5)
    # 1. Mosteiro et al 2023 (adult normal) - processed seperately as only raw sequencing data available
    # 2. Wang et al 2022 (adult normal) - processed seperately in this script, using Author's provided code to regenerate the same sc object
    # 3. Pu et al 2021 (adult normal + PTC)
    # 4. Hong et al 2023 (adult normal)
    # 6. Lu et al 2023 (adult normal + PTC)
    # 7. Peng et al 2020 (adult normal + PTC)
    # 8. Zheng et al 2025 (adult matched-paratumour + PTC) - 2 individuals both with BRAF mutations

##----------------##
##   Libraries  ####
##----------------##
library(Seurat)
library(tidyverse)
#source("R/utils/misc.R")
source("R/utils/sc_basicQC.R")


dataset = c('Wang22','Pu21','Hong23','Lu23','Mosteiro23','Peng21','Zheng25')
output_objects = c('Wang22'='Data/published_scRNAseq/Wang_etal_2022/Wang_etal_2022.RDS',
                   'Pu21'='Data/published_scRNAseq/Pu_etal_2021/Pu_etal_2021.RDS',
                   'Hong23'='Data/published_scRNAseq/Hong_etal_2023/Hong_etal_2023.RDS',
                   'Mosteiro23'='Data/published_scRNAseq/Mosteiro_etal_2023/Mosteiro_etal_2023.RDS',
                   'Lu23' = 'Data/published_scRNAseq/Lu_etal_2023/Lu_etal_2023.RDS',
                   'Peng21'='Data/published_scRNAseq/Peng_etal_2021/Peng_etal_2021.RDS',
                   'Zheng25'='Data/published_scRNAseq/Zheng_etal_2025/Zheng_etal_2025.RDS')
dataset_toProcess = names(output_objects[!file.exists(output_objects)])
for(dataset in dataset_toProcess){
  if(dataset == 'Zheng25'){
    message('Please run R/fetal_thyrocytes_2n/02.0_process_Zheng25.R')
  }
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
    T1L <- subset(T1L, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)
    
    T2L <- CreateSeuratObject(counts = T2L.data, project = "Thyroid_L2", min.cells = 3, min.features = 200)
    T2L$stim <- "LEFT2"
    T2L[["percent.mt"]] <- PercentageFeatureSet(T2L, pattern = "^MT-")
    #VlnPlot(T2L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    T2L <- subset(T2L, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)
    
    T3L <- CreateSeuratObject(counts = T3L.data, project = "Thyroid_L3", min.cells = 3, min.features = 200)
    T3L$stim <- "LEFT3"
    T3L[["percent.mt"]] <- PercentageFeatureSet(T3L, pattern = "^MT-")
    #VlnPlot(T3L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    T3L <- subset(T3L, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)
    
    NT <- CreateSeuratObject(counts = NT.data, project = "Thyroid_L4", min.cells = 3, min.features = 200)
    NT$stim <- "Normal"
    NT[["percent.mt"]] <- PercentageFeatureSet(NT, pattern = "^MT-")
    #VlnPlot(NT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    NT <- subset(NT, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)
    
    T1R <- CreateSeuratObject(counts = T1R.data, project = "Thyroid_R1", min.cells = 3, min.features = 200)
    T1R$stim <- "RIGHT1"
    T1R[["percent.mt"]] <- PercentageFeatureSet(T1R, pattern = "^MT-")
    #VlnPlot(T1R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    T1R <- subset(T1R, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)
    
    T2R <- CreateSeuratObject(counts = T2R.data, project = "Thyroid_R2", min.cells = 3, min.features = 200)
    T2R$stim <- "RIGHT2"
    T2R[["percent.mt"]] <- PercentageFeatureSet(T2R, pattern = "^MT-")
    #VlnPlot(T2R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    T2R <- subset(T2R, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)
    
    T3R <- CreateSeuratObject(counts = T3R.data, project = "Thyroid_R3", min.cells = 3, min.features = 200)
    T3R$stim <- "T3R"
    T3R[["percent.mt"]] <- PercentageFeatureSet(T3R, pattern = "^MT-")
    #VlnPlot(T3R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    T3R <- subset(T3R, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)
    
    
    ##Perform integration using CCA method
    
    thyroid.anchors <- FindIntegrationAnchors(object.list = list(T1L, T2L, T3L, NT, T1R, T2R, T3R), dims = 1:20)
    thyroid.combined <- IntegrateData(anchorset = thyroid.anchors, dims = 1:20)
    thyroid.combined
    
    #Perform an integrated analysis
    
    DefaultAssay(thyroid.combined) <- "integrated"
    #NormalizeData
    thyroid.combined <- NormalizeData(thyroid.combined, normalization.method = "LogNormalize", scale.factor = 10000)
    #FindVariableFeatures
    thyroid.combined <- FindVariableFeatures(thyroid.combined, selection.method = "vst", nfeatures = 2000)
    
    #to check nFeature_RNA", "nCount_RNA", "percent.mt
    #VlnPlot(thyroid.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    # Run the standard workflow for visualization and clustering
    thyroid.combined <- ScaleData(thyroid.combined, verbose = FALSE)
    thyroid.combined <- RunPCA(thyroid.combined, npcs = 30, verbose = FALSE)
    
    
    #Determine the 'dimensionality' of the dataset 
    thyroid.combined <- JackStraw(thyroid.combined, num.replicate = 100)
    thyroid.combined <- ScoreJackStraw(thyroid.combined, dims = 1:20)
    # JackStrawPlot(thyroid.combined, dims = 1:15)
    # ElbowPlot(thyroid.combined)
    
    
    # t-SNE and Clustering
    thyroid.combined <- RunUMAP(thyroid.combined, reduction = "pca", dims = 1:15)
    thyroid.combined <- FindNeighbors(thyroid.combined, reduction = "pca", dims = 1:15)
    thyroid.combined <- FindClusters(thyroid.combined, resolution = 0.8)
    thyroid.combined <- RunTSNE(thyroid.combined, reduction = "pca", dims = 1:15)
    thyroid.combined$cellID = colnames(thyroid.combined)
    thyroid.combined$sampleID = thyroid.combined$stim
    thyroid.combined@meta.data$stim <- factor(thyroid.combined@meta.data$stim, 
                                              levels = c("LEFT1", "LEFT2", "LEFT3", "Normal", "RIGHT1", "RIGHT2", "T3R"), 
                                              labels = c("T1L", "T2L", "T3L", "NT", "T1R", "T2R", "T3R"))
    
    # DimPlot(thyroid.combined, reduction = "umap", group.by = "sampleID",cols = col25)
    # DimPlot(thyroid.combined, reduction = "umap", group.by = "seurat_clusters",cols = c(col25,pal34H),label = T,repel = T,label.box = T) + NoLegend()
    
    #name the clusters
    thyroid.combined <- RenameIdents(thyroid.combined, `0` = "Follicular_cells", `2` = "Follicular_cells",`4` = "Follicular_cells", `5` = "Follicular_cells",`9` = "Follicular_cells", `7` = "Follicular_cells",`11` = "Follicular_cells",`22` = "Follicular_cells",`24` = "Follicular_cells",  
                                     `1` = "T_cells", `3` = "T_cells", `10` = "T_cells",`13` = "T_cells",`17` = "T_cells",`15` = "T_cells",
                                     `8` = "Myeloid_cells", `27` = "Myeloid_cells", 
                                     `14` = "Endothelial_cells", `26` = "Endothelial_cells", `25` = "Endothelial_cells", `31` = "Endothelial_cells", `23` = "Endothelial_cells", `30` = "Endothelial_cells", `16` = "Endothelial_cells", 
                                     `18` = "Pericyte",`21` = "Pericyte",`6` = "Pericyte",`12` = "Pericyte",`20` = "Fibroblast",
                                     `19` = "B_cells",`28` = "B_cells", `29` = "B_cells")
    thyroid.combined$celltype = as.character(Idents(thyroid.combined))
    
    # View(thyroid.combined@meta.data)
    
    
    # Visualization for Figure 1C and 1D
    #DimPlot(thyroid.combined, reduction = "umap",group.by = 'celltype',label = T,repel = T,label.box = T) + NoLegend()
    
    
    #VlnPlot for Figure 1E
    DefaultAssay(thyroid.combined) = 'RNA'
    thyroid.combined <- NormalizeData(thyroid.combined)
    thyroid.combined <- ScaleData(thyroid.combined)
    # VlnPlot(thyroid.combined, features = c("TG", "HIGD1B","CD3D", "CD68", "VWF", "CD79A","PDGFRA","TPSAB1"), stack = T, flip = T)
    
    
    # # # find markers for every cluster compared to all remaining cells, report only the positive ones
    # # thyroid.combined.markers <- FindAllMarkers(thyroid.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    # # write.csv(thyroid.combined.markers,"Supplementary_File_1_DEGs_among_8_known_cell_clusters.csv")
    # 
    # #heatmap Figure 1F
    # DoHeatmap(thyroid.combined, features = c("TG","KRT18","KRT19","TSHR","KRT7"
    #                                          ,"HIGD1B","CSRP2","CACNB2","COL25A1","RGS5"
    #                                          ,"CD3G","GNLY","PDCD1","CD8B","FOXP3"
    #                                          ,"CD86","LYZ","HLA-DRA","APOC1","S100A8"
    #                                          ,"VWF","ARL15","PLPP1","PTPRG","STC1"
    #                                          ,"CD79A","DERL3","BANK1","MS4A1","LY9"
    #                                          ,"PDGFRA","COL3A1","COL1A1","MGP","COL1A2"
    #                                          ,"KIT","CSF1","S100B"), angle = 90) + NoLegend()
    
    
    # ##Fraction of cell types in each sample for Figure 1G
    # table(thyroid.combined$stim,thyroid.combined$celltype)
    # write.csv(table(thyroid.combined$stim,thyroid.combined$celltype),"celltypefrequency.csv")
    
    ## Add metadata
    thyroid.combined$mutation = ifelse(thyroid.combined$stim == 'T1L','RET_FARP1_fusion',
                                       ifelse(thyroid.combined$stim == 'NT','Normal','BRAF_V600E'))
    thyroid.combined$PTC_subtype = ifelse(thyroid.combined$stim %in% c('T2L','T2R'),'Follicular',
                                          ifelse(thyroid.combined$stim %in% c('NT'),'Normal','Classical'))
    
    mdat = cbind(thyroid.combined@meta.data,as.data.frame(thyroid.combined@reductions$umap@cell.embeddings))
    write.csv(mdat,gsub('.RDS','_mdat.csv',output_objects[dataset]))
    saveRDS(thyroid.combined,output_objects[dataset])
  }else if(dataset == 'Pu21'){
    # As this dataset did not published annotation, please run R/fetal_thyrocytes_2n/02.0_annotate_Pu21.R
    message('Please run R/fetal_thyrocytes_2n/02.0_annotate_Pu21.R')
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
    
    saveRDS(srat,output_objects[dataset])
  }else if(dataset == 'Peng21'){
    # As this dataset did not published annotation, please run R/fetal_thyrocytes_2n/02.0_annotate_Peng21.R
    message('Please run R/fetal_thyrocytes_2n/02.0_annotate_Peng21.R')
  }
}
























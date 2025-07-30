## Preprocessing and annotation of published dataset Peng et al, 2021

outDir = 'Data/published_scRNAseq/Peng_etal_2021/'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

##----------------##
##   Libraries  ####
##----------------##
library(Seurat)
library(tidyverse)
source("R/utils/sc_basicQC.R")
source("R/utils/logisticRegression.R")
source("R/utils/runLR.R")
source("R/utils/misc.R")

out_fp = file.path(outDir,'Peng_etal_2021.RDS')


##---------------------------##
##   1. Preprocessing      ####
##---------------------------##


# Get sample metdata
gse <- GEOquery::getGEO("GSE158291", GSEMatrix = TRUE)
# If multiple platforms, choose the first one
gse <- gse[[1]]
sample_metadata <- pData(gse)

# --- Define your files and sample names
files <- list.files('Data/published_scRNAseq/Peng_etal_2021',pattern = '.st.txt.gz$',full.names = T)
sample_ids <- gsub('_.*$','',basename(files))

# --- Function to process a single file into sparse matrix
process_rhapsody_file <- function(file, sample_id) {
  df <- read_tsv(file, comment = "#") %>%
    select(Cell_Index, Gene, RSEC_Adjusted_Molecules) %>%
    rename(cell = Cell_Index, gene = Gene, count = RSEC_Adjusted_Molecules) %>%
    mutate(cell = paste0(sample_id, "_", cell))  # Prefix cell IDs by sample
  df = pivot_wider(df,id_cols = gene,names_from = cell,values_from = count, values_fill = 0)
  return(df)
}

# --- Read and combine all samples
combined_df <- reduce(map2(files, sample_ids, process_rhapsody_file), left_join, by = "gene")
combined_df <- combined_df %>% mutate(across(-gene, ~replace_na(.x, 0)))
combined_df <- column_to_rownames(combined_df,'gene')

# --- Create sparse count matrix
sparse_counts <- combined_df %>%
  as.matrix() %>%
  Matrix::Matrix(sparse = TRUE)

# --- Create Seurat object
seurat_obj <- CreateSeuratObject(counts = sparse_counts)

# --- Add sample info as metadata
seurat_obj$sampleID <- str_extract(colnames(seurat_obj), "^[^_]+")
seurat_obj$sex <- ifelse(seurat_obj$sampleID == 'GSM4796528','male',
                         ifelse(seurat_obj$sampleID == 'GSM4796529','female','-'))
seurat_obj$cellID = rownames(seurat_obj@meta.data)

sample_metadata$sampleID = sample_metadata$geo_accession
seurat_obj@meta.data <- left_join(seurat_obj@meta.data,sample_metadata,by='sampleID')
rownames(seurat_obj@meta.data) = seurat_obj@meta.data$cellID

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
#quantile(seurat_obj$percent.mt[seurat_obj$nFeature_RNA > 200])
seurat_obj <- subset(seurat_obj, subset = cellID %in% seurat_obj$cellID[seurat_obj$nFeature_RNA > 200 & seurat_obj$percent.mt < 30])
seurat_obj = standard_clustering(seurat_obj)

DimPlot(seurat_obj,group.by = 'sampleID')

# Save the Seurat object
saveRDS(seurat_obj,output_objects[dataset])




##-------------------------------------------------------------##
##   2. LRv1: LR with REF = fThyroid_2n, tgtSrat = Pu21      ####
##-------------------------------------------------------------##

seed = 2397
numPCs = 75
geneFilter = c('minGeneFilter','maxGeneFilter')
annot_column = c('celltype','cluster')

geneFilter = 'minGeneFilter'
annot_column = 'cluster'


##---- Import REF dataset ------##
REF.srat = readRDS('Data/fThyroid_2n_atlas.RDS')
REF.srat$cellID = rownames(REF.srat@meta.data)
REF.srat$age = REF.srat$pcw
if(annot_column == 'cluster'){
  REF.srat$annot = REF.srat[[annot_column]]
  REF.srat$annot[REF.srat$annot == 'Thyrocytes'] = REF.srat$celltype[REF.srat$annot == 'Thyrocytes']
}

REF.srat = subset(REF.srat,subset = cellID %in% REF.srat$cellID[!grepl('_Cycling',REF.srat$celltype)])
message('1. REF.srat loaded')



##---- Import tgt dataset ------##
#tgt.srat = readRDS(out_fp)
tgt.srat = seurat_obj
##---- Configure REF.srat and tgt.srat for LRv1   ------##

## Keep genes which are commonly expressed in both REF.srat and tgt.srat
genesToKeep = rowSums(tgt.srat@assays$RNA@counts)

# Also expressed in REF.srat
genesToKeep = genesToKeep[names(genesToKeep) %in% rownames(REF.srat)]
genesToKeep = names(genesToKeep)
# Remove rubbish genes
genesToKeep = genesToKeep[!grepl('^MT-|^RPL|^RPS|^MALAT1$|^NEAT1$|^AC\\d+|^AL\\d+',genesToKeep)]

# if(geneFilter == 'maxGeneFilter'){
#   # Remove soupy genes
#   soupyGenes = read.csv('~/lustre_mt22/Thyroid/Results/5_LR_fThyREF_on_pThy/fThy2nREF_downsampled_Thyrocytes_markersExprinpThy.csv')
#   soupyGenes = soupyGenes$gene[soupyGenes$frac_cellExpr >= 0.7]
#   
#   thy_prog = c('PAX8','NKX2-1','FOXE1','HHEX','UROD','DIO2','DUOX1','DUOX2','DUOXA1','DUOXA2',
#                'SLC5A5','ANO1','SLC26A4','TPO','IYD','TG','PAX8','GLIS3','TSHR','SLC16A2','SLC16A10','EXOC4','ELMO1','VPS13C','STON2','SPG11')
#   genesToKeep = genesToKeep[!genesToKeep %in% c(soupyGenes,thy_prog)]
# }


message(sprintf('Keep %d genes for LRv1 training',length(genesToKeep)))

## Subset REF.srat and tgt.srat to keep genes of interest only
REF.srat@assays$RNA@counts = REF.srat@assays$RNA@counts[genesToKeep,]
tgt.srat@assays$RNA@counts = tgt.srat@assays$RNA@counts[genesToKeep,]



##----------------------------##
##   Run Logistic Regression  ##
##----------------------------##
skipIfExists=T
ref_annot='annot'
maxCells=4000
outDir = file.path(outDir,paste0("LRv1_fThyroid2n.REF.",annot_column,"_aThyPeng21.tgt_",geneFilter))

#outDir = "~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/LRv1_fThy2n.REF/LRv1_fThy2n.REF.celltype_pPTC.tgt_maxGeneFilter"
model_fp = file.path(outDir,paste0('fThy2n.REF.',annot_column,'_trainModel_',geneFilter,'_4kmaxcells_70perc_2505.RDS'))
#model_fp = '~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/LRv1_fThy2n.REF/fThy2n.REF.celltype_trainModel_maxGeneFilter_4kmaxcells_70perc_240403.RDS'

LR_level='both'
srat_annot='seurat_clusters'
minGeneMatch = 0.99
maxCells=4000
tissue = 'thyroid'

out_prefix = paste0('fThy2n.REF.',annot_column,'_on_Peng21.tgt_maxCells_70perc_2505_')
plot_prefix = paste0('fThy2n.REF.',annot_column,'_on_Peng21.tgt_maxCells_70perc_2505_')
plot_prefix = NULL

outputs = runLR(REF.srat,ref_annot=ref_annot,srat=tgt.srat,LR_level=LR_level,srat_annot=srat_annot,
                model_fp = model_fp,outDir=outDir,plot_prefix=plot_prefix,out_prefix=out_prefix,
                minGeneMatch=minGeneMatch,maxCells = maxCells,tissue=tissue,
                scLR_TGTtype='',scLR_REFtype='',skipIfExists = skipIfExists)

message(sprintf('3. LR completed for tissue %s',tissue))


## Do heatmap
# model = readRDS(model_fp)
# 
# # Predict similarity
# logit_mtx = predictSimilarity(model, tgt.srat@assays$RNA@counts,logits = T,minGeneMatch = 0.95)

# Do heatmap
in_mtx = outputs[[2]][['scLR_tgt']]

type = tgt.srat$seurat_clusters[match(rownames(in_mtx),rownames(tgt.srat@meta.data))]
table(is.na(type))

pdf(file.path(outDir,'Peng21_LRv1_fThyroid2n.REF_scLR.pdf'),width = 10,height = 50)
show_row_names=F
hm = similarityHeatmap(in_mtx,
                       column_order = colnames(in_mtx)[order(colnames(in_mtx))],
                       row_title_rot = 0,
                       row_title_gp = gpar(fontsize=5),#row_names_gp = gpar(fontsize=5),row_names_max_width = unit(6,'cm'),
                       column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
                       split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = T)
ht = draw(hm)

dev.off()


## Add annotation based on LR
annot = do.call(rbind,apply(in_mtx,1,function(x){
  i = which(x==max(x))
  tmp = data.frame(max_LR = max(x),celltype = colnames(in_mtx)[i])
  return(tmp)
}))

tgt.srat = seurat_obj

tgt.srat$cellID = rownames(tgt.srat@meta.data)
tgt.srat@meta.data = cbind(tgt.srat@meta.data,annot[match(tgt.srat$cellID,rownames(annot)),])
tgt.srat$celltype[tgt.srat$max_LR < 1] = 'unknown'
tgt.srat$celltype[tgt.srat$celltype == 'thy_Lumen-forming'] = 'fTFC2'
tgt.srat$celltype[tgt.srat$celltype == 'thy_TH_processing'] = 'fTFC1'

DimPlot(tgt.srat,group.by = 'celltype',cols = col25,repel = T,label = T,label.box = T)
DimPlot(tgt.srat,cells.highlight = tgt.srat$cellID[tgt.srat$celltype == 'fTFC2' & tgt.srat$max_LR < 2])
DimPlot(tgt.srat,group.by = 'seurat_clusters',label = T,label.box = T) + NoLegend()
tgt.srat$tissue_type = ifelse(tgt.srat$`tissue source:ch1` == 'goiters','goiters','PTC')
DimPlot(tgt.srat,group.by = 'sex',label = T,label.box = T,cols = col25)



## Dot plot
markers = c(
  "TSHR", "NKX2-1", "PAX8", "GLIS3", "TG","SLC26A4","IYD", "HHEX", "FOXE1", "DUOXA1", "DUOXA2", "DUOX1","DUOX2", "SLC5A5", "ZNF804B", 
  "TPO", "COL23A1", "PPARGC1A", "SLC5A8", "DIO2", "TFF3", 
  "LRRK2", "HMGA2", "LMO3", "MET", "VAV3", "FN1", "LGALS3", "SERPINA1",  
  # Endothelium
  "FLT1", "PLVAP", "MECOM", "VWF","ANO2", "SLCO2A1", "CDH5", "ECSCR", "TBX1", "FLT4", "PROX1", 
  # Mesenchyme
  "CDH11","LAMA2", "COL6A3", "PDGFRA", "FBLN1", "BICC1", "MGP", 
  # Smooth muscle cells
  "GJC1", "PDGFRB","CLMN", 
  # Immune cells
  "MS4A1", "CD3E", "CD8A", "CD14", "TPSAB1", "CPA3", "KIT", "PTPRC"
)
tgt.srat$celltype[is.na(tgt.srat$celltype)] = 'unknown'
Idents(tgt.srat) = tgt.srat$seurat_clusters
DotPlot(tgt.srat,#group.by = 'celltype',#idents = c('fTFC1:2','fTFC2:2'),
        features = markers)+RotatedAxis()+
  theme(axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size=11),
        #panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),axis.line = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position = 'bottom') + xlab('') + ylab('')


tgt.srat$annot = '-'

tgt.srat$annot[tgt.srat$seurat_clusters %in% c(10,13,18)] = 'Tumour'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(15,16)] = 'Thyrocytes'

tgt.srat$annot[tgt.srat$seurat_clusters %in% c(30,12,26,28,1,32)] = 'VECs'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(34)] = 'LECs'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(23,21,2)] = 'SMCs'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(24)] = 'Mesenchymal'

tgt.srat$annot[tgt.srat$seurat_clusters %in% c(9,7)] = 'B_cells_GermlinalCentre'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(20,11,3,0)] = 'B_cells_naive'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(17,25,31)] = 'Plasma_cells'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(35)] = 'Follicular DCs'

tgt.srat$annot[tgt.srat$seurat_clusters %in% c(14)] = 'Monocytes'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(27)] = 'Macrophage'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(33)] = 'DC1'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(19)] = 'DC2'

tgt.srat$annot[tgt.srat$seurat_clusters %in% c(8,6,4,5,22)] = 'T_cells'

tgt.srat$annot[tgt.srat$seurat_clusters %in% c(29)] = 'doublets'
tgt.srat$annot[tgt.srat$annot == 'B_cells_GermlinalCentre' & tgt.srat$celltype == 'T_cells'] = 'T_cells'

# library(SoupX)
# qm = quickMarkers(tgt.srat@assays$RNA@counts[,tgt.srat$cellID[tgt.srat$seurat_clusters %in% c(14,19,27,35,33)]],tgt.srat$seurat_clusters[tgt.srat$seurat_clusters %in% c(14,19,27,35,33)])
# DimPlot(tgt.srat,cells.highlight = tgt.srat$cellID[tgt.srat$seurat_clusters==35])
a = as.data.frame(table(tgt.srat$annot,tgt.srat$celltype))
a = a[a$Freq >0,]
a = a[as.character(a$Var1) != as.character(a$Var2),]

DimPlot(tgt.srat,cells.highlight = tgt.srat$cellID[tgt.srat$annot == 'B_cells_GermlinalCentre' & tgt.srat$celltype == 'T_cells'])

DimPlot(tgt.srat,group.by = 'annot',label = T,label.box = T,cols = col25)


## Match it to fetal thyroid annotation
tgt.srat$cluster = '?'
tgt.srat$cluster[tgt.srat$annot %in% c('B_cells_GermlinalCentre','B_cells_naive','Plasma_cells','Follicular DCs')] = 'B_cells'
tgt.srat$cluster[tgt.srat$annot %in% c('T_cells')] = 'T_cells'
tgt.srat$cluster[tgt.srat$annot %in% c('DC1','DC2','Monocytes','Macrophage')] = 'Monocytes'
tgt.srat$cluster[tgt.srat$annot %in% c('LECs','VECs','Mesenchymal','SMCs','Thyrocytes','Tumour','doublets')] = 
  tgt.srat$annot[tgt.srat$annot %in% c('LECs','VECs','Mesenchymal','SMCs','Thyrocytes','Tumour','doublets')]



mdat = cbind(tgt.srat@meta.data[,colnames(tgt.srat@meta.data) !='celltype'],as.data.frame(tgt.srat@reductions$umap@cell.embeddings))
write.csv(mdat,gsub('.RDS','_mdat.csv',out_fp))
saveRDS(tgt.srat,out_fp)


## Sub-cluster thyrocytes / tumour cluster
thy = subset(tgt.srat,subset = annot %in% c('Thyrocytes','Tumour'))
thy = standard_clustering(thy)
DimPlot(thy,group.by = 'seurat_clusters',label = T,label.box = T)
qm = quickMarkers(thy@assays$RNA@counts,thy$seurat_clusters)
FeaturePlot(thy,'TPO')
DotPlot(thy,#group.by = 'celltype',#idents = c('fTFC1:2','fTFC2:2'),
        features = markers)+RotatedAxis()+
  theme(axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size=11),
        #panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),axis.line = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.position = 'bottom') + xlab('') + ylab('')

thy$annot = '-'
thy$annot[thy$seurat_clusters %in% c(0,2)] = 'Thyrocytes'
thy$annot[thy$seurat_clusters %in% c(1,3,4,5,6,7)] = 'Tumour'
thy$annot[thy$seurat_clusters %in% c(8,9)] = 'doublets'
thy$annot[thy$seurat_clusters %in% c(12)] = 'Fibroblasts'
thy$annot[thy$seurat_clusters %in% c(10)] = 'Endothelial_cells'
thy$annot[thy$seurat_clusters %in% c(11)] = 'SMCs'
DimPlot(tgt.srat,cells.highlight = thy$cellID[thy$seurat_clusters == 12])

DimPlot(thy,group.by = 'annot',label = T,label.box = T)

# Save Thyrocytes / Tumour only object
mdat = cbind(thy@meta.data[,colnames(thy@meta.data) !='celltype'],as.data.frame(thy@reductions$umap@cell.embeddings))
write.csv(mdat,gsub('.RDS','_thyOnly_mdat.csv',out_fp))
saveRDS(thy,gsub('.RDS','_thyOnly.RDS',out_fp))


## Add new annotation back to tgt.srat
tgt.srat$annot[tgt.srat$cellID %in% thy$cellID] = thy$annot[match(tgt.srat$cellID[tgt.srat$cellID %in% thy$cellID],thy$cellID)]
mdat = cbind(tgt.srat@meta.data[,colnames(tgt.srat@meta.data) !='celltype'],as.data.frame(tgt.srat@reductions$umap@cell.embeddings))
write.csv(mdat,gsub('.RDS','_mdat.csv',out_fp))
saveRDS(tgt.srat,out_fp)


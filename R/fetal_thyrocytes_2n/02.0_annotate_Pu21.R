## Preprocessing and annotation of published dataset Pu et al, 2021

outDir = 'Data/published_scRNAseq/Pu_etal_2021/'
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

out_fp = file.path(outDir,'Pu_etal_2021.RDS')


##---------------------------##
##   1. Preprocessing      ####
##---------------------------##


# Cannot run soupx and scrublet on just count matrices
dataDir = list.dirs('Data/published_scRNAseq/Pu_etal_2021/',full.names = T,recursive = F)
names(dataDir) = basename(dataDir)
names(dataDir) = gsub('_','.',names(dataDir))
dataDir = dataDir[!names(dataDir) %in% c('.ipynb.checkpoints','cleanCounts','LRv1.fThyroid2n.REF.cluster.aThyPu21.tgt.minGeneFilter')]

tgt.srat = Read10X(dataDir)
tgt.srat = CreateSeuratObject(tgt.srat)
tgt.srat[["percent.mt"]] <- PercentageFeatureSet(tgt.srat, pattern = "^MT-")

tgt.srat = subset(tgt.srat,subset = nCount_RNA >= 500 & nFeature_RNA >= 300 & percent.mt <= 30)
tgt.srat = standard_clustering(tgt.srat)
#DimPlot(tgt.srat,group.by = 'orig.ident',cols = col25,label = T,repel = T,label.box = T)

#FeaturePlot(tgt.srat,c('TG','TPO'))

mdat = cbind(tgt.srat@meta.data,as.data.frame(tgt.srat@reductions$umap@cell.embeddings))
write.csv(mdat,gsub('.RDS','_mdat.csv',output_objects[dataset]))
saveRDS(tgt.srat,out_fp)




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
tgt.srat = readRDS(out_fp)

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
outDir = file.path(outDir,paste0("LRv1_fThyroid2n.REF.",annot_column,"_aThyPu21.tgt_",geneFilter))
  
#outDir = "~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/LRv1_fThy2n.REF/LRv1_fThy2n.REF.celltype_pPTC.tgt_maxGeneFilter"
model_fp = file.path(outDir,paste0('fThy2n.REF.',annot_column,'_trainModel_',geneFilter,'_4kmaxcells_70perc_2505.RDS'))
#model_fp = '~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/LRv1_fThy2n.REF/fThy2n.REF.celltype_trainModel_maxGeneFilter_4kmaxcells_70perc_240403.RDS'

LR_level='both'
srat_annot='seurat_clusters'
minGeneMatch = 0.99
maxCells=4000
tissue = 'thyroid'

out_prefix = paste0('fThy2n.REF.',annot_column,'_on_Pu21.tgt_maxCells_70perc_2505_')
plot_prefix = paste0('fThy2n.REF.',annot_column,'_on_Pu21.tgt_maxCells_70perc_2505_')
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

pdf(file.path(outDir,'Pu21_LRv1_fThyroid2n.REF_scLR.pdf'),width = 10,height = 50)
show_row_names=F
hm = similarityHeatmap(in_mtx,
                       column_order = colnames(in_mtx)[order(colnames(in_mtx))],
                       row_title_rot = 0,
                       row_title_gp = gpar(fontsize=5),#row_names_gp = gpar(fontsize=5),row_names_max_width = unit(6,'cm'),
                       column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
                       split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = T)
ht = draw(hm)

dev.off()


# ## Add annotation based on LR
# annot = do.call(rbind,apply(in_mtx,1,function(x){
#   i = which(x==max(x))
#   tmp = data.frame(max_LR = max(x),celltype = colnames(in_mtx)[i])
#   return(tmp)
# }))
# tgt.srat$cellID = rownames(tgt.srat@meta.data)
# tgt.srat@meta.data = cbind(tgt.srat@meta.data,annot[match(tgt.srat$cellID,rownames(annot)),])
# tgt.srat$celltype[tgt.srat$max_LR < 1] = 'unknown'
# DimPlot(tgt.srat,group.by = 'celltype',cols = col25,repel = T,label = T,label.box = T)
# tgt.srat$celltype[tgt.srat$celltype == 'thy_Lument-forming'] = 'fTFC2'
# tgt.srat$celltype[tgt.srat$celltype == 'thy_TH_processing'] = 'fTFC1'
# DimPlot(tgt.srat,cells.highlight = tgt.srat$cellID[tgt.srat$celltype == 'fTFC1' & tgt.srat$max_LR < 2])
# DimPlot(tgt.srat,group.by = 'seurat_clusters',label = T)
# 
# tgt.srat$tissue_type=ifelse(grepl('\\.P$',tgt.srat$orig.ident),'para-tumour',
#                             ifelse(grepl('\\.T$',tgt.srat$orig.ident),'tumour',
#                                    ifelse(grepl('LN$',tgt.srat$orig.ident),'lymphnode',
#                                           ifelse(grepl('SC$',tgt.srat$orig.ident),'subcutaneous_metastase','others'))))
# DimPlot(tgt.srat,group.by = 'tissue_type',label = T,label.box = T,cols = col25)
# 
# mdat = cbind(tgt.srat@meta.data,as.data.frame(tgt.srat@reductions$umap@cell.embeddings))
# write.csv(mdat,gsub('.RDS','_mdat.csv',out_fp))
# saveRDS(tgt.srat,out_fp)
# 
# 
# 
# 
# 
# 
# 
# 
# ## Sub clustering thyrocytes/Tumour only
# srat.thy = subset(tgt.srat, subset = celltype %in% c('fTFC1','fTFC2'))
# srat.thy = standard_clustering(srat.thy)
# DimPlot(srat.thy,group.by = 'orig.ident',cols = col25)
# DimPlot(srat.thy,cells.highlight = srat.thy$cellID[srat.thy$orig.ident == 'PTC9.P'])
# DimPlot(srat.thy,cells.highlight = tgt.srat$cellID[tgt.srat$seurat_clusters == 2 & tgt.srat$celltype=='fTFC1'])
# DimPlot(srat.thy,cells.highlight = tgt.srat$cellID[tgt.srat$seurat_clusters == 2 & tgt.srat$celltype=='fTFC2'])
# srat.thy$group = tgt.srat$group[match(srat.thy$cellID,srat$cellID)]
# Idents(srat.thy) = srat.thy$group
# 
# thyroid_markers = c('TSHR','PAX8','GLIS3','TG', # follicular cells markers
#                     'NKX2-1','IYD','TPO','HHEX','FOXE1','DUOXA1','DUOXA2','DUOX1','DUOX2','DIO',
#                     'SLC26A4','ANO','SLC5A5','SLC16A2','SLC16A10',
#                     'DPP6','ZNF804B','KCNQ1' # thyroid epithelial cells
# )
# 
# DotPlot(srat.thy,scale = T,idents = c('fTFC1:2','fTFC2:2'),
#         features = thyroid_markers
#         #features = genes_toPlot
# )+RotatedAxis()+
#   theme(axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
#         axis.text.y = element_text(size=11),
#         #panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),axis.line = element_blank(),
#         legend.title = element_text(size=10),
#         legend.text = element_text(size=10),
#         legend.position = 'bottom') + xlab('') + ylab('')
# 
# tgt.srat$group = ifelse(tgt.srat$celltype %in% c('fTFC1','fTFC2') & tgt.srat$seurat_clusters == 2, paste0(tgt.srat$celltype,':2'),as.character(tgt.srat$celltype))
# tgt.srat$group = as.factor(tgt.srat$group)
# Idents(tgt.srat) = tgt.srat$group
# ncell_perCT = table(tgt.srat$group)
# DotPlot(tgt.srat,scale = T,idents = names(ncell_perCT[ncell_perCT > 30 & !names(ncell_perCT) %in% c('fTFC1','fTFC2')]),
#         features = thyroid_markers
#         #features = genes_toPlot
# )+RotatedAxis()+
#   theme(axis.text.x = element_text(size=9,angle = 90,vjust = 0.5,hjust = 1),
#         axis.text.y = element_text(size=11),
#         #panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),axis.line = element_blank(),
#         legend.title = element_text(size=10),
#         legend.text = element_text(size=10),
#         legend.position = 'bottom') + xlab('') + ylab('')
# 
# 
# 
# 
# 
# ##----------------------------------------------------------------##
# ##    LR_v1; REF = fTFC1 / fTFC2 / SCPs; tgt = Pu21_normalThy   ####
# ##----------------------------------------------------------------##
# source("R/fetal_thyrocytes_2n/logisticRegression.R")
# 
# # Using same model as LR for PTC
# model_fp = '~/lustre_mt22/Thyroid/Results_v2/05_LRv1_fThy.REF/LRv1_fThy.REF_pPTC.aPTC.tgt_alpha.0.5_maxCell.4k_maxGeneFilter/LRv1_fThy.REF_trainModel_alpha.0.5_maxCell.4k_maxGeneFilter.RDS'
# model = readRDS(model_fp)
# 
# # Predict similarity
# logit_mtx = predictSimilarity(model, srat@assays$RNA@counts,logits = T)
# # pp = predictSimilarity(model, tgt.srat@assays$RNA@counts,logits = F)
# 
# # Do heatmap
# in_mtx = logit_mtx
# 
# type = srat$celltype[match(rownames(in_mtx),rownames(srat@meta.data))]
# table(is.na(type))
# 
# pdf(file.path(outDir,paste0(plot_prefix,'scLR_all.pdf')),width = 10,height = 50)
# pdf('~/lustre_mt22/Thyroid/Data/published_scRNAseq/Pu_etal_2021/LRv1_fThy.REF_alpha.0.5_maxCell.4k_maxGeneFilter_scLR.pdf',width = 10,height = 50)
# show_row_names=F
# hm = similarityHeatmap(in_mtx,
#                        column_order = colnames(in_mtx)[order(colnames(in_mtx))],
#                        row_title_rot = 0,
#                        row_title_gp = gpar(fontsize=5),#row_names_gp = gpar(fontsize=5),row_names_max_width = unit(6,'cm'),
#                        column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
#                        split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = T)
# ht = draw(hm)
# 
# dev.off()
# 
# ##--- Ellie's style plot
# use_logit = F
# mtx = in_mtx
# if(!use_logit){
#   mtx = (1+exp(-mtx))**-1
# }
# 
# mtx = as.data.frame(mtx)
# mtx$cellID = rownames(mtx)
# df = pivot_longer(mtx,cols = colnames(mtx)[colnames(mtx) != 'cellID'],names_to = 'REF_celltype',values_to = 'LRpp_score')
# df$group = srat$celltype[match(df$cellID,srat$cellID)]
# 
# if(use_logit){
#   # df$LR_score_group = ifelse(df$LRpp_score > 3, 'high',
#   #                            ifelse(df$LRpp_score > 0, 'mid2',
#   #                                   ifelse(df$LRpp_score < -3,'low','mid')))
#   
# }else{
#   df$LR_score_group = ifelse(df$LRpp_score >= 0.9, 'high',
#                              ifelse(df$LRpp_score >= 0.6, 'mid2',
#                                     ifelse(df$LRpp_score >= 0.4, 'white',
#                                            ifelse(df$LRpp_score <= 0.1,'low','mid'))))  
# }
# 
# 
# df$LR_score_group = factor(df$LR_score_group,levels = rev(c('high','mid2','white','mid','low')))
# 
# df$toKeep = ifelse(grepl('fTFC',df$group),T,F)
# data_toPlot = df[df$toKeep == T,]
# data_toPlot$REF_celltype = factor(data_toPlot$REF_celltype,levels = c('fTFC1','fTFC2','SCPs','B_cells'))
# 
# # source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')
# # plotDir='~/lustre_mt22/Thyroid/Figures'
# 
# plotFun_LRv1 = function(noFrame=FALSE,noPlot=FALSE){
#   p1=ggplot(data_toPlot,aes(y=group,fill = LR_score_group))+
#     geom_bar(position = 'fill',col='black',linewidth=0.3,width=0.7)+
#     facet_grid(.~(REF_celltype),scales = 'free_y',space = 'free')+
#     #scale_fill_manual(values = rev(c('#EA2335','#F06571','white','#2D3F90'))) +
#     scale_fill_manual(values = rev(c('#EA2335','#F7B2B8','white','#B3B9DC','#2D3F90'))) +
#     scale_x_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1))+
#     theme_classic(base_size = 13) + xlab('Fraction of cells') + ylab('')+
#     theme(#panel.border = element_rect(fill=F,color='black',linewidth = 1),
#       panel.border = element_blank(),
#       axis.line.y = element_blank(),axis.ticks.y = element_blank(),
#       strip.background = element_blank(),axis.text = element_text(colour = 'black'))
#   
#   print(p1)
# }
# 
# 
# DimPlot(srat.thy,cells.highlight = srat$cellID[srat$celltype == 'fTFC2' & srat$seurat_clusters == 2])
# 
# 

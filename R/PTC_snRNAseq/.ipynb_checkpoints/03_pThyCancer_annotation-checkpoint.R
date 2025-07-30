##--- Annotate 10X single-nuclei paediatric thyroid cancer dataset ---##
##    1. LR_v1; REF = fThy_2n; tgtSrat = pPTC
##    2. Harmony to integrator Tumour and normal biopsies


outDir = "~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation"
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

setwd(outDir)


##----------------##
##   Libraries  ####
##----------------##
library(Seurat)
library(tidyverse)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")

##------------------------------------##
##   Cluster snRNASeq pPTC dataset  ####
##------------------------------------##



##---------------------------------------------------------##
##   1. LRv1: LR with REF = fThy_2n, tgtSrat = pPTC      ####
##---------------------------------------------------------##

source("~/lustre_mt22/generalScripts/utils/logisticRegression.R")
source("~/lustre_mt22/generalScripts/utils/runLR.R")
seed = 2397
numPCs = 75
geneFilter = c('minGeneFilter','maxGeneFilter')
annot_column = c('celltype','cluster')

geneFilter = 'maxGeneFilter'
annot_column = 'cluster'


##---- Import REF dataset ------##
REF.srat = readRDS('~/lustre_mt22/Thyroid/Data/fetalThyroid/fThyroid_annotated_fromHassan_jul23.RDS')
REF.srat$cellID = rownames(REF.srat@meta.data)
REF.srat$age = REF.srat$pcw
if(annot_column == 'cluster'){
  REF.srat$annot = REF.srat[[annot_column]]
  REF.srat$annot[REF.srat$annot == 'Thyrocytes'] = REF.srat$celltype[REF.srat$annot == 'Thyrocytes']
}

REF.srat = subset(REF.srat,subset = cellID %in% REF.srat$cellID[!grepl('_Cycling',REF.srat$celltype)])
message('1. REF.srat loaded')


##---- Import tgt dataset ------##
tgt.srat.fp = '~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/pPTC_clean_soupedRhoNone_may24_HARM_annotated_noDoublets.RDS'
if(file.exists(tgt.srat.fp)){
  tgt.srat = readRDS(tgt.srat.fp)
  tgt.srat$group = tgt.srat$finalAnn_broad
  tgt.srat$group[tgt.srat$group == 'Tumour'] = ifelse(tgt.srat$etiology[tgt.srat$group == 'Tumour'] == 'Met', 'Met:Y46',
                                                      paste0(tgt.srat$group[tgt.srat$group == 'Tumour'],':',tgt.srat$donor[tgt.srat$group == 'Tumour']))
  
  print(table(tgt.srat$donor))
}else{
  tgt.srat = readRDS('~/lustre_mt22/Thyroid/Results_v2/02_pThyCancer_snPreprocessing/may24/scCancerThyroid_10X_rhoLimNone__clean_noMTCells.RDS')
  
  ##----   Add metadata to the object   -----##
  tgt.srat$cellID = rownames(tgt.srat@meta.data)
  tgt.srat$donor = '?'
  tgt.srat$donor[tgt.srat$orig.ident %in% c('Thy.R13236839', 'Thy.R13236840', 'Thy.R13236841', 
                                            'Thy.R13236842', 'Thy.R13236843','Thy.R13236844')] = 'Y24'
  tgt.srat$donor[tgt.srat$orig.ident %in% c('NB14695466','NB14695467',
                                            'NB14664105','NB14664106','NB14664107','NB14664108',
                                            'NB14664109','NB14664110')] = 'Y46'
  
  tgt.srat$age = ifelse(tgt.srat$donor == 'Y24','4yrs','13yrs')
  tgt.srat$karyotype = 'diploid'
  tgt.srat$genotype = 'diploid'
  tgt.srat$organ = 'Thyroid'
  tgt.srat$sample_name = ifelse(tgt.srat$orig.ident %in% c('Thy.R13236839','Thy.R13236840','Thy.R13236841'),'Y24-TYR-0-FT-4',
                                ifelse(tgt.srat$orig.ident %in% c('Thy.R13236842', 'Thy.R13236843','Thy.R13236844'),'Y24-TYR-0-FT-1',
                                       ifelse(tgt.srat$orig.ident %in% c('NB14695466','NB14695467'),'Y46-TUM-0-FT-2d',
                                              ifelse(tgt.srat$orig.ident %in% c('NB14664105','NB14664106'),'Y46-TYR-0-FT-1',
                                                     ifelse(tgt.srat$orig.ident %in% c('NB14664107','NB14664108'),'Y46-TUM-0-FT-2',
                                                            ifelse(tgt.srat$orig.ident %in% c('NB14664109','NB14664110'),'Y46-TUM-0-FT-6',
                                                                   'others'))))))
  tgt.srat$sample_assignment = ifelse(tgt.srat$orig.ident %in% c('Thy.R13236839','Thy.R13236840','Thy.R13236841'),'Y24_4','Y24_1')
  tgt.srat$section = ifelse(tgt.srat$orig.ident %in% c('Thy.R13236839','Thy.R13236840','Thy.R13236841'),'thyroid_left_upper',
                            ifelse(tgt.srat$orig.ident %in% c('Thy.R13236842', 'Thy.R13236843','Thy.R13236844'),'thyroid_left_inferior',
                                   ifelse(tgt.srat$orig.ident %in% c('NB14695466','NB14695467','NB14664107','NB14664108'),'thyroid_tum',
                                          ifelse(tgt.srat$orig.ident %in% c('NB14664105','NB14664106'),'thyroid_normal',
                                                 ifelse(tgt.srat$orig.ident %in% c('NB14664109','NB14664110'),'lymph_met',
                                                        'others')))))
  tgt.srat$tissue_source = 'paediatric'
  tgt.srat$tissue_type = 'gland'
  tgt.srat$etiology = ifelse(tgt.srat$orig.ident %in% c('Thy.R13236839','Thy.R13236840','Thy.R13236841','NB14664105','NB14664106'),'Normal',
                             ifelse(tgt.srat$orig.ident %in% c('NB14664109','NB14664110'),'Met','Tumour'))
  tgt.srat$seq_chemistry = 'scNuc_5DUAL'
  tgt.srat$seq_aligner = 'cellranger'
  
  tgt.srat$gender = 'female'
  ## Determine sex
  # avg.expr.sexMarkers = AverageExpression(tgt.srat,group.by = 'donor',features = c('XIST','RPS4Y1'))
  # avg.expr.sexMarkers = as.data.frame(t(avg.expr.sexMarkers[[1]]))
  # avg.expr.sexMarkers$sex = tgt.srat@meta.data$sex[match(rownames(avg.expr.sexMarkers),tgt.srat@meta.data$donorID)]
  # avg.expr.sexMarkers$pred.sex = ifelse(avg.expr.sexMarkers$XIST > avg.expr.sexMarkers$RPS4Y1,'F',
  #                                       ifelse(avg.expr.sexMarkers$XIST < avg.expr.sexMarkers$RPS4Y1,'M','??'))
  
  
  
  
  tgt.srat.harm = standard_clustering(tgt.srat,runHarmony = T,harmonyVar = c('donor'))
  saveRDS(tgt.srat.harm,'~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/pPTC_clean_soupedRhoNone_may24_HARM.RDS')
  
  tgt.srat = standard_clustering(tgt.srat)
  saveRDS(tgt.srat,'~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/pPTC_clean_soupedRhoNone_may24.RDS')
}



# DimPlot(tgt.srat.harm,group.by = 'orig.ident')
# DimPlot(tgt.srat,group.by = 'orig.ident')
# DimPlot(tgt.srat,label = T,label.box = T,repel = T) + NoLegend() 
# DimPlot(tgt.srat,group.by = 'annot',label = T,label.box = T,repel = T) + NoLegend() 

message(('2. tgt.srat loaded'))





##---- Configure REF.srat and tgt.srat for LRv1   ------##

## Keep genes which are commonly expressed in both REF.srat and tgt.srat
# Genes expressed in > 50 cells in tgt.srat
genesToKeep = rowSums(tgt.srat@assays$RNA@counts)
genesToKeep = genesToKeep[genesToKeep>50]
# Also expressed in REF.srat
genesToKeep = genesToKeep[names(genesToKeep) %in% rownames(REF.srat)]
genesToKeep = names(genesToKeep)
# Remove rubbish genes
genesToKeep = genesToKeep[!grepl('^MT-|^RPL|^RPS|^MALAT1$|^NEAT1$|^AC\\d+|^AL\\d+',genesToKeep)]

if(geneFilter == 'maxGeneFilter'){
  # Remove soupy genes
  soupyGenes = read.csv('~/lustre_mt22/Thyroid/Results/5_LR_fThyREF_on_pThy/fThy2nREF_downsampled_Thyrocytes_markersExprinpThy.csv')
  soupyGenes = soupyGenes$gene[soupyGenes$frac_cellExpr >= 0.7]

  thy_prog = c('PAX8','NKX2-1','FOXE1','HHEX','UROD','DIO2','DUOX1','DUOX2','DUOXA1','DUOXA2',
               'SLC5A5','ANO1','SLC26A4','TPO','IYD','TG','PAX8','GLIS3','TSHR','SLC16A2','SLC16A10','EXOC4','ELMO1','VPS13C','STON2','SPG11')
  genesToKeep = genesToKeep[!genesToKeep %in% c(soupyGenes,thy_prog)]
}


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
outDir = paste0("~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/LRv1_fThy2n.REF/LRv1_fThy2n.REF.",annot_column,"_pPTC.tgt_",geneFilter)
#outDir = "~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/LRv1_fThy2n.REF/LRv1_fThy2n.REF.celltype_pPTC.tgt_maxGeneFilter"
model_fp = paste0('~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/LRv1_fThy2n.REF/fThy2n.REF.',annot_column,'_trainModel_',geneFilter,'_4kmaxcells_70perc_240403.RDS')
#model_fp = '~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/LRv1_fThy2n.REF/fThy2n.REF.celltype_trainModel_maxGeneFilter_4kmaxcells_70perc_240403.RDS'

LR_level='both'
srat_annot='group'
minGeneMatch = 0.99
maxCells=4000
tissue = 'thyroid'

out_prefix = paste0('fThy2n.REF.',annot_column,'_on_pThy.tgt_maxCells_70perc_240514_')
plot_prefix = paste0('fThy2n.REF.',annot_column,'_on_pThy.tgt_maxCells_70perc_240514_')
plot_prefix = NULL

outputs = runLR(REF.srat,ref_annot=ref_annot,srat=tgt.srat,LR_level=LR_level,srat_annot=srat_annot,
                model_fp = model_fp,outDir=outDir,plot_prefix=plot_prefix,out_prefix=out_prefix,
                minGeneMatch=minGeneMatch,maxCells = maxCells,tissue=tissue,
                scLR_TGTtype='',scLR_REFtype='',skipIfExists = skipIfExists)

message(sprintf('3. LR completed for tissue %s',tissue))




##---------------------------------##
##   Do some LR_similarity plots   ##
##---------------------------------##
model_fp = paste0('~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/LRv1_fThy2n.REF/fThy2n.REF.',annot_column,'_trainModel_',geneFilter,'_4kmaxcells_70perc_240403.RDS')
outputs = readRDS(file.path(paste0("~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/LRv1_fThy2n.REF/LRv1_fThy2n.REF.",annot_column,"_pPTC.tgt_",geneFilter),
                            paste0('fThy2n.REF.',annot_column,'_on_pThy.tgt_maxCells_70perc_240511_raw_LR_outputs.RDS')))

## annotated Cluster level #
if(length(outputs) >2){
  output = outputs[[1]]
}else{
  output = outputs[[2]][[1]]
}

type = ifelse(grepl('ref_',rownames(output)),'REF','TGT')
show_row_names = T
column_order = colnames(output)[order(colnames(output))]
row_order = rownames(output)[!grepl('Tumour',rownames(output))]
row_order = row_order[order(row_order)]
row_order = c(row_order,rownames(output)[grepl('Tumour',rownames(output))])


#plot_prefix = 'fThy2n.REF.celltype_on_pThy.tgt_maxCells_70perc_240511_'
plot_prefix = paste0('fThy2n.REF.',annot_column,'_on_pThy.tgt_maxCells_70perc_240511_')

pdf(file.path(outDir,paste0(plot_prefix,'clusterLR.pdf')),width = 12,height = 15)

hm = similarityHeatmap(output,
                       row_order=row_order,
                       column_order = column_order,
                       row_title_rot = 0,
                       row_title_gp = gpar(fontsize=10),row_names_gp = gpar(fontsize=10),row_names_max_width = unit(6,'cm'),
                       column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
                       split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = F)
draw(hm)
dev.off()


## annotated Single-cell level #
if(length(outputs) > 2){
  output = outputs[['scLR_all']]
}else{
  output = outputs[[2]][['scLR_all']]
}

#tgt.srat$annot = tgt.srat$finalAnn
#tgt.srat$annot2 = ifelse(tgt.srat$etiology == 'left_inferior_tumour',paste0('tum_',tgt.srat$annot),tgt.srat$annot)

in_mtx = output
type = paste0(#tgt.srat.harm$donor[match(rownames(output),tgt.srat.harm$cellID)],':',
              tgt.srat.harm$etiology[match(rownames(output),tgt.srat.harm$cellID)],':',
              as.character(tgt.srat.harm$seurat_clusters[match(rownames(output),tgt.srat.harm$cellID)]))
  # ifelse(as.character(tgt.srat$etiology[match(rownames(output),tgt.srat$cellID)]) == 'left_inferior_tumour',
  #             paste0('tum_',as.character(tgt.srat$seurat_clusters[match(rownames(output),tgt.srat$cellID)])),
  #             paste0('norm_',as.character(tgt.srat$seurat_clusters[match(rownames(output),tgt.srat$cellID)])))
type[is.na(type)] = paste0('ref_',as.character(REF.srat$annot[match(rownames(output)[is.na(type)],paste0('ref_',REF.srat$cellID))]))
type[type == 'NA:NA:NA'] = paste0('ref_',as.character(REF.srat$annot[match(rownames(output)[type == 'NA:NA:NA'],paste0('ref_',REF.srat$cellID))]))
type[type == 'NA:NA'] = paste0('ref_',as.character(REF.srat$annot[match(rownames(output)[type == 'NA:NA'],paste0('ref_',REF.srat$cellID))]))
#type = factor(type, levels = )

pdf(file.path(outDir,paste0(plot_prefix,'scLR_all.pdf')),width = 10,height = 50)

show_row_names=F
hm = similarityHeatmap(in_mtx,
                       column_order = colnames(output)[order(colnames(output))],
                       row_title_rot = 0,
                       row_title_gp = gpar(fontsize=5),#row_names_gp = gpar(fontsize=5),row_names_max_width = unit(6,'cm'),
                       column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
                       split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = F)
draw(hm)

dev.off()




##----- Add annotation    ------####
DimPlot(tgt.srat.harm,group.by = 'section',label = T,repel = T,label.box = T) + NoLegend()
## For every tumour cell, assign it to the reference cell type of best match....
mtx = outputs[['scLR_tgt']]
mtx = mtx[,!is.na(colSums(mtx))]
bestMatch = do.call(c,lapply(1:nrow(mtx),function(i){
  n = colnames(mtx)[which(mtx[i,] == max(mtx[i,]))]
  if(length(n) > 1){
    n = paste(n,collapse = ':')
  }
  return(n)
}))
names(bestMatch) = rownames(mtx)

table(tgt.srat$cellID %in% rownames(mtx))
tgt.srat$LRv1_fThy2n_pred = as.character(bestMatch[match(tgt.srat$cellID,names(bestMatch))])
tgt.srat.harm$LRv1_fThy2n_pred = as.character(bestMatch[match(tgt.srat.harm$cellID,names(bestMatch))])

a = as.data.frame(table(tgt.srat.harm$LRv1_fThy2n_pred,tgt.srat.harm$seurat_clusters))
a = a[a$Freq >0,]

##----    Manual annotation ------####
library(SoupX)
tgt.srat = tgt.srat.harm
qm = quickMarkers(tgt.srat@assays$RNA@counts,tgt.srat$seurat_clusters)
DimPlot(tgt.srat, group.by = 'seurat_clusters',cols = c(col25,pal34H),label = T,repel = T,label.box = T) + NoLegend()
FeaturePlot(tgt.srat,c('percent.mt','scrubScore'))
FeaturePlot(tgt.srat,'MS4A1')
DimPlot(tgt.srat,cells.highlight = tgt.srat$cellID[tgt.srat$annot_jul23 == 'endo_Capillary'])
DimPlot(tgt.srat,cells.highlight = tgt.srat$cellID[tgt.srat$seurat_clusters == 16])
tgt.srat$annot = '?'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(0,1,2,3,5,6,7,9,11,13,15,23,25,27,29,39,19,21)] = 'Tumour'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(13,25)] = 'Tumour_fTFC1'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(4,38,26,16)] = 'Thyrocytes'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(30,31)] = 'B_cells'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(10,28)] = 'T_cells'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(40)] = 'Mast.cells'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(8,33,20,22,35,37)] = 'Monocytes'

tgt.srat$annot[tgt.srat$seurat_clusters %in% c(14)] = 'SMCs'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(17)] = 'Mesenchymal'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(32)] = 'LECs'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(12,18)] = 'VECs'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(36)] = 'end_Capillary'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(34)] = 'end_Venus'

tgt.srat$annot[tgt.srat$seurat_clusters %in% c(24)] = 'doublets'


#tgt.srat$annot[tgt.srat$seurat_clusters %in% c(34,29)] = 'imm_Plasma.cells'




#tgt.srat$annot[tgt.srat$seurat_clusters %in% c(16,28,14)] = 'imm_ILCs'



## Import current annotation
current_pPTC = readRDS('~/lustre_mt22/Thyroid/Results/1_thyroid_annotation/thyroid_clean_ann_soupedRho0.2_jul23.RDS')
currentAnno = current_pPTC@meta.data

tgt.srat$annot_jul23 = current_pPTC$finalAnn[match(tgt.srat$cellID,current_pPTC$cellID)]
tgt.srat$annot_jul23[is.na(tgt.srat$annot_jul23)] = '-'
# tgt.srat$annot_apr24 = as.character(tgt.srat$annot)
# tgt.srat$annot_apr24[tgt.srat$cellID %in% current_pPTC$cellID[grepl('endo_',current_pPTC$finalAnn)]] = tgt.srat$annot_jul23[tgt.srat$cellID %in% current_pPTC$cellID[grepl('endo_',current_pPTC$finalAnn)]]
# tgt.srat$annot_apr24 = gsub('^end_','endo_',tgt.srat$annot_apr24)
# tgt.srat$annot_apr24[tgt.srat$cellID %in% current_pPTC$cellID & tgt.srat$seurat_clusters %in% c(14,31)] = tgt.srat$annot_jul23[tgt.srat$cellID %in% current_pPTC$cellID & tgt.srat$seurat_clusters %in% c(14,31)]
# tgt.srat$annot_apr24 = gsub('cells','cell',tgt.srat$annot_apr24)
# tgt.srat$annot_apr24 = gsub('MonoMac','Macrophage',tgt.srat$annot_apr24)
# tgt.srat$annot_apr24[tgt.srat$annot_apr24 == 'thy_Thyrocytes_LRRK2high'] = '?'
# tgt.srat$annot_apr24[tgt.srat$annot_apr24 == 'mes_Fibroblast'] = 'mes'
# tgt.srat$annot_apr24[tgt.srat$annot_apr24 == 'endo_Venus'] = 'endo_Venous'

tgt.srat$annot_may24 = tgt.srat$annot
tgt.srat$finalAnn = tgt.srat$annot_may24
tgt.srat$finalAnn_broad = tgt.srat$annot_may24

tgt.srat$finalAnn_broad[tgt.srat$finalAnn_broad %in% c('end_Capillary','end_Venus')] = 'VECs'
tgt.srat$finalAnn_broad[tgt.srat$finalAnn_broad %in% c('Monocytes')] = 'Myeloid_cells'
tgt.srat$finalAnn_broad[tgt.srat$finalAnn_broad %in% c('Tumour_fTFC1')] = 'Tumour'






keyHaemMarkers = c('CD34','CD38','SPINK2','MLLT3','PRSS57', # HSC_MPP
               'FLT3','RUNX2', 'LTB',# ELP
               'SERPINB1', 'GATA1',	'GATA2', 'TESPA1',	'CTNNBL1',#MEMP
               'FCER1A', 'ITGA2B', 'HBD','KLF1','PLEK', # MEP
               #'KIT',
               'CTSG',	'PRTN3', # CMP
               'AZU1','MPO','FLT3','PTPRC', # GMP
               'ZBTB16','LTB', 'CD52',# ILC precursor
               'IL7R','CD3D','GZMA', #Early.lymphoid_T.lymphocyte
               'NKG7','KLRD1', #NK
               'IGLL1','CD99', # Pre-pro
               'DNTT','CD79B','VPREB1','EBF1','CD19','RAG1',# pro-b-cells
               'MME','CD79A',# pre-b-cells
               'TCL1A','MME','RAG1','MS4A1',  # B-cells
               'ITGB2','SELL','ITGAM','CD14','FCGR3A',# Monocytes
               'CD68','MSR1',#Macrophages
               'IRF8',	'CLEC10A', #DC.precursor
               'CD1C','CLEC9A',#DC1
               'IL3RA', 'CLEC4C',# DC2 #pDC
               'CD27',#plasma cell?
               'CSF2RB','HDC','SERPINB1','TPSAB1','KIT', # Mast.cell
               'PF4','ITGA2B', #Megakaryocyte
               'GATA1','KLF1', # early Erythroid
               'ALAS2', # mid.erythroid
               'HBA1','BPGM', # late.erythroid
               'ESAM','PECAM1', 'KDR', 'PTPRB', 'PLVAP', # Endothelium
               'DCN','SERPINF1','COL3A1','BGN','COL1A1','VIM','ECM1', # endosteal fibroblast / fibroblast
               'APOA1','SCD','ALB','TTR' # Hepatocyte
)

DotPlot(tgt.srat,features = unique(keyHaemMarkers))+
  RotatedAxis()+
  theme(axis.text.x = element_text(size=7,vjust = 0.5,hjust = 1,angle = 90))

DotPlot(tgt.srat,group.by = 'finalAnn_broad',
        features = unique(c('NKX2-1','FOXE1','HHEX','UROD',
                     'CD302','BAG3','BLOC1S6','DIO2','DUOX1','DUOX2','DUOXA1','DUOXA2',
                     'SLC5A5','ANO1','SLC26A4','TPO','IYD','TG','PAX8','GLIS3','TSHR','SLC16A2','SLC16A10','EXOC4','ELMO1','VPS13C','STON2','SPG11',

                     'COL1A1', 'COL1A2', 'COL3A1', 'ACTA2','PTPRC','PECAM1', 'CD34', 'CDH5', 'VWF','HBA2','EPCAM', 'KRT18', 'KRT19','TG','TPO','EPCAM','KRT18','KRT19','S100A4', 'FN1', 'IGFBP6','TMSB4X'
))) + RotatedAxis() + theme(axis.text.x = element_text(size=7,vjust = 0.5,hjust = 1,angle = 90))





df = as.data.frame(table(tgt.srat$LRv1_fThy2n_pred,tgt.srat$annot_may24))
df=df[df$Freq >0,]
View(df)



DimPlot(tgt.srat,group.by = 'finalAnn_broad',cols = c(col25,pal34H,col25),label = T,label.box = T,repel = T) + NoLegend()
DimPlot(tgt.srat,cells.highlight = tgt.srat$cellID[tgt.srat$annot_lvl2 == '?'])

##---- Save Annotation results ------####

mdat = cbind(tgt.srat@meta.data,tgt.srat@reductions$umap@cell.embeddings)
write.csv(mdat,'~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/pPTC_clean_soupedRhoNone_may24_HARM_annotated_mdat.csv')

## Remove doublet cluster
tgt.srat = subset(tgt.srat,subset = annot_may24 != 'doublets')
tgt.srat = standard_clustering(tgt.srat,runHarmony = T,harmonyVar = c('donor'))


## Refine the annotation ##
DimPlot(tgt.srat,group.by = 'finalAnn_broad',cols=col25,label = T,repel = T)
DimPlot(tgt.srat,cells.highlight = tgt.srat$cellID[tgt.srat$section == 'thyroid_left_upper'])
table(tgt.srat$annot_jul23[tgt.srat$seurat_clusters == 26])
tgt.srat$annot_may24[tgt.srat$seurat_clusters %in% c(26,21,23)] = 'Plasma_cell'
tgt.srat$annot_may24[tgt.srat$annot_may24 == 'B_cells'] = 'B_cell'
tgt.srat$annot_may24[tgt.srat$annot_may24 == 'T_cells'] = 'T_cell'
tgt.srat$annot_may24[tgt.srat$annot_may24 == 'Mast.cells'] = 'Mast_cell'

tgt.srat$finalAnn = tgt.srat$annot_may24

tgt.srat$finalAnn_broad = tgt.srat$annot_may24
tgt.srat$finalAnn_broad[tgt.srat$finalAnn_broad %in% c('end_Capillary','end_Venus')] = 'VECs'
tgt.srat$finalAnn_broad[tgt.srat$finalAnn_broad %in% c('Monocytes')] = 'Myeloid_cell'
tgt.srat$finalAnn_broad[tgt.srat$finalAnn_broad %in% c('Tumour_fTFC1')] = 'Tumour'

mdat = cbind(tgt.srat@meta.data,tgt.srat@reductions$umap@cell.embeddings)
write.csv(mdat,'~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/pPTC_clean_soupedRhoNone_may24_HARM_annotated_noDoublets_mdat.csv')

saveRDS(tgt.srat,'~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/pPTC_clean_soupedRhoNone_may24_HARM_annotated_noDoublets.RDS')
tgt.srat = readRDS('~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/pPTC_clean_soupedRhoNone_may24_HARM_annotated_noDoublets.RDS')


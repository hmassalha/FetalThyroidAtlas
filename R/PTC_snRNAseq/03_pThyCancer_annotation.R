##--- Annotate 10X single-nuclei paediatric thyroid cancer dataset ---##
##    1. LR_v1; REF = fThy_2n; tgtSrat = pPTC
##    2. Harmony to integrator Tumour and normal biopsies


setwd('~/FetalThyroidAtlas/')


outDir = "Results/2505/PTC_snRNAseq/03_pThyCancer_annotation"
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}




##----------------##
##   Libraries  ####
##----------------##
library(Seurat)
library(tidyverse)
source("R/utils/logisticRegression.R")
source("R/utils/runLR.R")
source("R/utils/sc_utils.R")
source("R/utils/misc.R")



##-------------------------------------------------------------##
##   1. LRv1: LR with REF = fThyroid_2n, tgtSrat = pPTC      ####
##-------------------------------------------------------------##
seed = 2397
numPCs = 75
geneFilter = c('minGeneFilter','maxGeneFilter')
annot_column = c('celltype','cluster')

geneFilter = 'maxGeneFilter'
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
tgt.srat.fp = 'Results/2505/PTC_snRNAseq/03_pThyCancer_annotation/pPTC_clean_soupedXrhoLimNone_2505_HARM.RDS'
#'Results/2505/PTC_snRNAseq/02_pThyCancer_snPreprocessing/may25/scCancerThyroid_10X_rhoLimNone__clean_noMTCells.RDS'

if(file.exists(tgt.srat.fp)){
  tgt.srat = readRDS(tgt.srat.fp)
  # tgt.srat$group = tgt.srat$finalAnn_broad
  # tgt.srat$group[tgt.srat$group == 'Tumour'] = ifelse(tgt.srat$etiology[tgt.srat$group == 'Tumour'] == 'Met', 'Met:Y46',
  #                                                     paste0(tgt.srat$group[tgt.srat$group == 'Tumour'],':',tgt.srat$donor[tgt.srat$group == 'Tumour']))
  
  print(table(tgt.srat$donor))
}else{
  tgt.srat = readRDS('Results/2505/PTC_snRNAseq/02_pThyCancer_snPreprocessing/may25/scCancerThyroid_10X_rhoLimNone__clean_noMTCells.RDS')
  
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
  saveRDS(tgt.srat.harm,'Results/2505/PTC_snRNAseq/03_pThyCancer_annotation/pPTC_clean_soupedXrhoLimNone_2505_HARM.RDS')
  
  tgt.srat = standard_clustering(tgt.srat)
  saveRDS(tgt.srat,'Results/2505/PTC_snRNAseq/03_pThyCancer_annotation/pPTC_clean_soupedXrhoLimNone_2505.RDS')
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
  # # Remove soupy genes
  # soupyGenes = read.csv('~/lustre_mt22/Thyroid/Results/5_LR_fThyREF_on_pThy/fThy2nREF_downsampled_Thyrocytes_markersExprinpThy.csv')
  # soupyGenes = soupyGenes$gene[soupyGenes$frac_cellExpr >= 0.7]

  thy_prog = c('PAX8','NKX2-1','FOXE1','HHEX','UROD','DIO2','DUOX1','DUOX2','DUOXA1','DUOXA2',
               'SLC5A5','ANO1','SLC26A4','TPO','IYD','TG','PAX8','GLIS3','TSHR','SLC16A2','SLC16A10','EXOC4','ELMO1','VPS13C','STON2','SPG11')
  genesToKeep = genesToKeep[!genesToKeep %in% c(thy_prog)]
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
outDir = file.path(outDir,paste0("LRv1_fThy2n.REF/LRv1_fThy2n.REF.",annot_column,"_pPTC.tgt_",geneFilter))
#outDir = "~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/LRv1_fThy2n.REF/LRv1_fThy2n.REF.celltype_pPTC.tgt_maxGeneFilter"
model_fp = file.path(outDir,paste0('fThy2n.REF.',annot_column,'_trainModel_',geneFilter,'_4kmaxcells_70perc_2505.RDS'))
#model_fp = '~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/LRv1_fThy2n.REF/fThy2n.REF.celltype_trainModel_maxGeneFilter_4kmaxcells_70perc_240403.RDS'

LR_level='both'
srat_annot='seurat_clusters'
minGeneMatch = 0.99
maxCells=4000
tissue = 'thyroid'

out_prefix = paste0('fThy2n.REF.',annot_column,'_on_pThy.tgt_maxCells_70perc_2505_')
plot_prefix = paste0('fThy2n.REF.',annot_column,'_on_pThy.tgt_maxCells_70perc_2505_')
plot_prefix = NULL

outputs = runLR(REF.srat,ref_annot=ref_annot,srat=tgt.srat,LR_level=LR_level,srat_annot=srat_annot,
                model_fp = model_fp,outDir=outDir,plot_prefix=plot_prefix,out_prefix=out_prefix,
                minGeneMatch=minGeneMatch,maxCells = maxCells,tissue=tissue,
                scLR_TGTtype='',scLR_REFtype='',skipIfExists = skipIfExists)

message(sprintf('3. LR completed for tissue %s',tissue))




##---------------------------------##
##   Do some LR_similarity plots   ##
##---------------------------------##
model_fp = file.path(outDir,paste0('fThy2n.REF.',annot_column,'_trainModel_',geneFilter,'_4kmaxcells_70perc_2505.RDS'))
outputs = readRDS(file.path(outDir, paste0('fThy2n.REF.',annot_column,'_on_pThy.tgt_maxCells_70perc_2505_raw_LR_outputs.RDS')))

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
plot_prefix = paste0('fThy2n.REF.',annot_column,'_on_pThy.tgt_maxCells_70perc_2505_')

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
              tgt.srat$etiology[match(rownames(output),tgt.srat$cellID)],':',
              as.character(tgt.srat$seurat_clusters[match(rownames(output),tgt.srat$cellID)]))
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
DimPlot(tgt.srat,group.by = 'donor',label = T,repel = T,label.box = T) + NoLegend()
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
#tgt.srat.harm$LRv1_fThy2n_pred = as.character(bestMatch[match(tgt.srat.harm$cellID,names(bestMatch))])

a = as.data.frame(table(tgt.srat$LRv1_fThy2n_pred,tgt.srat$seurat_clusters))
a = a[a$Freq >0,]

##----    Manual annotation ------####
library(SoupX)
#tgt.srat = tgt.srat.harm
qm = quickMarkers(tgt.srat@assays$RNA@counts,tgt.srat$annot)
DimPlot(tgt.srat, group.by = 'annot',cols = c(col25,pal34H),label = T,repel = T,label.box = T) + NoLegend()
FeaturePlot(tgt.srat,c('percent.mt','scrubScore'))
FeaturePlot(tgt.srat,'MS4A1')
DimPlot(tgt.srat,cells.highlight = tgt.srat$cellID[tgt.srat$annot_jul23 == 'endo_Capillary'])
DimPlot(tgt.srat,cells.highlight = tgt.srat$cellID[tgt.srat$annot %in% c('34')])
tgt.srat$annot = '?'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(26,38,35,4,18,0,1,2,6,12,5,25,15,8,9,24,17,19)] = 'Tumour'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(39,3)] = 'Thyrocytes'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(29)] = 'B_cells'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(23,28)] = 'Plasma_cells'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(16,13)] = 'T_cells'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(11,14,36,33,41)] = 'VECs'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(31)] = 'LECs'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(27)] = 'Mesenchymal'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(10)] = 'SMCs'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(40)] = 'Mast_cells'
tgt.srat$annot[tgt.srat$seurat_clusters %in% c(20,22,7,30,32)] = 'Monocytes'
tgt.srat$annot[tgt.srat$annot == '?'] = as.character(tgt.srat$seurat_clusters[tgt.srat$annot == '?'])
tgt.srat$annot[tgt.srat$annot == '37'] = 'lowQual'


#tgt.srat$annot[tgt.srat$seurat_clusters %in% c(28)] = 'doublets'

# Subclustering non-cancer cells
s = subset(tgt.srat,annot != 'Tumour')
s = standard_clustering(s,runHarmony = T,harmonyVar='donor')
DimPlot(s, group.by = 'annot_tmp',cols = c(col25,pal34H),label = T,repel = T,label.box = T) + NoLegend()
s$annot_tmp = s$annot
s$annot_tmp[s$seurat_clusters %in% c(18,13,29,26,23,28,22,24,25)] = as.character(s$seurat_clusters[s$seurat_clusters %in% c(18,13,29,26,23,28,22,24,25)])
qm = quickMarkers(s@assays$RNA@counts,s$annot_tmp)

DimPlot(tgt.srat,cells.highlight = s$cellID[s$annot_tmp == '29'])
table(tgt.srat$seurat_clusters[tgt.srat$cellID %in% s$cellID[s$annot_tmp == '21']])

s$annot_tmp[s$annot_tmp == '13'] = 'B_cells'
s$annot_tmp[s$annot_tmp %in% c('21','34','37')] = 'doublets'
s$annot_tmp[s$annot_tmp == '22'] = 'endo_perivascular'
s$annot_tmp[s$annot_tmp == '24'] = 'capillary_VECs'
s$annot_tmp[s$annot_tmp == '25'] = 'vein_VECs'
s$annot_tmp[s$annot_tmp == '26'] = 'DC1'
s$annot_tmp[s$annot_tmp == '28'] = 'Fibroblasts'
s$annot_tmp[s$annot_tmp == '29'] = 'unknown'
s$annot_tmp[s$annot_tmp == '23'] = s$annot[s$annot_tmp == '23']

## Add this back to tgt.srat
table(tgt.srat$annot[tgt.srat$cellID %in% s$cellID[s$annot_tmp == 'unknown']])
DimPlot(tgt.srat,cells.highlight = tgt.srat$cellID[tgt.srat$cellID %in% s$cellID[s$annot_tmp == 'unknown']])
tgt.srat$annot[tgt.srat$cellID %in% s$cellID[s$annot_tmp == 'B_cells'] & tgt.srat$annot == '21'] = 'B_cells'
tgt.srat$annot[tgt.srat$cellID %in% s$cellID[s$annot_tmp == 'doublets']] = 'doublets'
tgt.srat$annot[tgt.srat$cellID %in% s$cellID[s$annot_tmp == 'DC1'] & tgt.srat$annot %in% c('21','34')] = 'DC1'
tgt.srat$annot[tgt.srat$cellID %in% s$cellID[s$annot_tmp == 'Fibroblasts'] & tgt.srat$annot %in% c('SMCs')] = 'Fibroblasts'
tgt.srat$annot[tgt.srat$cellID %in% s$cellID[s$annot_tmp == 'unknown']] = 'unknown'
tgt.srat$annot[tgt.srat$annot %in% c('21','34')] = 'doublets'

tgt.srat$celltype_detailed = as.character(tgt.srat$annot)
tgt.srat$celltype_detailed[tgt.srat$cellID %in% s$cellID[s$annot_tmp == 'endo_perivascular']] ='endo_perivascular'
tgt.srat$celltype_detailed[tgt.srat$cellID %in% s$cellID[s$annot_tmp == 'capillary_VECs']] ='capillary_VECs'
tgt.srat$celltype_detailed[tgt.srat$cellID %in% s$cellID[s$annot_tmp == 'vein_VECs']] ='vein_VECs'



## Check with dotplot



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


Idents(tgt.srat) = tgt.srat$annot
DotPlot(tgt.srat,features = unique(keyHaemMarkers))+
  RotatedAxis()+
  theme(axis.text.x = element_text(size=7,vjust = 0.5,hjust = 1,angle = 90))

DotPlot(tgt.srat,group.by = 'seurat_clusters',
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
# Harmony version
mdat = cbind(tgt.srat@meta.data,tgt.srat@reductions$umap@cell.embeddings)
mdat$annot[mdat$annot == 'Thyrocytes' & mdat$donor == 'Y24' & mdat$etiology =='Tumour'] = 'Tumour'
write.csv(mdat,file.path('Results/2505/PTC_snRNAseq/03_pThyCancer_annotation','pPTC_clean_soupedXrhoLimNone_annotated_2505_HARM_mdat.csv'))

# non-integrated version
srat = readRDS('Results/2505/PTC_snRNAseq/03_pThyCancer_annotation/pPTC_clean_soupedXrhoLimNone_2505.RDS')
srat@meta.data = cbind(srat@meta.data,mdat[match(srat$cellID,mdat$cellID),!colnames(mdat) %in% c(colnames(srat@meta.data),'UMAP_1','UMAP_2')])
DimPlot(srat,group.by = 'annot',cols = c(col25,pal34H),label = T,repel = T,label.box = T) + NoLegend()
mdat = cbind(srat@meta.data,srat@reductions$umap@cell.embeddings)
write.csv(mdat,file.path('Results/2505/PTC_snRNAseq/03_pThyCancer_annotation','pPTC_clean_soupedXrhoLimNone_annotated_2505_mdat.csv'))
saveRDS(srat,'Results/2505/PTC_snRNAseq/03_pThyCancer_annotation/pPTC_clean_soupedXrhoLimNone_annotated_2505.RDS')


# ## Remove doublet cluster
# srat = subset(srat,subset = annot != 'doublets')
# tgt.srat = standard_clustering(tgt.srat,runHarmony = T,harmonyVar = c('donor'))
# 
# 
# ## Refine the annotation ##
# DimPlot(tgt.srat,group.by = 'finalAnn_broad',cols=col25,label = T,repel = T)
# DimPlot(tgt.srat,cells.highlight = tgt.srat$cellID[tgt.srat$section == 'thyroid_left_upper'])
# table(tgt.srat$annot_jul23[tgt.srat$seurat_clusters == 26])
# tgt.srat$annot_may24[tgt.srat$seurat_clusters %in% c(26,21,23)] = 'Plasma_cell'
# tgt.srat$annot_may24[tgt.srat$annot_may24 == 'B_cells'] = 'B_cell'
# tgt.srat$annot_may24[tgt.srat$annot_may24 == 'T_cells'] = 'T_cell'
# tgt.srat$annot_may24[tgt.srat$annot_may24 == 'Mast.cells'] = 'Mast_cell'
# 
# tgt.srat$finalAnn = tgt.srat$annot_may24
# 
# tgt.srat$finalAnn_broad = tgt.srat$annot_may24
# tgt.srat$finalAnn_broad[tgt.srat$finalAnn_broad %in% c('end_Capillary','end_Venus')] = 'VECs'
# tgt.srat$finalAnn_broad[tgt.srat$finalAnn_broad %in% c('Monocytes')] = 'Myeloid_cell'
# tgt.srat$finalAnn_broad[tgt.srat$finalAnn_broad %in% c('Tumour_fTFC1')] = 'Tumour'
# 
# mdat = cbind(tgt.srat@meta.data,tgt.srat@reductions$umap@cell.embeddings)
# write.csv(mdat,'~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/pPTC_clean_soupedRhoNone_may24_HARM_annotated_noDoublets_mdat.csv')
# 
# saveRDS(tgt.srat,'~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/pPTC_clean_soupedRhoNone_may24_HARM_annotated_noDoublets.RDS')
# tgt.srat = readRDS('~/lustre_mt22/Thyroid/Results_v2/03_pThyCancer_annotation/pPTC_clean_soupedRhoNone_may24_HARM_annotated_noDoublets.RDS')
# 

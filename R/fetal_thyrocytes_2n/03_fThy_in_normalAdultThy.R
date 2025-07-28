## Finding evidence of TFC1 and TFC2 in publicly available adult normal thyrocyte datasets ##

main_outDir = "Results/2505/fThy_in_normal_aThy"
if(!dir.exists(main_outDir)){
  dir.create(main_outDir,recursive = T)
}

#setwd(main_outDir)
plotDir = 'Figures/2505'
if(!dir.exists(plotDir)){
  dir.create(plotDir,recursive = T)
}

##----------------##
##   Libraries  ####
##----------------##
library(Seurat)
library(tidyverse)
#source("~/lustre_mt22/Thyroid/scripts/final_script/forPublication/v1/helperFunctions.R")
source("R/utils/logisticRegression.R")
source("R/utils/runLR.R")


##-------------------------------------------------##
##    Import different adult scRNAseq datasets   ####
##-------------------------------------------------##
datasets_fp = c('Wang22'='Data/published_scRNAseq/Wang_etal_2022/Wang_etal_2022.RDS',
                   'Pu21'='Data/published_scRNAseq/Pu_etal_2021/Pu_etal_2021.RDS',
                   'Hong23'='Data/published_scRNAseq/Hong_etal_2023/Hong_etal_2023.RDS',
                   'Mosteiro23'='Data/published_scRNAseq/Mosteiro_etal_2023/Mosteiro_etal_2023.RDS',
                   'Lu23' = 'Data/published_scRNAseq/Lu_etal_2023/Lu_etal_2023.RDS')
if(!all(file.exists(datasets_fp))){
  warning(sprintf('These datasets do not exist, please check: %s',names(datasets_fp[!file.exists(datasets_fp)])))
  datasets_fp = datasets_fp[file.exists(datasets_fp)]
}

sratObj_list = list()


if('Mosteiro23' %in% names(datasets_fp)){
  aThy_Mosteiro23 = readRDS(datasets_fp['Mosteiro23'])
  aThy_Mosteiro23$dataset = 'aThy_Mosteiro23'
  sratObj_list[['aThy_Mosteiro23']] = aThy_Mosteiro23
}

if('Hong23' %in% names(datasets_fp)){
  #aThy_Hong23 = readRDS('~/lustre_mt22/Thyroid/Data/published_scRNAseq/Hong_etal_2023/GSE182416_sratObj.RDS')
  aThy_Hong23 = readRDS(datasets_fp['Hong23'])
  aThy_Hong23$annot = as.character(aThy_Hong23$celltype_sub)
  aThy_Hong23$dataset = 'aThy_Hong23'
  sratObj_list[['aThy_Hong23']] = aThy_Hong23
}

if('Pu21' %in% names(datasets_fp)){
  # This dataset contains tumour and para-tumour (normal), so need to subset to just keep "normal" cells
  aThy_Pu21 = readRDS(datasets_fp['Pu21'])
  aThy_Pu21 <- subset(aThy_Pu21, subset = orig.ident %in% c('PTC1.P','PTC2.P','PTC3.P','PTC8.P','PTC9.P'))
  aThy_Pu21 = standard_clustering(aThy_Pu21)
  aThy_Pu21$celltype[aThy_Pu21$celltype == 'thy_Lumen-forming'] = 'fTFC2'
  aThy_Pu21$annot = as.character(aThy_Pu21$celltype)
  aThy_Pu21$dataset = 'aThy_Pu21'
  sratObj_list[['aThy_Pu21']] = aThy_Pu21
}


if('Wang22' %in% names(datasets_fp)){
  # This dataset contains tumour and para-tumour (normal), so need to subset to just keep "normal" cells
  aThy_Wang22 = readRDS(datasets_fp['Wang22'])
  aThy_Wang22 = subset(aThy_Wang22,stim %in% c('NT'))
  aThy_Wang22$annot = as.character(aThy_Wang22$celltype)
  aThy_Wang22$dataset = 'aThy_Wang22'
  sratObj_list[['aThy_Wang22']] = aThy_Wang22
}



if('Lu23' %in% names(datasets_fp)){
  
  aThy_Lu23 = readRDS(datasets_fp['Lu23'])
  
  # This dataset contains tumour and para-tumour (normal), so need to subset to just keep "normal" cells
  aThy_Lu23$condition = ifelse(grepl('NORM',aThy_Lu23$sample.ID),'normal',
                               ifelse(grepl('PTC',aThy_Lu23$sample.ID),'PTC','others'))
  aThy_Lu23 = subset(aThy_Lu23,condition == 'normal')
  aThy_Lu23 = standard_clustering(aThy_Lu23)
  aThy_Lu23$dataset = 'aThy_Lu23'
  aThy_Lu23$annot = as.character(aThy_Lu23$celltype)
  
  sratObj_list[['aThy_Lu23']] = aThy_Lu23
}






##-------------------------------------------##
##    Train LR model on fetal thyrocytes   ####
##-------------------------------------------##
seed = 2397
numPCs = 75
geneFilter = 'minGeneFilter'
annot_column = 'celltype'
ref_celltypes = 'thyrocytes_only'
with_SCPs = T

##---- Import REF dataset ------##
## Import fetal 2n thyrocytes only
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
    # Add UMAP_coords
    fThyroid = standard_clustering(fThyroid)
    umap = read.csv('/nfs/team292/Thyroid_hm11_mt22/fThyroid_2n/fThyroid_2n_UMAP.csv.gz')
    umap = column_to_rownames(umap,'X')
    fThyroid@reductions$umap@cell.embeddings[,'UMAP_1'] = umap$UMAP_1[match(rownames(fThyroid@meta.data),rownames(umap))]
    fThyroid@reductions$umap@cell.embeddings[,'UMAP_2'] = umap$UMAP_2[match(rownames(fThyroid@meta.data),rownames(umap))]
    DimPlot(fThyroid,group.by = 'cluster',label = T,label.box = T,repel = T) + NoLegend()
    saveRDS(fThyroid,'Data/fThyroid_2n_atlas.RDS')
    
    
    ## Similarly for 2n-T21 object
    # Create a seurat object of just thyrocytes from diploid fetuses
    fThyroid_2nT21 = Read10X('/nfs/team292/Thyroid_hm11_mt22/fThyroid_2nT21/mtx/')
    fThyroid = CreateSeuratObject(fThyroid)
    # Add metadata
    mdat = read.csv('/nfs/team292/Thyroid_hm11_mt22/fThyroid_2n/fThyroid_2n_obs.csv.gz',row.names = 1)
    table(rownames(mdat) %in% rownames(fThyroid@meta.data))
    table(colnames(fThyroid@meta.data) %in% colnames(mdat))
    fThyroid@meta.data =  cbind(fThyroid@meta.data,mdat[match(rownames(fThyroid@meta.data),rownames(mdat)),!colnames(mdat) %in% colnames(fThyroid@meta.data)])
    # Add UMAP_coords
    fThyroid = standard_clustering(fThyroid)
    umap = read.csv('/nfs/team292/Thyroid_hm11_mt22/fThyroid_2n/fThyroid_2n_UMAP.csv.gz')
    umap = column_to_rownames(umap,'X')
    fThyroid@reductions$umap@cell.embeddings[,'UMAP_1'] = umap$UMAP_1[match(rownames(fThyroid@meta.data),rownames(umap))]
    fThyroid@reductions$umap@cell.embeddings[,'UMAP_2'] = umap$UMAP_2[match(rownames(fThyroid@meta.data),rownames(umap))]
    DimPlot(fThyroid,group.by = 'cluster',label = T,label.box = T,repel = T) + NoLegend()
    saveRDS(fThyroid,'Data/fThyroid_2n_atlas.RDS')
    
    
    
    
    fThy = subset(fThyroid,subset = cluster == 'Thyrocytes')
    fThy = standard_clustering(fThy)
    saveRDS(fThy,fThy_fp)
}else{
  fThy = readRDS(fThy_fp)  
  fThy$cellID = rownames(fThy@meta.data)
}

fThy$finalAnn = fThy[[annot_column]]
fThy$cellID = colnames(fThy)
fThy = subset(fThy,subset = cellID %in% fThy$cellID[!fThy$celltype %in% c('thy_Cycling')])
fThy$finalAnn[fThy$celltype == 'thy_TH_processing'] = 'fTFC1'
fThy$finalAnn[fThy$celltype == 'thy_Lumen-forming'] = 'fTFC2'
table(fThy$finalAnn)

if(with_SCPs){
  ## Import fAdr
  fAdr = readRDS('~/reference_datasets/fetalAdrenal_Kildisuite21/adr_all.rds')# 10X indexed GRCh38 1.2.0 reference
  fAdr@meta.data$cellID = rownames(fAdr@meta.data)
  # Subset to just SCPs
  fAdr = subset(fAdr,Annotation == 'SCPs')
  fAdr$annot = fAdr$Annotation
  fAdr$finalAnn = fAdr$Annotation
  
  ## Combine to create REF.srat
  REF.srat = merge_seurat_objects(fThy,fAdr,keepAllGenes = F,genomeVersions = c('v38','v38'))
  ## Fix up metadata
  REF.srat$dataset = ifelse(REF.srat$cellID %in% fThy$cellID,'fThy',
                            ifelse(REF.srat$cellID %in% fAdr$cellID,'fAdr','others'))
  
}else{
  REF.srat = fThy
}



REF.srat$annot = REF.srat$finalAnn
table(REF.srat$annot)

message('REF.srat loaded')





##---- Configure REF.srat and tgt.srat for LRv1 ------##
# Set tgt.srat temporarily as one of the adult thyrocytes object 
tgt.srat = sratObj_list[[1]]

tgt.genes = do.call(c,lapply(1:length(sratObj_list),function(i){return(rownames(sratObj_list[[i]]))}))
tgt.genes = table(tgt.genes)
tgt.genes = tgt.genes[tgt.genes == length(sratObj_list)]
tgt.genes = names(tgt.genes)
# Also expressed in REF.srat
genesToKeep = intersect(tgt.genes,rownames(REF.srat))

# Remove rubbish genes
genesToKeep = genesToKeep[!grepl('^MT-|^RPL|^RPS|^MALAT1$|^NEAT1$|^AC\\d+|^AL\\d+',genesToKeep)]
length(genesToKeep)

message(sprintf('Keep %d genes for LRv1 training',length(genesToKeep)))

## Subset REF.srat and tgt.srat to keep genes of interest only
REF.srat@assays$RNA@counts = REF.srat@assays$RNA@counts[genesToKeep,]
tgt.srat@assays$RNA@counts = tgt.srat@assays$RNA@counts[genesToKeep,]





##---- Run Logistic Regression ------##

skipIfExists=F
ref_annot='annot'
maxCells=4000

LR_level='both'
srat_annot='annot'
minGeneMatch = 0.99
maxCells=4000
alpha=0.1
tissue = 'thyroid'

outDir = file.path(main_outDir,sprintf("LRv1_fThy.REF_aThy.tgt_alpha.%s_maxCell.%dk_%s",alpha,maxCells/1000,geneFilter))

#[USE THIS] REF classess: fTFC1/2 + SCPs, but genes_toKeep only based on aThy
model_fp = file.path(outDir,sprintf('LRv1_fThy.REF_trainModel_alpha.%s_maxCell.%dk_%s.RDS',alpha,maxCells/1000,geneFilter))


out_prefix = sprintf('fThy.REF_on_aThy.tgt_alpha.%s_maxCell.%dk_%s_',alpha,maxCells/1000,geneFilter)
plot_prefix = out_prefix
plot_prefix = NULL

outputs = runLR(REF.srat,ref_annot=ref_annot,srat=tgt.srat,LR_level=LR_level,srat_annot=srat_annot,
                model_fp = model_fp,outDir=outDir,plot_prefix=plot_prefix,out_prefix=out_prefix,
                minGeneMatch=minGeneMatch,maxCells = maxCells,tissue=tissue,skipIfExists = F,
                scLR_TGTtype='',scLR_REFtype='',alpha=alpha)

message(sprintf('LR completed for tissue %s',tissue))













##----------------------------------------------------------##
##    Compute similarity scores in adult thyrocytes data  ####
##----------------------------------------------------------##
outDir = file.path(main_outDir,sprintf("LRv1_fThy.REF_aThy.tgt_alpha.%s_maxCell.%dk_%s",alpha,maxCells/1000,geneFilter))
model_fp = file.path(outDir,sprintf('LRv1_fThy.REF_trainModel_alpha.%s_maxCell.%dk_%s.RDS',alpha,maxCells/1000,geneFilter))

#model_fp = '~/lustre_mt22/Thyroid/Results_v3/03_fThy_in_normalAdultThy/LRv1_fThy.REF_aThy.tgt_alpha.0.1_maxCell.4k_minGeneFilter/LRv1_fThy.REF_trainModel_alpha.0.1_maxCell.4k_minGeneFilter.RDS'
model = readRDS(model_fp)

## Predict similarity across all adult datasets
logit_mtx = lapply(1:length(sratObj_list),function(i){ 
  print(i)
  mtx = predictSimilarity(model, sratObj_list[[i]]@assays$RNA@counts,logits = T,minGeneMatch = 0.8)
  return(mtx)})
names(logit_mtx) = names(sratObj_list)
logit_mtx = do.call(rbind,logit_mtx)
colnames(logit_mtx)

## Subset for just TFC1/2 reference
tfc_mtx = logit_mtx[,c('fTFC1','fTFC2','SCPs')]
tfc_mtx = as.data.frame(tfc_mtx)
tfc_mtx$cellID = rownames(tfc_mtx)

## Add values back to seurat_objects
for(i in 1:length(sratObj_list)){
  sratObj_list[[i]]$LRv1_fThyfull_fTFC1 = tfc_mtx$fTFC1[match(sratObj_list[[i]]$cellID,tfc_mtx$cellID)]
  sratObj_list[[i]]$LRv1_fThyfull_fTFC2 = tfc_mtx$fTFC2[match(sratObj_list[[i]]$cellID,tfc_mtx$cellID)]
  sratObj_list[[i]]$LRv1_fThyfull_SCP = tfc_mtx$SCPs[match(sratObj_list[[i]]$cellID,tfc_mtx$cellID)]
  sratObj_list[[i]]$LRv1_fThyfull_diff = (sratObj_list[[i]]$LRv1_fThyfull_fTFC1 - sratObj_list[[i]]$LRv1_fThyfull_fTFC2)  
  
  # if(i == 2){
  #   diff = (sratObj_list[[i]]$LRv1_fThyfull_fTFC1 - sratObj_list[[i]]$LRv1_fThyfull_fTFC2)
  #   sratObj_list[[i]]$LRv1_fThyfull_diff = diff[match(sratObj_list[[i]]$cellID,names(diff))]
  #   
  # }else{
  #   sratObj_list[[i]]$LRv1_fThyfull_diff = (sratObj_list[[i]]$LRv1_fThyfull_fTFC1 - sratObj_list[[i]]$LRv1_fThyfull_fTFC2)  
  # }
  # 
}


## Subclustering only thyrocytes
thyrocytesObj_list = list()
for(i in 1:length(sratObj_list)){
  s = subset(sratObj_list[[i]],subset = annot %in% c('aTFC1','aTFC2','aTFC3','aTFC4','aTFC5',
                                                     'C0','C1','C2','C3','C4',
                                                     'fTFC1','fTFC2',
                                                     'Follicular_cells','Epithelial cell'))
  s = standard_clustering(s,runHarmony = T,harmonyVar = 'orig.ident')
  thyrocytesObj_list[[i]] = s
}

names(thyrocytesObj_list) = names(sratObj_list)
DimPlot(thyrocytesObj_list[[5]],group.by = 'annot',cols = col25,label = T,repel = T,label.box = T)
DimPlot(thyrocytesObj_list[[4]],group.by = 'seurat_clusters',cols = col25,label = T,repel = T,label.box = T)
DimPlot(thyrocytesObj_list[[4]],cells.highlight = thyrocytesObj_list[[4]]$cellID[thyrocytesObj_list[[4]]$seurat_clusters %in% c(6,13,9,8,12)])

# Pu21: remove clustes 6,13,9,8,12 - as these look like immune / thyrocytes doublets
clusters_toKeep = unique(thyrocytesObj_list[[which(names(thyrocytesObj_list) == 'aThy_Pu21')]]$seurat_clusters)
clusters_toKeep = clusters_toKeep[!clusters_toKeep %in% c(6,13,9,8,12)]
thyrocytesObj_list[[which(names(thyrocytesObj_list) == 'aThy_Pu21')]] = subset(thyrocytesObj_list[[which(names(thyrocytesObj_list) == 'aThy_Pu21')]],subset=seurat_clusters %in% clusters_toKeep)
thyrocytesObj_list[[which(names(thyrocytesObj_list) == 'aThy_Pu21')]] = standard_clustering(thyrocytesObj_list[[which(names(thyrocytesObj_list) == 'aThy_Pu21')]],runHarmony = T,harmonyVar = 'orig.ident')

# Lu23: regenerate UMAP with less tight clustering
s = thyrocytesObj_list[[which(names(thyrocytesObj_list) == 'aThy_Lu23')]]
DimPlot(s,group.by = 'seurat_clusters',cols=col25,label = T,label.box = T)
s = FindClusters(s, resolution = 0.2, reduction = "harmony")
s = RunUMAP(s, dims=1:20, reduction ="harmony",min.dist=0.45)
m = quickMarkers(s@assays$RNA@counts,s$seurat_clusters)

thyrocytesObj_list[[which(names(thyrocytesObj_list) == 'aThy_Lu23')]] = s

saveRDS(thyrocytesObj_list,file.path(main_outDir,'published_aThyrocytesOnly_5datasets.RDS'))
thyrocytesObj_list = readRDS(file.path(main_outDir,'published_aThyrocytesOnly_5datasets.RDS'))


##-------------------------------##
##       Generating plots      ####
##-------------------------------##
source('R/utils/misc.R')

data_toPlot = do.call(rbind,lapply(1:length(thyrocytesObj_list),function(i){
  tmp = cbind(thyrocytesObj_list[[i]]@meta.data[,c('cellID','annot','LRv1_fThyfull_fTFC1','LRv1_fThyfull_fTFC2','LRv1_fThyfull_diff','dataset')],
              thyrocytesObj_list[[i]]@reductions$umap@cell.embeddings)
  return(tmp)
}))

# data_toPlot$score = ifelse(data_toPlot$LRv1_fThyfull_fTFC2 > 5, 5,
#                            ifelse(data_toPlot$LRv1_fThyfull_fTFC2 < -5,-5,data_toPlot$LRv1_fThyfull_fTFC2))
data_toPlot$score = ifelse(data_toPlot$LRv1_fThyfull_diff > 5, 5,
                           ifelse(data_toPlot$LRv1_fThyfull_diff < -5,-5,data_toPlot$LRv1_fThyfull_diff))
data_toPlot$cell_assignment = ifelse(data_toPlot$score >= 1.5,'TFC1',
                                     ifelse(data_toPlot$score <= -1.5,'TFC2','ambiguous'))

write.csv(data_toPlot,file.path(main_outDir,'adultThyrocytes_5datasets_LRv1.fThy.full.csv'))



## Plot barplot showing number of cells
## Logits are log-of-odds. Odds = p/(1-p) --> probability of success over failures
## Difference in logits --> log of Odd1/Odd 2 --> log of odd ratio (OR). 
## This implies logit_difference of 1.6 means that the Odd1 is 5 times more likely than Odd2
## Compare odds might be more useful than comparing probabilities because of the 0-1 limits of probability
## The probability might not be too disimilar, but the OR might be close to Infinity or -Inf
## --> high "score" --> odds of being a fTFC1 is much higher than the odds of the cell being fTFC2, and vice versa 

table(data_toPlot$cell_assignment,data_toPlot$dataset)
dd = data_toPlot %>% group_by(dataset,cell_assignment) %>% summarise(nCell = n()) %>% 
  group_by(dataset) %>% mutate(totalCell = sum(nCell),frac = nCell/totalCell)
dd$cell_assignment = factor(dd$cell_assignment,rev(c('TFC1','TFC2','ambiguous')))



plotFun_aThy_LRv1.fThy_frac.Cell = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(0.1,0.1,1,0.1))
  p1 = ggplot(dd,aes(x=frac,y=dataset,fill=cell_assignment))+
    geom_col(width = 0.75) +
    scale_x_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1))+
    scale_fill_manual(values = c('TFC2'='#64b9c6','ambiguous'=grey(0.8),'TFC1'='#356f87'))+
    theme_classic() + 
    theme(#panel.border = element_rect(fill = F,colour = 'black',linewidth = 1),
      axis.line = element_blank(),
      strip.background = element_blank(),
      axis.text = element_text(colour='black'),
      axis.ticks = element_blank())+xlab('') + ylab('')
  
  print(p1)
}

saveFig(file.path(plotDir,'Fig2_aThy_LRv1.fThy_fracCell'),plotFun_aThy_LRv1.fThy_frac.Cell,rawData=dd,width = 4.5,height = 2.8,res = 500)  
  


plotFun_aThy_UMAP_LRv1.fThy = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(0.1,0.1,1,0.1))
  p1 = ggplot(data_toPlot[data_toPlot$cell_assignment %in% c('TFC1'),],aes(UMAP_1,UMAP_2,col=score))+
    geom_point(size=0.1,alpha=0.9) + 
    geom_point(data = data_toPlot[data_toPlot$cell_assignment %in% c('TFC2','ambiguous'),],size=0.1,alpha=0.9) + 
    facet_wrap(vars(dataset),scales = 'free',nrow=1)+
    #scale_color_gradient2(low='#64b9c6',mid='white',high='#356f87')+
    scale_color_gradient2(low='black',mid='white',high='#ED820E')+
    theme_classic() + 
    theme(panel.border = element_rect(fill = F,colour = 'black',linewidth = 0.5),
      axis.line = element_blank(),
      strip.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank())+xlab('') + ylab('')
  
  if(noFrame & !noPlot){
    p1 = p1 + theme(strip.text = element_blank())
  }else if(!noFrame & noPlot){
    dd = rbind(data_toPlot[sample(1:nrow(data_toPlot),50),],
               data_toPlot[data_toPlot$score == max(data_toPlot$score),][1:50,],
               data_toPlot[data_toPlot$score == min(data_toPlot$score),][1:50,])
    p1 = ggplot(dd,aes(UMAP_1,UMAP_2,col=score))+
      geom_point(size=0.1,alpha=0.9) + 
      facet_wrap(vars(dataset),scales = 'free',nrow=1)+
      #scale_color_gradient2(low='#64b9c6',mid='white',high='#356f87')+
      scale_color_gradient2(low='black',mid='white',high='#ED820E')+
      theme_classic() + 
      theme(panel.border = element_rect(fill = F,colour = 'black',linewidth = 0.5),
        axis.line = element_blank(),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+xlab('') + ylab('')
  }
  
  print(p1)
}

saveFig(file.path(plotDir,'Fig2_aThy_UMAP_LRv1.fThy'),plotFun_aThy_UMAP_LRv1.fThy,rawData=data_toPlot,width = 15,height = 3.3,res = 500)  



## Dot plots of marker genes per cell_type assignment

for(i in 1:length(thyrocytesObj_list)){
  thyrocytesObj_list[[i]]$cell_assignment = data_toPlot$cell_assignment[match(thyrocytesObj_list[[i]]$cellID,data_toPlot$cellID)]
  Idents(thyrocytesObj_list[[i]]) = thyrocytesObj_list[[i]]$cell_assignment
  thyrocytesObj_list[[i]]$cell_assignment = factor(thyrocytesObj_list[[i]]$cell_assignment,c('TFC1','TFC2','ambiguous'))
}

genes_toPlot = c('NKX2-1','HHEX','FOXE1','PAX8','GLIS3','TSHR')
library(gridExtra)
plotFun_dotPlot_hor = function(noFrame=FALSE,noPlot=FALSE){
  plot_list = list()
  for(i in 1:length(thyrocytesObj_list)){
    p = DotPlot(thyrocytesObj_list[[i]],group.by = 'cell_assignment',
                 features = genes_toPlot,scale.min = 10,scale.max = 90)+
      RotatedAxis()+
      scale_y_discrete(position = "left")+
      scale_x_discrete(position = "bottom")+
      scale_color_gradientn(colours = c('#F1F6EB','#85BFB2','#1B4E89'))+
      theme(axis.text.y = element_text(size=11),
            axis.text.x = element_text(size=7,angle = 90,vjust = 0.5,hjust = 1,face='italic'),
            legend.title = element_text(size=8),
            legend.text = element_text(size=8),
            legend.position = 'right',plot.margin = unit(c(0,0,0,0),'cm')) + xlab('') + ylab('')
    plot_list[[names(thyrocytesObj_list)[i]]] = p
  }
  
  do.call("grid.arrange", c(plot_list, ncol=5))

}
saveFig(file.path(plotDir,'Fig2_aThy_LRv1.fThy_DotPlot_hor'),plotFun_dotPlot_hor,width = 18,height = 1.7,res = 500)
saveFig(file.path(plotDir,'Fig2_aThy_LRv1.fThy_DotPlot_hor_big'),plotFun_dotPlot_hor,width = 16,height = 3,res = 500)





































## Plot PAX8 expression level 


## Plot expression of PAX8 / GLIS3 in adult thyrocytes
library(SoupX)
thyrocytesObj_list = readRDS(file.path(main_outDir,'published_aThyrocytesOnly_5datasets.RDS'))


data_toPlot = do.call(rbind,lapply(1:length(thyrocytesObj_list),function(i){
  s = thyrocytesObj_list[[i]]
  df = do.call(cbind,list(s@meta.data[,c('cellID','dataset','annot','cell_assignment')],
            s@reductions$umap@cell.embeddings,
            t(as.matrix(s@assays$RNA@data[c('PAX8','GLIS3'),]))))
  return(df)
}))


dd = data_toPlot[data_toPlot$dataset == 'aThy_Mosteiro23',]
plotFun_aThy_UMAP_LRv1.fThy_markerExpr_Mosteiro23 = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(0.1,0.1,1,0.1))
  p1 = ggplot(dd,aes(UMAP_1,UMAP_2,col=PAX8))+
    geom_point(size=0.1,alpha=0.9) + 
    facet_wrap(vars(dataset),scales = 'free',nrow=1)+
    #scale_color_gradient2(low='#64b9c6',mid='white',high='#356f87')+
    #scale_color_gradientn(colours = c(grey(0.8),'#7ec0ee','#2e4472'))+
    scale_color_gradientn(colours = c(grey(0.85),'#85bfb2','#102E52'))+
    theme_classic() + 
    theme(#panel.border = element_rect(fill = F,colour = 'black',linewidth = 1),
      axis.line = element_blank(),
      strip.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank())+xlab('') + ylab('')
  
  
  if(noFrame & !noPlot){
    p1 = p1 + theme(strip.text = element_blank())
  }else if(!noFrame & noPlot){
    dd = rbind(dd[sample(1:nrow(dd),50),],
               dd[dd$PAX8 == max(dd$PAX8),][1:50,],
               dd[dd$PAX8 == min(dd$PAX8),][1:50,])
    p1 = ggplot(dd,aes(UMAP_1,UMAP_2,col=PAX8))+
      geom_point(size=0.1,alpha=0.9) + 
      facet_wrap(vars(dataset),scales = 'free',nrow=1)+
      scale_color_gradientn(colours = c(grey(0.85),'#85bfb2','#102E52'))+
      theme_classic() + 
      theme(#panel.border = element_rect(fill = F,colour = 'black',linewidth = 1),
        axis.line = element_blank(),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+xlab('') + ylab('')
  }
  
  print(p1)
}

saveFig(file.path(plotDir,'FigXX_aThy_UMAP_LRv1.fThy_PAX8expr_Mosteiro23'),plotFun_aThy_UMAP_LRv1.fThy_markerExpr_Mosteiro23,rawData=dd,width = 3.7,height = 3,res = 500)  



dd = data_toPlot
plotFun_aThy_UMAP_LRv1.fThy_markerExpr_all = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(0.1,0.1,1,0.1))
  p1 = ggplot(dd,aes(UMAP_1,UMAP_2,col=PAX8))+
    geom_point(size=0.1,alpha=0.9) + 
    facet_wrap(vars(dataset),scales = 'free',nrow=1)+
    #scale_color_gradient2(low='#64b9c6',mid='white',high='#356f87')+
    #scale_color_gradientn(colours = c(grey(0.8),'#7ec0ee','#2e4472'))+
    scale_color_gradientn(colours = c(grey(0.85),'#85bfb2','#102E52'))+
    theme_classic() + 
    theme(#panel.border = element_rect(fill = F,colour = 'black',linewidth = 1),
      axis.line = element_blank(),
      strip.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank())+xlab('') + ylab('')
  
  
  if(noFrame & !noPlot){
    p1 = p1 + theme(strip.text = element_blank())
  }else if(!noFrame & noPlot){
    dd = rbind(dd[sample(1:nrow(dd),50),],
               dd[dd$PAX8 == max(dd$PAX8),][1:50,],
               dd[dd$PAX8 == min(dd$PAX8),][1:50,])
    p1 = ggplot(dd,aes(UMAP_1,UMAP_2,col=PAX8))+
      geom_point(size=0.1,alpha=0.9) + 
      facet_wrap(vars(dataset),scales = 'free',nrow=1)+
      scale_color_gradientn(colours = c(grey(0.85),'#85bfb2','#102E52'))+
      theme_classic() + 
      theme(#panel.border = element_rect(fill = F,colour = 'black',linewidth = 1),
        axis.line = element_blank(),
        strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+xlab('') + ylab('')
  }
  
  print(p1)
}

saveFig(file.path(plotDir,'FigXX_aThy_UMAP_LRv1.fThy_PAX8expr_5datasets'),plotFun_aThy_UMAP_LRv1.fThy_markerExpr_all,rawData=dd,width = 15,height = 3,res = 500)  


p2 = ggplot(data_toPlot[data_toPlot$dataset == 'aThy_Wang22',],aes(UMAP_1,UMAP_2,col=PAX8))+
  geom_point(size=0.3,alpha=0.9) + 
  #geom_point(data = data_toPlot[data_toPlot$cell_assignment %in% c('TFC2','ambiguous'),],size=0.1,alpha=0.9) + 
  facet_wrap(vars(dataset),scales = 'free')+
  #scale_color_gradient2(low='#64b9c6',mid='white',high='#356f87')+
  #scale_color_gradient2(low='white',mid=grey(0.5),high='black')+
  scale_color_gradientn(colours = c(grey(0.85),'#85bfb2','#102E52'))+
  theme_classic() + 
  theme(#panel.border = element_rect(fill = F,colour = 'black',linewidth = 1),
    axis.line = element_blank(),
    strip.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank())+xlab('') + ylab('')

p2

















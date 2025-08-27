##--- Scoring for the enrichment of fetal thyrocyte cell states (fTFC1 and fTFC2) signals in adult and paediatric PTC ---##
##    0. Process inhouse bulk RNA-seq data
##    1. Define fetal cell state specific gene module
##    2. Score for the enrichment of these modules in adult + children PTC

setwd('~/FetalThyroidAtlas/')

outDir = "Results/2505/PTC_bulkRNAseq"
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

plotDir = 'Figures/2505'




##----------------##
##   Libraries  ####
##----------------##
library(tidyverse)
library(GenomicFeatures)
source("R/utils/misc.R")



##----------------------------##
##   Set Global parameters  ####
##----------------------------##

#Define genomic coordinates
gtf = '/nfs/cellgeni/STAR/human/2020A-full/GRCh38_v32_modified.gtf'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)

gene_map <- rtracklayer::import(gtf)  %>% 
  GenomeInfoDb::keepSeqlevels(
    paste0('chr',c(1:22, "X", "Y", "M")),
    pruning.mode = "coarse"
  ) %>% 
  subset(type == "gene")

names(gene_map) <- gene_map$gene_id
write.csv(as.data.frame(gene_map),'~/FetalThyroidAtlas/Results/geneMap.csv')


##---------------------------------------------##
##   0. Process inhouse bulk RNA-seq data    ####
##---------------------------------------------##
bulk_raw_counts = read.delim('Data/inhouse_bulk/results/expression_tables/salmon_reads.counts.tsv',sep = '\t')
rownames(bulk_raw_counts) = bulk_raw_counts$ensemblID
bulk_raw_TPM = read.delim('Data/inhouse_bulk/results/expression_tables/salmon_reads.TPM.tsv',sep = '\t')
rownames(bulk_raw_TPM) = bulk_raw_TPM$ensemblID

row_data = bulk_raw_counts[,c("ensemblID", "geneSymbol", "geneLength","effLength")]

col_data = readxl::read_excel('~/projectManifest.xlsx',sheet = 'bulk_thyroid')
dim(col_data)
col_data$sampleID = paste0('X',col_data$sangerSampleID)
col_data = col_data[col_data$sampleID %in% colnames(bulk_raw_counts),]
dim(col_data)
table(colnames(bulk_raw_counts)[!colnames(bulk_raw_counts) %in% c("ensemblID", "geneSymbol", "geneLength","effLength")] %in% paste0('X',col_data$sangerSampleID))

col_data = as.data.frame(col_data)
rownames(col_data) = col_data$sampleID
col_data$source = 'Sanger'
col_data$cancerType = dplyr::case_when(grepl('Normal',col_data$Tissue) ~ 'Normal',
                                       grepl('Tumour',col_data$Tissue) ~ 'PTC',TRUE ~ 'others')
colnames(col_data)[colnames(col_data) == 'Sex'] = 'sex'
colnames(col_data)[colnames(col_data) == 'Age'] = 'age'
col_data$age[grepl('fetus_',col_data$age)] = 'foetus'
col_data[,c('sampleID','source','sampleName','cancerType','age','sex')]


# Save this as a SummarizedExperiment object
bulk_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts_raw = bulk_raw_counts[,!colnames(bulk_raw_counts) %in% c("ensemblID", "geneSymbol", "geneLength","effLength")],
                counts_tpm = bulk_raw_TPM[,!colnames(bulk_raw_TPM) %in% c("ensemblID", "geneSymbol", "geneLength","effLength")]),
  rowData = row_data,
  colData = col_data,
  metadata = list(genomeVersion = 'GRCh38_v32_2020A-full',
                  gtf_filepath = '/nfs/cellgeni/STAR/human/2020A-full/GRCh38_v32_modified.gtf')
)
saveRDS(bulk_se,'Data/inhouse_bulk/inhouse_bulkRNA_fetalThyroid_paedPTC.RDS')




##----------------------------------------------------------##
##   Import relevant scRNA datasets: foetal_thyrocytes    ####
##----------------------------------------------------------##

##----- Normal foetal 
## Import foetal 2n thyrocytes only
fThy = readRDS('Data/fThyrocytes_2n_atlas.RDS')

fThy$finalAnn = fThy$cluster
fThy$finalAnn[fThy$finalAnn == 'Thyrocytes' & fThy$celltype == 'thy_TH_processing'] = 'fTFC1'
fThy$finalAnn[fThy$finalAnn == 'Thyrocytes' & fThy$celltype == 'thy_Lumen-forming'] = 'fTFC2'
fThy$annot = fThy$finalAnn
fThy$cellID = rownames(fThy@meta.data)

##----- Normal adult
aThy_Mosteiro23 = readRDS('Data/published_scRNAseq/Mosteiro_etal_2023/Mosteiro_etal_2023.RDS')
aThy_Mosteiro23$annot[aThy_Mosteiro23$annot %in% c('aTFC1','aTFC2','aTFC4','aTFC5')] = 'aTFC1'
aThy_Mosteiro23$annot[aThy_Mosteiro23$annot %in% c('aTFC3')] = 'aTFC2'
aThy_Mosteiro23$dataset = 'aThy_Mosteiro23'
Idents(aThy_Mosteiro23) = aThy_Mosteiro23$annot

##----- PTC adult
source('R/utils/sc_utils.R')

adult_PTC_dataset_list = c('Pu21'='Data/published_scRNAseq/Pu_etal_2021/Pu_etal_2021.RDS',
                      'Lu23'='Data/published_scRNAseq/Lu_etal_2023/Lu_etal_2023.RDS',
                      'Peng21'='Data/published_scRNAseq/Peng_etal_2021/Peng_etal_2021_thyOnly.RDS')
adult_PTC_sratObj_list = list()
for(dataset in names(adult_PTC_dataset_list)){
  srat = readRDS(adult_PTC_dataset_list[[dataset]])
  if(dataset == 'Pu21'){
    ## Pu et al., 2021
    srat = subset(srat,subset = cellID %in% srat$cellID[srat$celltype %in% c('fTFC1','fTFC2','thy_Lumen-forming') & 
                                                          srat$tissue_type %in% c('para-tumour','tumour')])
    srat$annot = paste0(srat$celltype,'-',srat$tissue_type)
    srat$annot = gsub('thy_Lumen-forming','fTFC2',srat$annot)
    srat$dataset = 'Pu_2021'
    
  }else if(dataset == 'Lu23'){
    ## Lu et al., 2023
    srat = subset(srat,subset = cellID %in% srat$cellID[srat$celltype %in% c('Epithelial cell','Malignant cell')])
    srat$annot = srat$celltype
    srat$annot[srat$annot == 'Epithelial cell'] = 'Thyrocytes'
    srat$annot[srat$annot == 'Malignant cell'] = 'Tumour'
    srat$dataset <- 'Lu_2023'
    
  }else if(dataset == 'Peng21'){
    ## Peng et al., 2021
    srat = subset(srat,subset = cellID %in% srat$cellID[srat$annot %in% c('Thyrocytes','Tumour')])
    srat$dataset = 'Peng_2021'
  }
  
  srat = standard_clustering(srat)
  DimPlot(srat,group.by = 'annot')
  
  adult_PTC_sratObj_list[[dataset]] = srat
}








##---------------------------------------------------##
##   1. Define the fTFC modules - based on DEGs    ####
##---------------------------------------------------##

## Assess the fraction of cells within each cell types (in the thyroid atlas) expressing each gene
fp = 'Results/2505/PTC_bulkRNAseq/fThyroid_pctExpressed_perCluster_perGene.csv'
if(file.exists(fp)){
  data_allCelltype = read.csv(fp)
}else{
  ## Import foetal 2n thyroid, with all cell types
  fThyroid = readRDS('Data/fThyroid_2n_atlas.RDS')
  fThyroid$finalAnn = fThyroid$cluster
  fThyroid$finalAnn[fThyroid$celltype == 'thy_TH_processing'] = 'fTFC1'
  fThyroid$finalAnn[fThyroid$celltype == 'thy_Lumen-forming'] = 'fTFC2'
  Idents(fThyroid) = fThyroid$finalAnn
  
  data_allCelltype = data.frame()
  for(clust in unique(fThyroid$finalAnn)){
    print(clust)
    nCell = rowSums(fThyroid@assays$RNA@counts[,fThyroid$cellID[fThyroid$finalAnn == clust]] > 0)
    pct = nCell/length(fThyroid$cellID[fThyroid$finalAnn == clust])
    tmp = data.frame(geneSym = names(nCell),nCell=nCell,pct=pct)
    colnames(tmp)[colnames(tmp) != 'geneSym'] = paste0(clust,'_',colnames(tmp)[colnames(tmp) != 'geneSym'])
    if(nrow(data_allCelltype) == 0){
      data_allCelltype = tmp
    }else{
      data_allCelltype = cbind(data_allCelltype,tmp[match(data_allCelltype$geneSym,tmp$geneSym),!colnames(tmp) %in% colnames(data_allCelltype)])
    }
  }
  
  
  write.csv(data_allCelltype,fp)
  
  
}





##----------------------------------------------------------------------##
##   1.-- Define the fetal thyrocyte cell states modules                ####      
##   based on DEGs + percentage_cellExpressed in other cell types       ##
##----------------------------------------------------------------------##

## Import the list of DEGs between 2n fTFC1 and fTFC2
deg = read.csv('~/lustre_mt22/Thyroid/Results_v2/2.1_DEG_2n_fTFC1.vs.fTFC2/oct24/celltype_donorID_allAgeGroup/allGenes_log2FC_allAgeGroup.csv')
deg = deg[deg$FDR < 0.05,]

## add pct.expressed in other cell types
dd = data_allCelltype[data_allCelltype$geneSym %in% deg$geneSym,!grepl('fTFC|Thyrocytes',colnames(data_allCelltype))]
dd2 = do.call(rbind,lapply(1:nrow(dd),function(i){
  ct = names(dd[i,grepl('_pct$',colnames(dd))])[dd[i,grepl('_pct$',colnames(dd))] == max(dd[i,grepl('_pct$',colnames(dd))])]
  pct = dd[i,ct]
  tmp = data.frame(ct = paste(ct,collapse = ':'),pct=unique(as.numeric(pct)),gene=dd$geneSym[i])
  if(ncol(tmp)!=3){
    print(i)
    print(dim(tmp))  
  }
  
  return(tmp)}))


deg = cbind(deg[deg$FDR < 0.05 & deg$geneSym %in% dd2$gene,],
            dd2[match(deg[deg$FDR < 0.05 & deg$geneSym %in% dd2$gene,]$geneSym,dd2$gene),])
# Only keep DEGs which are NOT highly expressed in other cell types (<30% of cells with expression)
deg = deg[deg$pct < 0.5,]
deg = deg[order(abs(deg$pct_diff),decreasing = T),]

# Further prioritise top DEGs by setting additional criteria:
# 1. abs(logFC) > 0.5
# 2. expressed in > 20% cells
# 3. top 100 genes for each cell state, order by pct_cellExpressed_difference between two cell states (maximise cell-state specificity)
deg2 = rbind(deg[abs(deg$logFC) > 0.5 & deg$pct_fTFC2 >= 30 & deg$direction == 'fTFC2_up',][1:100,],
            deg[abs(deg$logFC) > 0.5 & deg$pct_fTFC1 >= 30 & deg$direction == 'fTFC2_down',][1:100,])
table(deg2$direction)
quantile(abs(deg2$pct_diff),na.rm=T)
deg = deg2
##---> This is the final fTFC1/2 gene module
## Write this module
geneModule = deg[,c('ensID','geneSym','chr','logFC','logCPM','F','PValue','FDR','pct_fTFC1', 'pct_fTFC2','direction')]
geneModule$module = ifelse(geneModule$direction == 'fTFC2_down','fTFC1 signature','fTFC2 signature')
geneModule = geneModule[order(abs(geneModule$logFC),decreasing = T),]
geneModule = rbind(geneModule[geneModule$direction == 'fTFC2_up',],
                   geneModule[geneModule$direction == 'fTFC2_down',])
colnames(geneModule) = c('Ensembl_ID','Gene_symbol','Chromosome','log2FC fTFC2-vs-fTFC1','logCPM','F','PValue','FDR','percentage cells expressed - fTFC1','percentage cells expressed - fTFC2','DE direction','Gene signature')
write.csv(geneModule,'SupplementaryTables/SupplementaryTableS8_fTFC1.2_geneSignatures.csv',row.names = F)


deg = read.csv('SupplementaryTables/SupplementaryTableS8_fTFC1.2_geneSignatures.csv')
colnames(deg) = c('ensID','geneSym','chr','logFC','logCPM','F','PValue','FDR','pct_fTFC1', 'pct_fTFC2','direction','module')
write.csv(deg,file.path(outDir,'fTFC1.2_top100_geneSignatures.csv'))

##--------------------------------------------------------##
##    1.--Plot expression of the gene modules (DotPlot) ####
##--------------------------------------------------------##

genes_toPlot = c(deg$geneSym[deg$direction == 'fTFC2_down'],
                 deg$geneSym[deg$direction == 'fTFC2_up'])

Idents(fThy) = fThy$annot
DotPlot(fThy,
        features = genes_toPlot
)+
  RotatedAxis()+
  scale_y_discrete(position = "left")+
  scale_x_discrete(position = "bottom")+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=7,angle = 90,vjust = 0.5,hjust = 1,face='italic'),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = 'top') + xlab('') + ylab('')




##----------------------------------------------##
## 2. Module score for fTFC1/fTFC2 signatures ####
##----------------------------------------------##

## 2.-- Module score in single-cell RNAseq data -----

library(UCell)
## Define gene module
geneList = list('fTFC1' = deg$geneSym[deg$direction == 'fTFC2_down'],
                'fTFC2' = deg$geneSym[deg$direction == 'fTFC2_up'],
                'fTFC2_combined' = c(paste0(deg$geneSym[deg$direction == 'fTFC2_down'],'-'),
                                     paste0(deg$geneSym[deg$direction == 'fTFC2_up'],'+')),
                'fTFC1_combined' = c(paste0(deg$geneSym[deg$direction == 'fTFC2_down'],'+'),
                                     paste0(deg$geneSym[deg$direction == 'fTFC2_up'],'_')))


##----- Normal foetal
fThy <- UCell::AddModuleScore_UCell(fThy, features = geneList,ncores = 3)
fThy$dataset = 'fThy'
FeaturePlot(fThy,'fTFC2_UCell') 
DimPlot(fThyrocytes,group.by = 'celltype')

##----- Normal foetal thyroid
# fThyroid <- UCell::AddModuleScore_UCell(fThyroid, features = geneList,ncores = 1)
# fThyroid$dataset = 'fThyroid'
# FeaturePlot(fThyroid,'fTFC2_UCell') 
# write.csv(fThyroid@meta.data,'~/lustre_mt22/Thyroid/Results_v2/06_fThyrocytes_module_in_pPTC/UCell_fTFC1.2_signature_inFThyroid_2412.csv')

##----- Normal adult (Mosteiro)
aThy_Mosteiro23 <- UCell::AddModuleScore_UCell(aThy_Mosteiro23, features = geneList,ncores = 10)
FeaturePlot(aThy_Mosteiro23,'fTFC2_UCell') 

##----- Normal + PTC adult
adult_PTC_sratObj_list = lapply(adult_PTC_sratObj_list,function(x){
  srat = UCell::AddModuleScore_UCell(x, features = geneList,ncores = 10)
  return(srat)
})

# ##----- Normal + PTC paediatric
# pThy = readRDS('Results/2505/PTC_snRNAseq/03_pThyCancer_annotation/pPTC_clean_soupedXrhoLimNone_annotated_2505.RDS')
# pThy = subset(pThy,subset = annot %in% c('Thyrocytes','Tumour'))
# pThy = standard_clustering(pThy)
# pThy$annot = pThy$celltype
# pThy$cellID
# pThy$dataset <- 'pPTC_Sanger'
# 
# pThy <- UCell::AddModuleScore_UCell(pThy, features = geneList,ncores = 20)
# FeaturePlot(pThy,'fTFC2_combined_UCell')


##----- Aggregate data -------##
columns = c('cellID','annot','fTFC1_UCell','fTFC2_UCell','fTFC1_combined_UCell','fTFC2_combined_UCell','dataset')

## Adult PTC datasets
ucell_data_adult = do.call(rbind,lapply(adult_PTC_sratObj_list,function(x){
  x@meta.data[,columns]
}))

ucell_data_adult$annot_2 = gsub('fTFC1-|fTFC2-','',ucell_data_adult$annot)
ucell_data_adult$annot_2[ucell_data_adult$annot_2 %in% c('Thyrocytes')] = 'Normal'
ucell_data_adult$annot_2[ucell_data_adult$annot_2 %in% c('Malignant cell','tumour')] = 'Tumour'
ucell_data_adult$annot_2[ucell_data_adult$annot_2 %in% c('Epithelial cell','para-tumour')] = 'Normal'

write.csv(ucell_data_adult,file.path(outDir,'fTFC1.2_top100_geneSignatures_UCELL_scRNAseq_adultPTC.csv'))


## Normal fetal + adult datasets
ucell_data_normal = do.call(rbind,list('fThy'=fThy@meta.data[,columns],
                                       'aThy_Mosteiro23'=aThy_Mosteiro23@meta.data[,columns]
))
ucell_data_normal$annot[ucell_data_normal$annot %in% c('aTFC1')] = 'fTFC1'
ucell_data_normal$annot[ucell_data_normal$annot %in% c('aTFC2','thy_Lumen-forming')] = 'fTFC2'

write.csv(ucell_data_normal,file.path(outDir,'fTFC1.2_top100_geneSignatures_UCELL_scRNAseq_fetal.adult.Thy.csv'))

ucell_data_adult = read.csv(file.path(outDir,'fTFC1.2_top100_geneSignatures_UCELL_scRNAseq_adultPTC.csv'))
ucell_data_normal = read.csv(file.path(outDir,'fTFC1.2_top100_geneSignatures_UCELL_scRNAseq_fetal.adult.Thy.csv'))

## Combine them
ucell_data_adult$annot = ucell_data_adult$annot_2
ucell_data = rbind(ucell_data_normal,ucell_data_adult[,colnames(ucell_data_normal)])

## Do some plots
celltypes_toKeep = unique(ucell_data$annot[grepl('Thyrocyte|TFC1|TFC2|Tumour|Epithelial|Malignant|Normal',ucell_data$annot)])

ggplot(ucell_data[ucell_data$annot %in% celltypes_toKeep &
                    ucell_data$dataset %in% c('Lu_2023','Peng_2021','Pu_2021'),],aes(annot,fTFC1_UCell))+
  geom_hline(yintercept = 0,linetype=2,col=colAlpha(grey(0.4),0.4))+
  geom_boxplot(outlier.shape = NA,aes(fill = annot))+
  geom_quasirandom(size = 0.1,alpha = 0.1)+
  #scale_y_log10()+
  scale_fill_manual(values = c('Normal' = grey(0.8),'Tumour'=colAlpha('#511378',0.8)))+
  facet_grid(.~dataset,scales = 'free_x',space = 'free_x')+
  # theme_classic()+
  # theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))+
  theme_classic()+
  #ylim(-0.1,0.1)+
  #ggtitle(title)+
  xlab('')+ylab('fTFC2 signature score')+
  theme(panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
        strip.background=element_rect(linewidth=0),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1,colour = 'black'),legend.position = 'none')


ggplot(ucell_data[grepl('aTFC|Thyrocyte|Tum|Met|fTFC|Malig|Epi',ucell_data$annot),],aes(fTFC2_UCell,fTFC1_UCell,col=annot_2))+
  geom_point(size=0.2,alpha=0.2)+
  geom_hline(yintercept = 0.3)+
  geom_vline(xintercept = 0.2)+
  scale_color_manual(values = col25)+
  facet_grid(annot~dataset)+
  theme_classic()+
  theme(panel.border = element_rect(fill=F),axis.line = element_blank())


ucell_data = pivot_longer(ucell_data,cols = c('fTFC1_UCell','fTFC2_UCell','fTFC1_combined_UCell','fTFC2_combined_UCell'),names_to = 'module',values_to = 'score')

ggplot(ucell_data[ucell_data$annot %in% celltypes_toKeep,],aes(annot,score,fill=module))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = col25)+
  xlab('')+
  facet_grid(.~dataset,scales = 'free_x',space = 'free_x')+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        panel.border = element_rect(fill=F),axis.line = element_blank())


## Latest version in figures.R
fig4b_fTFC1.2_moduleScore = function(){
  ucell_data = read.csv(file.path(outDir,'fTFC1.2_top100_geneSignatures_UCELL_scRNAseq_fetal.adult.Thy.csv'),row.names = 1)
  
  columns = c('cellID','annot','fTFC1_UCell','fTFC2_UCell','fTFC1_combined_UCell','fTFC2_combined_UCell','dataset')
  
  dd = ucell_data[ucell_data$dataset %in% c('fThy','aThy_Mosteiro23') & ucell_data$annot %in% c('fTFC1','fTFC2','aTFC1','aTFC2'),]
  dd$annot[dd$annot %in% c('aTFC1','fTFC1')] = 'TFC1'
  dd$annot[dd$annot %in% c('aTFC2','fTFC2')] = 'TFC2'
  dd$dataset = factor(dd$dataset,c('fThy','aThy_Mosteiro23'))
  dd$annot = factor(dd$annot,c('TFC2','TFC1'))
  
  plotFun_fTFC1.2_moduleScore_fThy.aThy = function(noFrame=FALSE,noPlot=FALSE){
    if(noPlot & !noFrame){
      dd.sub = dd[sample(1:nrow(dd),100),]
      p1 = ggplot(dd.sub,aes(fTFC1_UCell,fTFC2_UCell,col=annot))+
        geom_point(alpha=0.5,size=0.01)+
        scale_color_manual(values = c('orange','black'))+
        facet_grid(dataset~.)+
        theme_classic(base_size = 11)+
        theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 0.8),axis.line = element_blank(),
              axis.text = element_text(colour = 'black',size=9.5),
              axis.ticks = element_line(colour = 'black'),
              strip.background = element_blank()) +
        xlim(min(dd$fTFC1_UCell),max(dd$fTFC1_UCell)) + ylim(min(dd$fTFC2_UCell),max(dd$fTFC2_UCell))+
        xlab('fTFC1 module score') + ylab('fTFC2 module score')
    }
    
    
    if(!noPlot){
      p1 = ggplot(dd,aes(fTFC1_UCell,fTFC2_UCell,col=annot))+
        geom_point(alpha=0.5,size=0.2)+
        scale_color_manual(values = c('orange','black'))+
        facet_grid(dataset~.)+
        theme_classic(base_size = 11)+
        theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 0.8),axis.line = element_blank(),
              axis.text = element_text(colour = 'black',size=9.5),
              axis.ticks = element_line(colour = 'black'),
              strip.background = element_blank()) +
        xlab('fTFC1 module score') + ylab('fTFC2 module score')
    }
    
    print(p1)
    
    # ## Plot with marginal density curves
    # library(ggplot2)
    # library(ggExtra)
    # library(patchwork)
    # 
    # # Split data by dataset
    # split_data <- split(dd, dd$dataset)
    # 
    # # Create ggMarginal plots for each subset
    # plot_list <- lapply(names(split_data), function(ds_name) {
    #   data_subset <- split_data[[ds_name]]
    #   
    #   # Base ggplot
    #   p <- ggplot(data_subset, aes(fTFC1_combined_UCell, fTFC2_combined_UCell, col = annot)) +
    #     geom_point(alpha = 0.5, size = 0.2) +
    #     scale_color_manual(values = c('orange', 'black')) +
    #     theme_classic(base_size = 11) +
    #     theme(
    #       panel.border = element_rect(fill = NA, colour = 'black', linewidth = 0.8),
    #       axis.line = element_blank(),
    #       axis.text = element_text(colour = 'black', size = 9.5),
    #       axis.ticks = element_line(colour = 'black'),
    #       strip.background = element_blank()
    #     ) +
    #     xlab('fTFC1 module score') +
    #     ylab('fTFC2 module score') +
    #     ggtitle(ds_name)
    #   
    #   # Add marginal densities
    #   ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)
    # })
    # 
    # # Step 3: Combine all plots into a grid
    # wrap_plots(plotlist = plot_list, ncol = 2)
    
  }
  
  saveFig(file.path(plotDir,'Fig4e_fTFC1.2_moduleScore_fThy_aThy'),plotFun_fTFC1.2_moduleScore_fThy.aThy,rawData=dd,width = 4.5,height = 5.6,res = 500)
  
  
  
  plotFun_fTFC1.2_moduleScoreDifference_fThy.aThy = function(noFrame=FALSE,noPlot=FALSE){
    library(ggbeeswarm)
    dd$annot = factor(dd$annot,c('TFC1','TFC2'))
    if(noPlot & !noFrame){
      dd.sub = dd[sample(1:nrow(dd),100),]
      p1 = ggplot(dd.sub,aes(annot,fTFC2_combined_UCell,fill=annot))+
        geom_quasirandom(width = 0.2,size=0.25,alpha=0.2,col=grey(0.7))+
        geom_boxplot(outlier.size = 0.01,alpha=0.7,width=0.24)+
        scale_fill_manual(values = c('TFC2'='orange','TFC1'=grey(0.1)))+
        facet_grid(dataset~.)+
        theme_classic(base_size = 11)+
        theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 0.8),axis.line = element_blank(),
              axis.text = element_text(colour = 'black',size=9.5),
              axis.ticks = element_line(colour = 'black'),
              strip.background = element_blank()) +
        ylim(min(dd$fTFC2_combined_UCell),max(dd$fTFC2_combined_UCell))+
        xlab('') + ylab('fTFC2 combined module score')
    }
    
    
    if(!noPlot){
      dd$score = dd$fTFC2_UCell / dd$fTFC1_UCell
      p1 = ggplot(dd,aes(annot,score,fill=annot))+
        geom_quasirandom(width = 0.2,size=0.25,alpha=0.2,col=grey(0.7))+
        geom_boxplot(outlier.size = 0.01,alpha=0.7,width=0.24)+
        scale_fill_manual(values = c('TFC2'='orange','TFC1'=grey(0.1)))+
        #geom_hline(yintercept = 0,linetype='dashed')+
        facet_grid(dataset~.)+
        theme_classic(base_size = 11)+
        theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 0.8),axis.line = element_blank(),
              axis.text = element_text(colour = 'black',size=9.5),
              axis.ticks = element_line(colour = 'black'),
              strip.background = element_blank()) +
        xlab('') + ylab('fTFC2 combined module score')
    }
    
    print(p1)
  }
  
  saveFig(file.path(plotDir,'Fig4e_fTFC1.2_combinedModuleScore_fThy_aThy'),plotFun_fTFC1.2_moduleScoreDifference_fThy.aThy,rawData=dd,width = 3,height = 5.6,res = 500)
  
}













## 2.-- Module score in bulk RNAseq data -----


source('R/helperFunctions.R')

##--- Define the gene module (convert geneSymbol --> ensID)
moduleList = lapply(geneList[names(geneList) %in% c('fTFC1','fTFC2')],function(i){
  o = gene_map$gene_id[match(i,gene_map$gene_name)]
  o = o[!is.na(o)]
  return(o)})

moduleList[['fTFC2_combined']] = list('up' = moduleList[['fTFC2']],
                                      'down' = moduleList[['fTFC1']])

moduleList[['fTFC1_combined']] = list('down' = moduleList[['fTFC2']],
                                      'up' = moduleList[['fTFC1']])

##--- import bulk counts and calculate cpmCnt in xx01_moduleScoring.R
#bulkRNA = import_bulkRNA_thyroid(bulk_sources = c('TCGA_Thyroid','He2021','inhouse'))
bulkRNA = import_bulkRNA_thyroid(bulk_sources = c('inhouse'='Data/inhouse_bulk/inhouse_bulkRNA_fetalThyroid_paedPTC.RDS',
                                                  'TCGA_Thyroid'='Data/published_bulkRNAseq/TCGA_Thyroid/TCGA_Thyroid_bulkRNA_se.RDS',
                                                  'REBC_THYR' = 'Data/published_bulkRNAseq/REBC_THYR_se.RDS',
                                                  'He2021'='Data/published_bulkRNAseq/He_etal_21/aPTC_He_2021_se.RDS',
                                                  'Lee2024' = 'Data/published_bulkRNAseq/Lee_etal_24/aPTC_Lee_2024_se.RDS'),gene_map = gene_map)
bulk_samples = bulkRNA[['bulk_samples']]
cpmCnt = bulkRNA[['cpmCnt']]
tpmCnt = bulkRNA[['tpm_count']]
rawCnt = bulkRNA[['raw_count']]

##---  Score the modules -----##

mtx = tpmCnt[,!colnames(tpmCnt) %in% c('ensID','geneLength')]
# apply the rankGenes method
bulk_ranked = rankGenes(mtx)


# apply the scoring function
allScore = data.frame()
for(i in 1:length(moduleList)){
  if(length(moduleList[[i]]) == 2){
    moduleScores = simpleScore(bulk_ranked,
                               upSet = moduleList[[i]][['up']],
                               downSet = moduleList[[i]][['down']])
    moduleScores = moduleScores[,c('TotalScore', 'TotalDispersion')]
  }else{
    moduleScores = simpleScore(bulk_ranked,upSet = moduleList[[i]])
  }
  
  # create a dataframe with the data required: scores and sample group
  scoredf = merge(bulk_samples,moduleScores,by.x=0,by.y=0)
  scoredf$moduleType = paste0(names(moduleList)[i])
  
  ## Add to allScore
  allScore = rbind(allScore,scoredf)
}


table(allScore$moduleType)
write.csv(allScore,file.path(outDir,'fTFC1.2_top100_geneSignatures_SingScore_bulkRNAseq.csv'))
allScore = read.csv(file.path(outDir,'fTFC1.2_top100_geneSignatures_SingScore_bulkRNAseq.csv'))



#plotDir = '~/lustre_mt22/Thyroid/Figures/manuscriptDraft_1124/Plots/'
#source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')


id_col = 'Sample'
title = 'foetalSig module'

fig4c_fThy_moduleScore_inBulkSamples = function(){
  library(ggbeeswarm)
  
  # create a dataframe with the data required: scores and sample group
  allScore$cancerType[grepl('FFPE',allScore$sampleName)] = paste0('FFPE_',allScore$cancerType[grepl('FFPE',allScore$sampleName)])
  allScore$cancerType[allScore$source == 'Sanger' & allScore$sampleName %in% c('PR66788b','PR66789b','PR66790b','PR66791b','PR66792b')] = 'Normal_foetal'
  allScore$sampleCol = ifelse(grepl('FFPE',allScore$sampleName),'FFPE','normal')
  allScore$sampleCol[allScore$source2 == 'scaThy'] = gsub(':.*$','',allScore$sampleID[allScore$source2 == 'scaThy'])
  allScore$source[allScore$source == 'Yoo_2021'] = 'Yoo_2016'
  allScore$source = factor(allScore$source,c('aPTC_Pu21','aPTC_Wang22','aThy_Hong23','aThy_Mosteiro23',
                                             'GTEx_Thyroid','TCGA_Thyroid','Yoo_2016','He_2021','Lee_2024',
                                             'stJudes_Thyroid','Lee_2021','Sanger','scRNAseq_fThy','snRNAseq_Y24.Y46'))
  
  allScore$group_facet_hor = allScore$source
  allScore$group_facet_ver = allScore$moduleType
  allScore$group_fill = allScore$cancerType
  
  allScore$group_facet_ver = allScore$moduleType
  dd = allScore[allScore$age == 'foetus' & allScore$source == 'Sanger' & allScore$moduleType %in% c('fTFC1','fTFC2'),]
  
  plotFun_sc.fThy.moduleScore_in_Sanger.Fetal.BulkSamples = function(noFrame=FALSE,noPlot=FALSE){
    
    p1 = ggplot(dd, aes(moduleType, TotalScore)) +
      geom_boxplot(aes(fill=group_fill),outlier.colour = 'white',position = 'dodge', alpha = 0.7,width=0.4,linewidth=0.3,fill=grey(0.7)) +
      geom_quasirandom(size=0.4,width = 0.15,alpha=0.6)+
      scale_y_continuous(limits = c(0,0.32))+
      theme_classic()+
      #ggtitle(title)+
      xlab('')+ylab('Module score')+
      theme(panel.border = element_rect(fill=F),axis.line = element_blank(),
            strip.background=element_rect(linewidth=0),
            axis.text = element_text(colour = 'black'),
            axis.ticks = element_line(colour = 'black'),
            axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1))
    
    print(p1)
  }
  
  saveFig(file.path(plotDir,'Fig4b_fTFC1.2_moduleScore_bulk.Foetal.Samples'),plotFun_sc.fThy.moduleScore_in_Sanger.Fetal.BulkSamples,rawData=dd,width = 1.6,height = 4,res = 500)
  

  
  moduleType_toUse = c('fTFC1','fTFC2') # not fTFC2_combined
  plotFun_fTFC_moduleScore = function(noFrame=FALSE,noPlot=FALSE){
    
    allScore$group_facet_ver = allScore$moduleType
    
    dd = allScore[grepl('FFPE|Thyrocytes|Tumour|fTFC|fThy|Metastatic|Normal|aTFC|C\\d|PTC',allScore$cancerType) & 
                    !grepl('follicular|tallCell|Metastatic|Primary',allScore$cancerType) &
                    !grepl('Primary',allScore$cancerType_details) &
                    allScore$source %in% c('Sanger','TCGA_Thyroid',#'Yoo_2016',
                                           'He_2021','Lee_2024') &
                    allScore$moduleType %in% moduleType_toUse,]
                    #allScore$moduleType %in% c('fTFC1','fTFC2'),]
    
    dd$ageGroup = ifelse(dd$source == 'Sanger',dd$ageCat,'adult')
    dd = dd[dd$ageGroup != 'foetus',]
    dd$cancerNormal = ifelse(dd$cancerType %in% c('Normal'),'Normal',
                             ifelse(dd$cancerType == 'Normal.adj','Normal.adj','Tumour'))
    dd$cancerNormal = factor(dd$cancerNormal,c('Normal','Normal.adj','Tumour'))
    dd$med_normal = NA
    for(dataset in unique(dd$source)){
      for(mod in unique(dd$moduleType)){
        tmp = dd[dd$source == dataset & dd$moduleType == mod,]
        med_normal = median(tmp$TotalScore[tmp$cancerNormal == 'Normal'])
        dd$med_normal[dd$source == dataset & dd$moduleType == mod] = med_normal
      }
    }
    
    dd$normalised_score = dd$TotalScore - dd$med_normal
    dd$source = factor(dd$source,c('Sanger','TCGA_Thyroid','He_2021','Lee_2024'))
    
    table(dd$cancerNormal,dd$cancerNormal,dd$source)
    
    
    p1 = ggplot(dd[!dd$sampleName %in% sample_metadata$geo_accession[grepl('_P|-P',sample_metadata$title)],], aes(cancerNormal, normalised_score)) +
      geom_hline(yintercept = 0,linetype=2,linewidth=0.3)+
      geom_quasirandom(size=0.4,width = 0.15,alpha=0.4)+
      geom_boxplot(aes(fill=cancerNormal),outlier.colour = 'white',position = 'dodge', alpha = 0.8,width=0.5,linewidth=0.3,colour='black') +
      #geom_point(data=dd[dd$sampleName %in% sample_metadata$geo_accession[grepl('_P|-P',sample_metadata$title)],],col='red',size=2)+
      scale_fill_manual(values =c('Normal'=grey(0.8),'Normal.adj'=grey(0.4),'Tumour'='#511378'))+
      #scale_fill_manual(values =c(col25,pal34H))+
      #scale_color_manual(values =c(col25,pal34H))+
      
      #scale_fill_manual(values = c(rep(col25[4],2),'#c7065a',col25[4],rep(colAlpha(col25[1],0.4),3),rep(grey(0.7),5))) +
      #scale_fill_manual(values = c(rep(col25[4],2),rep(grey(0.7),7))) +
      facet_grid(group_facet_ver~source,scales = 'free',space = 'free_x')+
      theme_classic()+
      #ylim(-0.1,0.1)+
      #ggtitle(title)+
      xlab('')+ylab('Centralised fTFC signature score')+
      theme(panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
            strip.background=element_rect(linewidth=0),
            axis.text = element_text(colour = 'black'),
            axis.ticks = element_line(colour = 'black'),
            axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1,colour = 'black'))
    
    print(p1)
  }
  
  saveFig(file.path(plotDir,'Fig4c_fTFC1.2_moduleScore_bulkSamples_sub'),plotFun_fTFC_moduleScore,rawData=allScore,width = 6,height = 7,res = 500)
  
}








## 2.-- Module score in Microarray data -----
moduleScore_microarray = read.csv('Results/2505/PTC_bulkRNAseq/published_pPTC_microarray/SingScore_microarray_GSE35570_fTFC1.2_signature.csv',row.names = 1)
# Remove PTC samples from individuals with age > 16 years old
moduleScore_microarray$donorID = gsub('normal thyroid-|PTC-radiation exposed-|PTC-radiation not exposed-','',moduleScore_microarray$Title)
moduleScore_microarray = moduleScore_microarray[!moduleScore_microarray$donorID %in% moduleScore_microarray$donorID[!is.na(moduleScore_microarray$Age) & moduleScore_microarray$Age > 16],]
moduleScore_microarray$dataset = 'GSE35570'
norm_medScore = median(moduleScore_microarray$TotalScore[moduleScore_microarray$Radiation == 'normal thyroid' & moduleScore_microarray$moduleType == 'fTFC2'])
moduleScore_microarray$score = moduleScore_microarray$TotalScore - norm_medScore
ggplot(moduleScore_microarray[moduleScore_microarray$moduleType == 'fTFC2',],aes(Radiation,score))+
  geom_boxplot(aes(fill=RET_PTC),outlier.colour = 'white',position = 'dodge', alpha = 1,width=0.5,linewidth=0.3,colour='black')+
  
  facet_wrap(vars(moduleType),nrow=1)+
  theme_classic(base_size = 14)+
  xlab('')+ylab('fTFC2 signature score')+
  geom_hline(yintercept = 0)+
  #scale_fill_manual(values =c('Normal'=grey(0.8),'Tumour'='#511378'))+
  theme(panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
        strip.background=element_rect(linewidth=0),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1,colour = 'black'))



a = pivot_wider(moduleScore_microarray,id_cols = c('RET_PTC','Radiation','Title'),names_from = 'moduleType',values_from = 'TotalScore')
a$fTFC2_minus_fTFC1 = a$fTFC2 - a$fTFC1
ggplot(a,aes(RET_PTC,fTFC2_minus_fTFC1))+
  geom_boxplot(aes(fill=RET_PTC))+
  facet_wrap(vars(Radiation),nrow=1)+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))



## Some tweaks to sample metadata
bulk_samples$celltype = ifelse(bulk_samples$source2 %in% c('scaThy'),gsub('^.*:','',bulk_samples$sampleID),
                               ifelse(bulk_samples$source2 %in% c('snRNAseq_Y24.Y46'),bulk_samples$sampleID,bulk_samples$source))
bulk_samples$cancerType[bulk_samples$source2 == 'scaThy'] = bulk_samples$celltype[bulk_samples$source2 == 'scaThy']
bulk_samples$cancerType[bulk_samples$source2 == 'snRNAseq_Y24.Y46'] = bulk_samples$celltype[bulk_samples$source2 == 'snRNAseq_Y24.Y46']
bulk_samples$sampleGroup = ifelse(bulk_samples$cancerType %in% c('Normal','aTFC1','aTFC2','C0','C1','C2','C3','C4','fTFC1','fTFC2','fThy_cycling','Thyrocytes','Thyrocytes:Y24','Thyrocytes:Y46'),'Normal',
                                  ifelse(bulk_samples$source2 %in% c('scaThy','snRNAseq_Y24.Y46') & bulk_samples$cancerType %in% c('Tumour','Tumour:Y24','Tumour:Y46','Tumour_BRAF_V600E', 'Tumour_RET_FARP1_fusion','Met:Y46'),'Tumour',
                                         ifelse(bulk_samples$source2 %in% c('scaThy','snRNAseq_Y24.Y46'),'other_celltypes','Tumour')))
bulk_samples$sampleGroup[bulk_samples$source == 'scRNAseq_fThy'] = bulk_samples$cancerType[bulk_samples$source == 'scRNAseq_fThy']
bulk_samples$sampleGroup[bulk_samples$source == 'Sanger' & bulk_samples$sampleName %in% c('PR66788b','PR66789b','PR66790b','PR66791b','PR66792b')] = 'Normal_foetal'


##--- Plot a heatmap of the genes in the modules
fig4b_heatmap_bulkRNA = function(){
  col_fun = circlize::colorRamp2(c(-2, 0, 3), c("white", "#f7e8e3", "#bf363a"))  
  library(ComplexHeatmap)
  cpmCnt.sub = cpmCnt[rownames(cpmCnt) %in% do.call(c,moduleList),]
  rownames(cpmCnt.sub) = geneMap$geneSym[match(rownames(cpmCnt.sub),geneMap$ensID)]
  #cpmCnt.sub = cpmCnt.sub[,colnames(cpmCnt.sub) %in% bulk_samples$sampleID[bulk_samples$source %in% c('GTEx_Thyroid','TCGA_Thyroid','Yoo_2021','stJudes_Thyroid','He_2021','Lee_2021')]]
  #cpmCnt.sub = cpmCnt.sub[,colnames(cpmCnt.sub) %in% bulk_samples$sampleID[!(bulk_samples$source2 %in% c('snRNAseq_Y24.Y46','scaThy') & bulk_samples$sampleGroup == 'other_celltypes') ]]
  cpmCnt.sub = cpmCnt.sub[,colnames(cpmCnt.sub) %in% bulk_samples$sampleID[!(bulk_samples$source2 %in% c('snRNAseq_Y24.Y46','scaThy') & bulk_samples$sampleGroup %in% c('other_celltypes')) & 
                                                                             bulk_samples$source %in% c('GTEx_Thyroid',#'TCGA_Thyroid','Yoo_2021',
                                                                                                        'stJudes_Thyroid','scRNAseq_fThy','Sanger') & !bulk_samples$cancerType %in% c('FA','others','fThy_cycling')]]
  #& bulk_samples$source %in% c('GTEx_Thyroid','TCGA_Thyroid','Yoo_2021','stJudes_Thyroid','He_2021','Lee_2021')
  #colnames(cpmCnt.sub) %in% bulk_samples$sampleID[bulk_samples$source %in% c('Lee_2021')]]
  
  # group = paste0(bulk_samples$cancerType[match(colnames(cpmCnt.sub),bulk_samples$sampleID)],
  #                bulk_samples$source[match(colnames(cpmCnt.sub),bulk_samples$sampleID)])
  dim(cpmCnt.sub)
  
  ## Plot heatmap
  mdat = bulk_samples[match(colnames(cpmCnt.sub),bulk_samples$sampleID),]
  mdat$colGroup = paste0(mdat$sampleGroup,':',mdat$source)
  mdat$colGroup = factor(mdat$colGroup,c('fTFC1:scRNAseq_fThy','fTFC2:scRNAseq_fThy','fThy_cycling:scRNAseq_fThy','Normal:GTEx_Thyroid','Normal:TCGA_Thyroid','Normal:Yoo_2021',
                                         'Tumour:TCGA_Thyroid','Tumour:Yoo_2021','Tumour:stJudes_Thyroid','Normal_foetal:Sanger','Normal:Sanger','Tumour:Sanger'))
  mdat$colGroup = factor(mdat$colGroup,c('fTFC1:scRNAseq_fThy','fTFC2:scRNAseq_fThy','fThy_cycling:scRNAseq_fThy','Normal:GTEx_Thyroid','Normal:TCGA_Thyroid','Tumour:TCGA_Thyroid',
                                         'Normal:Yoo_2021','Tumour:Yoo_2021','Tumour:stJudes_Thyroid','Normal_foetal:Sanger','Normal:Sanger','Tumour:Sanger'))
  table(mdat$sampleGroup,mdat$source)
  
  
  
  
  
  plotFun_bRNA.only = function(noFrame=FALSE,noPlot=FALSE){
    botAnno = HeatmapAnnotation('Cancer_normal' = mdat$sampleGroup,
                                col = list('Cancer_normal'=c('other_celltypes'=grey(0.7),'Tumour'='red','Normal' =grey(0.5),
                                                             'fThy_cycling' = grey(0.7),
                                                             'fTFC1'='#047d94','fTFC2'='#4CC7CF',
                                                             'Normal_foetal' = 'green')))
    
    col_fun = circlize::colorRamp2(c(-4, 0, 4), c("white",'#e6b5a8', "#bf363a"))
    
    hm = Heatmap(t(scale(t(scale(cpmCnt.sub)))),cluster_column_slices = F,
                 col = col_fun,
                 show_row_dend = F,show_column_dend = F,
                 row_names_gp = gpar(fontsize = 7),column_title_rot = 90,column_names_gp = gpar(fontsize=9),
                 show_row_names = F,show_column_names = F,
                 row_split = ifelse(rownames(cpmCnt.sub) %in% geneList[['fTFC1']],'fTFC1','fTFC2'),
                 column_split = mdat$colGroup,name = 'column then row z-scaled',
                 bottom_annotation = botAnno)
    hm
    col_fun = circlize::colorRamp2(c(0.1,0.9,1.3), c("black", "white", "#bf363a"))
    col_fun = circlize::colorRamp2(c(-2,0,2), c("black", "white", "#bf363a"))
    hm = Heatmap(t(scale(t(cpmCnt.sub),center = F,scale = T)),cluster_column_slices = F,
                 #col = col_fun,
                 show_row_dend = F,show_column_dend = F,
                 row_names_gp = gpar(fontsize = 7),column_title_rot = 90,column_names_gp = gpar(fontsize=9),
                 show_row_names = F,show_column_names = F,
                 row_split = ifelse(rownames(cpmCnt.sub) %in% geneList[['fTFC1']],'fTFC1','fTFC2'),
                 column_split = mdat$colGroup,name = 'column then row z-scaled',
                 bottom_annotation = botAnno)
    print(hm)
  }
  
  saveFig(file.path(plotDir,'Fig4b_heatmap_fThy.markers_bRNA'),plotFun_bRNA.only,rawData=cpmCnt.sub,width = 13,height = 5,res = 500)  
  
  
  plotFun_bRNA.only_v2 = function(noFrame=FALSE,noPlot=FALSE){
    rowAnno = rowAnnotation('Cancer_normal' = mdat$sampleGroup,
                            col = list('Cancer_normal'=c('other_celltypes'=grey(0.7),'Tumour'='red','Normal' = '#dec2b8',#grey(0.5),
                                                         'fThy_cycling' = grey(0.7),
                                                         'fTFC1'='#047d94','fTFC2'='#4CC7CF','Normal_foetal'='green')))
    
    
    
    
    # col_fun = circlize::colorRamp2(c(-4, 0, 4), c("black", "white", "#bf363a"))
    # 
    # hm = Heatmap((scale(t(scale(cpmCnt.sub)))),cluster_row_slices = F,
    #              col = col_fun,
    #              show_row_dend = F,show_column_dend = F,
    #              row_names_gp = gpar(fontsize = 7),row_title_rot = 0,column_names_gp = gpar(fontsize=9),
    #              row_title_side = 'right',
    #              show_row_names = F,show_column_names = F,
    #              column_split = ifelse(rownames(cpmCnt.sub) %in% geneList[['fTFC1']],'fTFC1','fTFC2'),
    #              row_split = mdat$colGroup,name = 'column then row z-scaled',
    #              left_annotation = rowAnno)
    
    
    
    col_fun = circlize::colorRamp2(c(-1, 0, 3), c("white", grey(0.9), 'black'))
    hm = Heatmap((scale(t(scale(cpmCnt.sub)),center = F)),cluster_row_slices = F,
                 col = col_fun,column_gap = unit(3,'mm'),
                 show_row_dend = F,show_column_dend = F,cluster_columns = T,
                 row_names_gp = gpar(fontsize = 7),row_title_rot = 0,column_names_gp = gpar(fontsize=7),
                 row_title_side = 'right',
                 show_row_names = F,show_column_names = F,
                 column_split = ifelse(rownames(cpmCnt.sub) %in% geneList[['fTFC1']],'fTFC1','fTFC2'),
                 row_split = mdat$colGroup,name = 'column then row z-scaled',
                 left_annotation = rowAnno)
    
    
    print(hm)
  }
  saveFig(file.path(plotDir,'Fig4b_heatmap_fThy.markers_bRNA_v2'),plotFun_bRNA.only_v2,rawData=cpmCnt.sub,width = 9,height = 9,res = 500)  
  
  
}


# 3.-- Additional dataset REBC-THYR -----------------
rebc_thyr_scores = read.delim('~/thyroid/data/REBC-THYR-bulk-TPMs.scores.tsv',sep = '\t')
rebc_mdat = read.delim('~/thyroid/data/REBC-THYR.metadata.tsv',sep='\t',header = F)
rebc_clinical = read.delim('~/FetalThyroidAtlas/Data/published_bulkRNAseq/Morton_21_TableS1_abg2538-data-s1.txt',sep='\t',header = T)
rebc_sample_sheet = read.delim('~/FetalThyroidAtlas/Data/published_bulkRNAseq/gdc_sample_sheet.2025-08-21.tsv',sep='\t')
rebc_sample_sheet$Case.ID = gsub('-YQ','',rebc_sample_sheet$Case.ID)

table(rebc_sample_sheet$Case.ID %in% rebc_clinical$REBC_ID)
table(rebc_clinical$REBC_ID %in% rebc_sample_sheet$Case.ID)
rebc_clinical$REBC_ID[!rebc_clinical$REBC_ID %in% rebc_sample_sheet$Case.ID]
rebc_sample_sheet$Case.ID[!rebc_sample_sheet$Case.ID %in% rebc_clinical$REBC_ID]

rebc_sample_sheet = left_join(rebc_sample_sheet,rebc_clinical,by=c('Case.ID' = 'REBC_ID'))
rebc_sample_sheet$file_name = gsub('\\.rna_seq.augmented_star_gene_counts.tsv$','',rebc_sample_sheet$File.Name)
table(rebc_thyr_scores$Row.names %in% rebc_sample_sheet$file_name)

colnames(rebc_mdat) = c('sample','tissue_type')
table(rebc_thyr_scores$Row.names %in% rebc_sample_sheet$file_name)
rebc_thyr_scores = left_join(rebc_thyr_scores,rebc_sample_sheet,by=c('Row.names'='file_name'))
rebc_thyr_scores$age_group = ifelse(rebc_thyr_scores$AGE_SURGERY <= 10,'<=10',
                                    ifelse(rebc_thyr_scores$AGE_SURGERY <= 16,'<=16','>16'))
rebc_thyr_scores$age_group = factor(rebc_thyr_scores$age_group,c('<=10','<=16','>16'))

ggplot(rebc_thyr_scores,aes(AGE_SURGERY,TotalScore,group=Tissue.Type,col=Tissue.Type))+
  geom_point(aes())+
  geom_smooth()+
  facet_wrap(vars(moduleType))

ggplot(rebc_thyr_scores,aes(Tissue.Type,TotalScore))+
  geom_boxplot(aes(fill=Tissue.Type))+
  geom_quasirandom(size=1,width = 0.2)+
  facet_grid(age_group ~ moduleType)


# 3.-- Additional dataset Junak. 2016 GSE35570 -----------------
junak16 = read.delim('~/thyroid/data/gse35570_scores.tsv')
junak16$group = ifelse(grepl('normal',junak16$Tumor_type),'normal',
                       ifelse(junak16$Age_operation <= 12,'<=10',
                              ifelse(junak16$Age_operation <= 16,'<=16','>16')))
junak16$group = factor(junak16$group,c('normal','<=10','<=16','>16'))
ggplot(junak16,aes(group,TotalScore))+
  geom_boxplot(aes(fill=group),outlier.shape = NA)+
  geom_quasirandom(aes(col=Tumor_type),size=1,width = 0.2)+
  facet_wrap(vars(moduleType),scales = 'free_y')

ggplot(junak16,aes(Age_operation,TotalScore))+
  geom_point(aes(co))+
  facet_wrap(vars(moduleType))


dom12 = read.delim('~/thyroid/data/gse33630_scores.tsv',sep = '\t')
dom12$donorID = paste0('UA',gsub('N$|T$','',dom12$Title))
dom12 = dom12[!grepl('ATC',dom12$Title),]

dom12_mdat = read_excel('~/FetalThyroidAtlas/Data/published_bulkRNAseq/Dom12_clinicaldata_bjc2012302x1.xls')
dom12_mdat = dom12_mdat[!is.na(dom12_mdat$`Exposed Paired samples`) & dom12_mdat$`Exposed Paired samples` != 'Non exposed Paired samples',]
dom12_mdat$donorID = gsub('N$|T$','',dom12_mdat$`Exposed Paired samples`)
table(dom12$donorID %in% dom12_mdat$donorID)

dom12 = left_join(dom12,dom12_mdat,by='donorID')
dom12$group = ifelse(grepl('N$',dom12$Title),'normal','tumour')
dom12$age_group = ifelse(dom12$`Age at operation` <=16,'<=16','>16')
dom12$age_group = factor(dom12$age_group,c('<=16','>16'))
ggplot(dom12,aes(age_group,TotalScore))+
  geom_boxplot(aes(fill=group),outlier.shape = NA)+
  geom_quasirandom(aes(col=Tumor_type),size=1,width = 0.2)+
  facet_wrap(vars(moduleType),scales = 'free_y')


ggplot(dom12,aes(`Age at operation`,TotalScore))+
  geom_point(aes(col=group),size=0.6)+
  facet_wrap(vars(moduleType))



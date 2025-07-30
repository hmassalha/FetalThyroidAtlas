##--- Scoring for the enrichment of fetal thyrocyte cell states (fTFC1 and fTFC2) signals in adult and paediatric PTC ---##
##    1. Define fetal cell state specific gene module
##    2. Score for the enrichment of these modules in adult + children PTC

outDir = "~/lustre_mt22/Thyroid/Results_v2/06_fThyrocytes_module_in_pPTC/oct24"
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

setwd(outDir)




##----------------##
##   Libraries  ####
##----------------##
library(tidyverse)
library(Seurat)
library(GenomicFeatures)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")


##----------------------------##
##   Set Global parameters  ####
##----------------------------##

#Define genomic coordinates
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)

## Generic gene map
geneMap = read.table('~/lustre_mt22/Data/thyroid_10X/cellranger710_count_46320_SB_Thy_R13236839_GRCh38-2020-A/filtered_feature_bc_matrix/features.tsv.gz',sep = '\t')
colnames(geneMap) = c('ensID','geneSym','library')
geneMap$geneSym[duplicated(geneMap$geneSym)] = paste0(geneMap$geneSym[duplicated(geneMap$geneSym)],'.1')
geneMap$geneSym = gsub('_','-',geneMap$geneSym)
geneMap$chr = as.character(seqnames(gns)[match(geneMap$ensID,gns$gene_id)])




##----------------------------------------------------------##
##   Import relevant scRNA datasets: foetal_thyrocytes    ####
##----------------------------------------------------------##

##----- Normal foetal 
## Import foetal 2n thyrocytes only
fThy = readRDS('/lustre/scratch125/casm/team274sb/mt22/Thyroid/Data/fetalThyroid/fThyrocytes_2n_jul23.RDS')

fThy$finalAnn = fThy$cluster
fThy$finalAnn[fThy$finalAnn == 'Thyrocytes' & fThy$celltype == 'thy_TH_processing'] = 'fTFC1'
fThy$finalAnn[fThy$finalAnn == 'Thyrocytes' & fThy$celltype == 'thy_Lumen-forming'] = 'fTFC2'
fThy$annot = fThy$finalAnn


##----- Normal adult
aThy_2 = readRDS('/lustre/scratch125/casm/team274sb/mt22/Thyroid/Data/published_scRNAseq/aThyroid_Mosteiro_2023/aThyroid_Mosteiro_2023_sratObj.RDS')
aThy_2$annot[aThy_2$annot %in% c('aTFC1','aTFC2','aTFC4','aTFC5')] = 'aTFC1'
aThy_2$annot[aThy_2$annot %in% c('aTFC3')] = 'aTFC2'
aThy_2$dataset = 'aThy_Mosteiro23'
Idents(aThy_2) = aThy_2$annot


# aPTC_Pu21 = readRDS('~/lustre_mt22/Thyroid/Data/published_scRNAseq/Pu_etal_2021/Pu21_annotated_sratObj.RDS')
# aPTC_Pu21$dataset = 'aPTC_Pu21'
# aPTC_Pu21$annot = aPTC_Pu21$celltype



##------------------------------------------------##
##   Define the fTFC modules - based on DEGs    ####
##------------------------------------------------##

## Assess the fraction of cells within each cell types (in the thyroid atlas) expressing each gene
fp = '~/lustre_mt22/Thyroid/Results_v2/06_fThyrocytes_module_in_pPTC/fThyroid_pctExpressed_perCluster_perGene.csv'
if(file.exists(fp)){
  data_allCelltype = read.csv(fp)
}else{
  ## Import foetal 2n thyroid, with all cell types
  fThyroid = readRDS('~/lustre_mt22/Thyroid/Data/fetalThyroid/fThyroid_annotated_fromHassan_jul23.RDS')
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
##   Define the fetal thyrocyte cell states modules                   ####      
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
write.csv(geneModule,'~/lustre_mt22/Thyroid/Figures/manuscriptDraft_1124/SupplementaryTableS8_fTFC1.2_geneSignatures.csv',row.names = F)

##------------------------------------------------------##
##    Plot expression of the gene modules (DotPlot)   ####
##------------------------------------------------------##

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




##-------------------------------------------##
## Module score for fTFC1/fTFC2 signatures ####
##-------------------------------------------##

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


##----- Normal adult
aThy_2 <- UCell::AddModuleScore_UCell(aThy_2, features = geneList,ncores = 3)
FeaturePlot(aThy_2,'fTFC1_UCell')
DimPlot(aThy_2,group.by = 'seurat_clusters')


# aPTC_Pu21 <- UCell::AddModuleScore_UCell(aPTC_Pu21, features = geneList,ncores = 3)
# FeaturePlot(aPTC_Pu21,'fTFC1_UCell')


##----- Aggregate data -------##
columns = c('cellID','annot','fTFC1_UCell','fTFC2_UCell','fTFC1_combined_UCell','fTFC2_combined_UCell','dataset')
ucell_data = do.call(rbind,list('fThy'=fThy@meta.data[,columns],
                                'aThy'=aThy_2@meta.data[,columns]
                                # 'aPTC_Pu21'=aPTC_Pu21@meta.data[,columns]
))
ucell_data$annot[ucell_data$annot %in% c('aTFC1')] = 'fTFC1'
ucell_data$annot[ucell_data$annot %in% c('aTFC2')] = 'fTFC2'

write.csv(ucell_data,'~/lustre_mt22/Thyroid/Results_v2/06_fThyrocytes_module_in_pPTC/UCell_fTFC1.2_signature_2412.csv')

## Do some plots
celltypes_toKeep = unique(ucell_data$annot[grepl('Thyrocyte|TFC1|TFC2|Tumour',ucell_data$annot)])
ggplot(ucell_data[ucell_data$annot %in% celltypes_toKeep,],aes(annot,fTFC2_combined_UCell))+
  geom_boxplot(outlier.shape = NA)+
  #scale_y_log10()+
  geom_hline(yintercept = 0)+
  facet_grid(.~dataset,scales = 'free_x',space = 'free_x')+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))

ggplot(ucell_data[grepl('aTFC|Thyrocyte|Tum|Met|fTFC',ucell_data$annot),],aes(fTFC2_UCell,fTFC1_UCell,col=annot))+
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
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))


## Latest version in figures.R
fig4b_fTFC1.2_moduleScore = function(){
  ucell_data = read.csv('~/lustre_mt22/Thyroid/Results_v2/06_fThyrocytes_module_in_pPTC/UCell_fTFC1.2_signature_2412.csv')
  
  columns = c('cellID','annot','fTFC1_UCell','fTFC2_UCell','fTFC1_combined_UCell','fTFC2_combined_UCell','dataset')
  
  dd = ucell_data[ucell_data$dataset %in% c('fThy','aThy_Mosteiro23') & ucell_data$annot %in% c('fTFC1','fTFC2','aTFC1','aTFC2'),]
  dd$annot[dd$annot %in% c('aTFC1','fTFC1')] = 'TFC1'
  dd$annot[dd$annot %in% c('aTFC2','fTFC2')] = 'TFC2'
  dd$dataset = factor(dd$dataset,c('fThy','aThy_Mosteiro23'))
  
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
  }
  
  saveFig(file.path(plotDir,'Fig4e_fTFC1.2_moduleScore_fThy_aThy'),plotFun_fTFC1.2_moduleScore_fThy.aThy,rawData=dd,width = 4.5,height = 5.6,res = 500,useDingbats = F)
  
  
  
  plotFun_fTFC1.2_moduleScoreDifference_fThy.aThy = function(noFrame=FALSE,noPlot=FALSE){
    library(ggbeeswarm)
    if(noPlot & !noFrame){
      dd.sub = dd[sample(1:nrow(dd),100),]
      p1 = ggplot(dd.sub,aes(annot,fTFC2_combined_UCell,fill=annot))+
        geom_quasirandom(width = 0.2,size=0.25,alpha=0.2,col=grey(0.7))+
        geom_boxplot(outlier.size = 0.01,alpha=0.7,width=0.24)+
        scale_fill_manual(values = c('orange',grey(0.1)))+
        facet_grid(dataset~.)+
        theme_classic(base_size = 11)+
        theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 0.8),axis.line = element_blank(),
              axis.text = element_text(colour = 'black',size=9.5),
              axis.ticks = element_line(colour = 'black'),
              strip.background = element_blank()) +
        ylim(min(dd$fTFC2_combined_UCell),max(dd$fTFC2_combined_UCell))+
        xlab('') + ylab('fTFC2 module score')
    }
    
    
    if(!noPlot){
      p1 = ggplot(dd,aes(annot,fTFC2_combined_UCell,fill=annot))+
        geom_quasirandom(width = 0.2,size=0.25,alpha=0.2,col=grey(0.7))+
        geom_boxplot(outlier.size = 0.01,alpha=0.7,width=0.24)+
        scale_fill_manual(values = c('orange',grey(0.1)))+
        #geom_hline(yintercept = 0,linetype='dashed')+
        facet_grid(dataset~.)+
        theme_classic(base_size = 11)+
        theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 0.8),axis.line = element_blank(),
              axis.text = element_text(colour = 'black',size=9.5),
              axis.ticks = element_line(colour = 'black'),
              strip.background = element_blank()) +
        xlab('') + ylab('fTFC2-to-fTFC1 module score')
    }
    
    print(p1)
  }
  
  saveFig(file.path(plotDir,'Fig4e_fTFC1.2_moduleScoreDifference_fThy_aThy'),plotFun_fTFC1.2_moduleScoreDifference_fThy.aThy,rawData=dd,width = 3,height = 5.6,res = 500,useDingbats = F)
  
}









##--------------------------------##
##      Score in bulk data      ####
##--------------------------------##
source('~/lustre_mt22/Thyroid/scripts/final_script/forPublication/v1/helperFunctions.R')

##--- Define the gene module (convert geneSymbol --> ensID)
moduleList = lapply(geneList[names(geneList) %in% c('fTFC1','fTFC2')],function(i){
  o = geneMap$ensID[match(i,geneMap$geneSym)]
  o = o[!is.na(o)]
  return(o)})

moduleList[['fTFC2_combined']] = list('up' = moduleList[['fTFC2']],
                                      'down' = moduleList[['fTFC1']])

moduleList[['fTFC1_combined']] = list('down' = moduleList[['fTFC2']],
                                      'up' = moduleList[['fTFC1']])

##--- import bulk counts and calculate cpmCnt in xx01_moduleScoring.R
bulkRNA = import_bulkRNA_thyroid(bulk_sources = c('TCGA_Thyroid','He2021','inhouse'))
bulk_samples = bulkRNA[['bulk_samples']]
cpmCnt = bulkRNA[['cpmCnt']]
tpmCnt = bulkRNA[['tpmCnt']]
rawCnt = bulkRNA[['rawCnt']]

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
write.csv(allScore,'~/lustre_mt22/Thyroid/Results_v2/06_fThyrocytes_module_in_pPTC/SingScore_bulkRNA_fTFC1.2_signature_2412.csv')




plotDir = '~/lustre_mt22/Thyroid/Figures/manuscriptDraft_1124/Plots/'
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')


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
                                             'GTEx_Thyroid','TCGA_Thyroid','Yoo_2016','He_2021',
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
      scale_y_continuous(breaks = c(0,0.1,0.2,0.3),labels = c(0.0,0.1,0.2,0.3),limits = c(0,0.2))+
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
  
  saveFig(file.path(plotDir,'Fig4b_fTFC1.2_moduleScore_bulk.Foetal.Samples'),plotFun_sc.fThy.moduleScore_in_Sanger.Fetal.BulkSamples,rawData=dd,width = 1.6,height = 4,res = 500,useDingbats = F)
  
  
  
  plotFun_fTFC2_combined_moduleScore = function(noFrame=FALSE,noPlot=FALSE){
    
    allScore$group_facet_ver = allScore$moduleType
    
    dd = allScore[grepl('FFPE|Thyrocytes|Tumour|fTFC|fThy|Metastatic|Normal|aTFC|C\\d|PTC',allScore$cancerType) & 
                    !grepl('follicular|tallCell|Metastatic|Primary',allScore$cancerType) &
                    !grepl('Primary',allScore$cancerType_details) &
                    allScore$source %in% c('Sanger','TCGA_Thyroid',#'Yoo_2016',
                                           'He_2021') &
                    allScore$moduleType == 'fTFC2_combined',]
    
    dd$ageGroup = ifelse(dd$source == 'Sanger',dd$ageCat,'adult')
    dd = dd[dd$ageGroup != 'foetus',]
    dd$cancerNormal = ifelse(dd$cancerType %in% c('Normal'),'Normal',
                             ifelse(dd$cancerType == 'Normal.adj','Normal.adj','Tumour'))
    dd$cancerNormal = factor(dd$cancerNormal,c('Normal','Normal.adj','Tumour'))
    dd$med_normal = NA
    for(dataset in unique(dd$source)){
      tmp = dd[dd$source == dataset,]
      med_normal = median(tmp$TotalScore[tmp$cancerNormal == 'Normal'])
      dd$med_normal[dd$source == dataset] = med_normal
    }
    
    dd$normalised_score = dd$TotalScore - dd$med_normal
    dd$source = factor(dd$source,c('Sanger','TCGA_Thyroid','He_2021'))
    
    table(dd$cancerNormal,dd$cancerNormal,dd$source)
    
    
    p1 = ggplot(dd, aes(cancerNormal, normalised_score)) +
      geom_boxplot(aes(fill=cancerNormal),outlier.colour = 'white',position = 'dodge', alpha = 0.7,width=0.5,linewidth=0.3,colour='black') +
      geom_quasirandom(size=0.4,width = 0.15,alpha=0.6)+
      scale_fill_manual(values =c(grey(0.8),grey(0.4),'#511378'))+
      #scale_fill_manual(values =c(col25,pal34H))+
      #scale_color_manual(values =c(col25,pal34H))+
      geom_hline(yintercept = 0)+
      #scale_fill_manual(values = c(rep(col25[4],2),'#c7065a',col25[4],rep(colAlpha(col25[1],0.4),3),rep(grey(0.7),5))) +
      #scale_fill_manual(values = c(rep(col25[4],2),rep(grey(0.7),7))) +
      facet_grid(group_facet_ver~source,scales = 'free',space = 'free_x')+
      theme_classic()+
      #ggtitle(title)+
      xlab('')+ylab('Centralised fTFC2 signature score')+
      theme(panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
            strip.background=element_rect(linewidth=0),
            axis.text = element_text(colour = 'black'),
            axis.ticks = element_line(colour = 'black'),
            axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1,colour = 'black'))
    
    print(p1)
  }
  saveFig(file.path(plotDir,'Fig3c_MLDS_moduleScore_bulkSamples'),plotFun_MLDS_moduleScore_inBulkSamples,rawData=allScore,width = 7,height = 5,res = 500,useDingbats = T)
  saveFig(file.path(plotDir,'Fig3c_fTFC1.2_moduleScore_bulkSamples_sub'),plotFun_fTFC2_combined_moduleScore,rawData=allScore,width = 4.8,height = 4,res = 500,useDingbats = F)
  
}















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



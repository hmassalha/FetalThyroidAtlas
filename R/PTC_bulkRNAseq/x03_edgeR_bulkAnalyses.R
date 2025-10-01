## bulkRNAseq DE Analyses with edgeR, to compare:
## 1. adult PTC vs adult normal thyroid (TCGA)
## 2. paediatric PTC vs paediatric normal thyroid (inhouse data)
## due to batch effects across dataset, cannot compare adult PTC vs paediatric PTC directly


outDir = '~/FetalThyroidAtlas/Results/2505/PTC_bulkRNAseq/DEG_tumour_vs_normal'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
setwd(outDir)



##-------------##
##  Libraries  ##
##-------------##

library(Seurat)
library(GenomicFeatures)
library(tidyverse)
library(edgeR)

source("~/FetalThyroidAtlas/R/utils/misc.R")
source("~/FetalThyroidAtlas/R/utils/sc_utils.R")
source("~/FetalThyroidAtlas/R/helperFunctions.R")



##-----------------------##
##        Params          #
##-----------------------##
tgtChrs = paste0('chr',c(1:22))

#Define genomic coordinates
geneMap_fp = '~/FetalThyroidAtlas/Results/geneMap.csv'
if(!file.exists(geneMap_fp)){
  gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
  txdb = makeTxDbFromGFF(gtf)
  gns = genes(txdb)
  
  geneMap = read.delim('/lustre/scratch126/cellgen/team292/hm11/with_Mi/T21_Oct24/HCA_GLNDrna14662854/output/GeneFull/filtered/features.tsv.gz',header = F,sep = '\t')
  colnames(geneMap) = c('ensID','geneSym','GEX')
  geneMap$geneSym[duplicated(geneMap$geneSym)] = paste0(geneMap$geneSym[duplicated(geneMap$geneSym)],'.1')
  geneMap$geneSym = gsub('_','-',geneMap$geneSym)
  geneMap$chr = as.character(seqnames(gns[match(geneMap$ensID,gns$gene_id)]))
  
  bioMart_annot = read.csv('~/lustre_mt22/generalResources/GRCh38_2020A_geneAnnotation.csv')
  geneMap = cbind(geneMap,bioMart_annot[match(geneMap$ensID,bioMart_annot$ensID),!colnames(bioMart_annot) %in% c('X',colnames(geneMap))])
  
}else{
  geneMap = read.csv(geneMap_fp,row.names = 1) %>% 
    dplyr::mutate(ensID = gene_id, geneSym = gene_name, gene_biotype = gene_type)
}



##----------------------##
##  Helper Functions  ####
##----------------------##


fit_model <- function(pb,colDat,formula,geneMap,groupID='group',MDS_groups = c('Genotype','donorID','sex'),pb_groupID='donorID',coef=NULL,mycontrast=NULL){
  
  if(!grepl('%s',formula,fixed=TRUE) || grepl('%s',sub('%s','',formula,fixed=TRUE),fixed=TRUE))
    stop("Invalid formula, must contain one instance of %s")
  #Convert it now in case it is still invalid in more subtle ways
  formula = as.formula(sprintf(formula,groupID))
  
  
  # create an edgeR object with counts and grouping factor
  y <- DGEList(pb, group = colDat[[groupID]])
  y$samples = cbind(y$samples,colDat[match(rownames(y$samples),colDat[[pb_groupID]]),!colnames(colDat) %in% colnames(y$samples)])
  
  ## Plot library sizes
  df = y$samples
  df[['pb_groupID']] = rownames(df)
  p = ggplot(df,aes(reorder(pb_groupID,`lib.size`),fill=group,y=`lib.size`))+
    geom_col()+
    scale_fill_manual(values = c(col25,pal34H))+
    theme_classic(base_size = 12)+xlab('')+
    theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
          axis.line = element_blank(),axis.text = element_text(color='black'),
          panel.border = element_rect(fill=F,colour = 'black'))
  
  
  
  print(p)
  
  
  
  # filter out genes with low counts
  print("Dimensions before subsetting:")
  print(dim(y))
  print("")
  all(y$samples$lib.size>1e5)
  
  
  
  keep <- filterByExpr(y)
  
  op = par(no.readonly = TRUE)
  par(mfrow = c(1, 2))
  hist(edgeR::cpm(y, log = TRUE), main = 'Unfiltered', xlab = 'logCPM')
  abline(v = log(1), lty = 2, col = 2)
  hist(edgeR::cpm(y[keep, ], log = TRUE), main = 'Filtered', xlab = 'logCPM')
  abline(v = log(1), lty = 2, col = 2)
  par(op)
  
  y <- y[keep, , keep.lib.sizes=FALSE]
  print("Dimensions after subsetting:")
  print(dim(y))
  print("")
  
  
  
  ## Normalization 
  y <- calcNormFactors(y)
  # MDS plot
  cluster<-as.factor(colDat[[groupID]][match(rownames(y$samples),colDat[[pb_groupID]])])
  df = plotMDS(y, pch=16,col=c(col25,pal34H)[cluster], main="MDS") 
  df = data.frame(x = df$x,y=df$y,pb_groupID = names(df$x))
  df = cbind(df[,c(1:3)],colDat[match(df$pb_groupID,colDat[[pb_groupID]]),!colnames(colDat) %in% colnames(df)])
  for(f in MDS_groups){
    df$group = df[[f]]
    p = ggplot(df,aes(x,y,col=group))  +
      geom_point()+
      theme_classic(base_size = 10)+
      scale_color_manual(values = c(col25,pal34H))+
      theme(axis.line = element_blank(),axis.text = element_text(color='black'),
            panel.border = element_rect(fill=F,colour = 'black'))
    print(p)
  }
  
  
  # create a design matrix: here we have multiple donors so also consider that in the design matrix
  design <- model.matrix(formula,data=y$samples)
  
  # estimate dispersion
  y <- estimateDisp(y, design = design,robust = T)
  # fit the model
  fit <- glmQLFit(y, design)
  plotQLDisp(fit)
  
  # Extract coefficients
  if(is.null(coef)){
    contrast = makeContrasts(mycontrast, levels=y$design)
    qlf<-glmQLFTest(fit, contrast=contrast)
  }else{
    qlf<-glmQLFTest(fit, coef = coef)
  }
  plotMD(qlf)
  abline(h=c(-1,1), col="blue")
  
  # get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
  tt <- topTags(qlf, n = Inf,p.value = 0.05,)
  tt <- tt$table
  tt = annotateGenes(tt,geneMap = geneMap)
  
  ## Calculate logCPM
  y$logCPM <- edgeR::cpm(y, log=TRUE, prior.count = 1)
  
  return(list("fit"=fit, "design"=design, "y"=y, 'tt'=tt))
}


plot_logCPM_byGroup = function(y,genes,geneMap,group='cancerType'){
  df = y$logCPM[geneMap$ensID[geneMap$geneSym %in% genes],]
  rownames(df) = geneMap$geneSym[match(rownames(df),geneMap$ensID)]
  df = as.data.frame(df)
  df$geneSym = rownames(df)
  df = pivot_longer(df,col=colnames(df)[colnames(df) != 'geneSym'],names_to='sampleID',values_to='logCPM')
  colDat = y$samples
  df = cbind(df,colDat[match(df$sampleID,colDat$sampleID),!colnames(colDat) %in% colnames(df)])
  df$group = paste0(df[[group]])
  p=ggplot(df,aes(group,logCPM,fill=group))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(size=0.7,width = 0.25)+
    scale_fill_manual(values = col25)+
    facet_wrap(vars(geneSym),scales = 'free_y')+
    theme_classic(base_size = 12)+xlab('')+
    theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
          axis.line = element_blank(),axis.text = element_text(color='black'),
          strip.background = element_blank(),
          panel.border = element_rect(fill=F,colour = 'black'))
  print(p)
  return(df)
}





##----------------------##
##  Import datasetes  ####
##----------------------##

library(TCGAbiolinks)

query <- TCGAbiolinks::GDCquery(project = "TCGA-THCA", 
                  data.category = "Simple Nucleotide Variation", 
                  data.type = "Masked Somatic Mutation")
TCGAbiolinks::GDCdownload(query)
maf <- TCGAbiolinks::GDCprepare(query)

##---- adult TCGA bulk RNAseq
tcga_sce_path = '~/FetalThyroidAtlas/Data/published_bulkRNAseq/TCGA_Thyroid/TCGA_Thyroid_bulkRNA_se.RDS'
tcga_sce = readRDS(tcga_sce_path)
tcga_rawCnt = assays(tcga_sce)[['unstranded']]
tcga_mdat = as.data.frame(colData(tcga_sce)) %>% dplyr::mutate(sampleID = barcode,cancerType=tissue_type)
#tcga_mdat = tcga_mdat[match(colnames(tcga_rawCnt),tcga_mdat$sampleID),]
#tcga_mdat$cancerType[tcga_mdat$cancerType == 'Classical/usual' & tcga_mdat$driver_classification == 'RET_fusion'] = 'RET_fusion_Classical'
tcga_mdat = tcga_mdat[tcga_mdat$tissue_type == 'Normal' | (!is.na(tcga_mdat$paper_fusionDriverGenes) & grepl('RET',tcga_mdat$paper_fusionDriverGenes)),]


# Determine sex
tpm_sex = assays(tcga_sce)[['tpm_unstranded']][geneMap$ensID[geneMap$geneSym %in% c('XIST','RPS4Y1')],tcga_mdat$sampleID]
rownames(tpm_sex) = geneMap$geneSym[match(rownames(tpm_sex),geneMap$ensID)]
sex = apply(tpm_sex,2,function(s){
  names(s)=rownames(tpm_sex)
  if(s['XIST'] > s['RPS4Y1']){
    return('F')
  }else{
    return('M')
  }
})
tcga_mdat$sex = sex[match(tcga_mdat$sampleID,names(sex))]

table(tcga_mdat$sex,tcga_mdat$gender)



##---- paediatric bulk RNAseq
sce_path = '~/FetalThyroidAtlas/Data/inhouse_bulk/inhouse_bulkRNA_fetalThyroid_paedPTC.RDS'
inhouse_sce = readRDS(sce_path)
inhouse_rawCnt = assays(inhouse_sce)[['counts_raw']]
inhouse_mdat = as.data.frame(colData(inhouse_sce))
inhouse_mdat = inhouse_mdat[match(colnames(inhouse_rawCnt),inhouse_mdat$sampleID),]
# Determine sex
tpm_sex = assays(inhouse_sce)[['counts_tpm']][geneMap$ensID[geneMap$geneSym %in% c('XIST','RPS4Y1')],]
rownames(tpm_sex) = geneMap$geneSym[match(rownames(tpm_sex),geneMap$ensID)]
sex = apply(tpm_sex,2,function(s){
  names(s)=rownames(tpm_sex)
  if(s['XIST'] > s['RPS4Y1']){
    return('F')
  }else{
    return('M')
  }
})
inhouse_mdat$sex = sex[match(inhouse_mdat$sampleID,names(sex))]

table(inhouse_mdat$sex)

##---- FThyrocytes - scRNAseq
fThy = readRDS('/lustre/scratch125/casm/team274sb/mt22/Thyroid/Data/fetalThyroid/fThyrocytes_2n_jul23.RDS')
fThy$annot = fThy$celltype
fThy$annot[fThy$annot == 'thy_Lumen-forming'] = 'TFC2'
fThy$annot[fThy$annot == 'thy_TH_processing'] = 'TFC1'


##---- REBC-THYR (Morton etal 2021) - validation cohort
rebc_fp = '~/FetalThyroidAtlas/Data/published_bulkRNAseq/REBC-THYR/REBC_THYR_2508.RData'
rebc_se = readRDS(rebc_fp)
rebc_rawCnt = assays(rebc_se)[['unstranded']]
rebc_mdat = as.data.frame(colData(rebc_se))
rebc_mdat = rebc_mdat[grepl('RET',rebc_mdat$WGS_CandidateDriverFusion),]


##--------------------------------------------------------##
##  ADULT - bulkRNAseq DEG between TCGA normal and PTC  ####
##--------------------------------------------------------##

##---- Import bulk RNAseq data
# subset to just normal and RET-fusion PTC samples
#tcga_mdat = tcga_mdat[tcga_mdat$cancerType %in% c('Normal','CCDC6_RET:Classical/usual','RET_fusion_Classical'),]
tcga_rawCnt = tcga_rawCnt[rownames(tcga_rawCnt) %in% geneMap$ensID[!is.na(geneMap$gene_biotype) & geneMap$gene_biotype %in% c('protein_coding')],
                          tcga_mdat$sampleID]

tcga_mdat$group = ifelse(tcga_mdat$cancerType == 'Normal','Normal','aPTC_RET')

##---- Perform edgeR
tcga_ptc_vs_normal = fit_model(pb=tcga_rawCnt,
                               colDat=tcga_mdat,
                               formula='~ %s + sex',
                               geneMap=geneMap,groupID='group',
                               MDS_groups = c('cancerType','group','sex'),
                               pb_groupID='sampleID',coef=2,mycontrast=NULL)

saveRDS(tcga_ptc_vs_normal,'tcga_retPTC_vs_normal_edgeR.RDS')
tcga_ptc_vs_normal = readRDS('tcga_retPTC_vs_normal_edgeR.RDS')

# Extract DEGs
deg = tcga_ptc_vs_normal[['tt']]
deg = deg[abs(deg$logFC) >= 1.5 & deg$FDR < 0.05 & 
            #deg$logCPM > -0.5 & 
            !grepl('^RPL|^RPS|^MT-|^MALAT1|^NEAT1|^GAPDH|^GABARAP|^MRPL|^MRPS|AC\\d+|LINC\\d+',deg$geneSym),]
deg$direction = ifelse(deg$logFC > 0,'aPTC_down','aPTC_up')
table(deg$direction)
deg = deg[order(abs(deg$logFC),decreasing = T),]

tcga_deg = deg
write.csv(tcga_deg,'tcga_retPTC_vs_normal_edgeR_DEG.csv')
tcga_deg = read.csv('tcga_retPTC_vs_normal_edgeR_DEG.csv')





##--------------------------------------------------------##
##  PAEDIATRIC - bulkRNAseq DEG between normal and PTC  ####
##--------------------------------------------------------##

##---- Import bulk RNAseq data
# subset to just normal foetal and paediatric samples
inhouse_mdat = inhouse_mdat[!grepl('foetal',inhouse_mdat$Tissue),]
# subset to just genes which are DEGs between ADULT PTC and normal tissues
inhouse_rawCnt = inhouse_rawCnt[rownames(inhouse_rawCnt) %in% geneMap$ensID[geneMap$geneSym %in% tcga_deg$geneSym],inhouse_mdat$sampleID]
inhouse_rawCnt = inhouse_rawCnt[,inhouse_mdat$sampleID]

##---- Perform edgeR
ptc_vs_normal = fit_model(pb=inhouse_rawCnt,
                          colDat=inhouse_mdat,
                          formula='~ %s + sex',
                          geneMap=geneMap,groupID='cancerType',
                          MDS_groups = c('sampleName','cancerType','sex'),
                          pb_groupID='sampleID',coef=2,mycontrast=NULL)
saveRDS(ptc_vs_normal,'Sanger_pPTC_vs_normal_TCGA.degs.only.RDS')
ptc_vs_normal = readRDS('Sanger_pPTC_vs_normal_TCGA.degs.only.RDS')

# Extract DEGs
deg = ptc_vs_normal[['tt']]
# deg = deg[abs(deg$logFC) >= 1 & deg$FDR < 0.05 & 
#             deg$logCPM > -0.5 & 
#             !grepl('^RPL|^RPS|^MT-|^MALAT1|^NEAT1|^GAPDH|^GABARAP|^MRPL|^MRPS|AC\\d+|LINC\\d+',deg$geneSym),]
deg$direction = ifelse(deg$logFC >0,'pPTC_up','pPTC_down')
table(deg$direction)
deg = deg[order(abs(deg$logFC),decreasing = T),]

# qlf<-glmQLFTest(ptc_vs_normal[['fit']], coef = 2)
# allGenes = topTags(qlf, n = Inf,p.value = 0.1)
# allGenes <- allGenes$table
# allGenes = annotateGenes(allGenes,geneMap = geneMap)


pPTC_deg = deg
write.csv(pPTC_deg,'Sanger_pPTC_vs_normal_TCGA.degs.only_DEGs.csv')
pPTC_deg = read.csv('Sanger_pPTC_vs_normal_TCGA.degs.only_DEGs.csv')




##----- Intersection between TCGA and paediatric Tum-vs-normal comparison ####
library(UpSetR)
upset(fromList(list('tcga_up' = tcga_deg$geneSym[tcga_deg$direction == 'aPTC_up'],
                    'tcga_down' = tcga_deg$geneSym[tcga_deg$direction == 'aPTC_down'],
                    'pPTC_up'=pPTC_deg$geneSym[pPTC_deg$direction == 'pPTC_up'],
                    'pPTC_down'=pPTC_deg$geneSym[pPTC_deg$direction == 'pPTC_down'])),text.scale = 2,nsets = 20)

## Unique pPTC genes
pPTC_unique = rbind(pPTC_deg[pPTC_deg$direction == 'pPTC_up' & !pPTC_deg$geneSym %in% tcga_deg$geneSym[tcga_deg$direction == 'aPTC_up'],],
                    pPTC_deg[pPTC_deg$direction == 'pPTC_down' & !pPTC_deg$geneSym %in% tcga_deg$geneSym[tcga_deg$direction == 'aPTC_down'],])

pPTC_opposite = rbind(pPTC_deg[pPTC_deg$direction == 'pPTC_up' & pPTC_deg$geneSym %in% tcga_deg$geneSym[tcga_deg$direction == 'aPTC_down'],],
                    pPTC_deg[pPTC_deg$direction == 'pPTC_down' & pPTC_deg$geneSym %in% tcga_deg$geneSym[tcga_deg$direction == 'aPTC_up'],])

## Unique TCGA genes
TCGA_unique = rbind(tcga_deg[tcga_deg$direction == 'aPTC_up' & !tcga_deg$geneSym %in% pPTC_deg$geneSym[pPTC_deg$direction == 'pPTC_up'],],
                    tcga_deg[tcga_deg$direction == 'aPTC_down' & !tcga_deg$geneSym %in% pPTC_deg$geneSym[pPTC_deg$direction == 'pPTC_down'],])

## Some sanity check plots
genes_toPlot = pPTC_unique$geneSym[pPTC_unique$direction == 'pPTC_down']
genes_toPlot = TCGA_unique$geneSym[TCGA_unique$direction == 'aPTC_down']
genes_toPlot = pPTC_opposite$geneSym
genes_toPlot = emt
length(genes_toPlot)
df = plot_logCPM_byGroup(ptc_vs_normal[['y']],genes = genes_toPlot[1:60],geneMap = geneMap,group='cancerType')
df = plot_logCPM_byGroup(tcga_ptc_vs_normal[['y']],genes = genes_toPlot[1:60],geneMap = geneMap,group='group')
df = plot_logCPM_byGroup(tcga_ptc_vs_normal[['y']],genes = c('TG','SLC26A7','IYD','ERG1','FOSB'),geneMap = geneMap,group='cancerType')
df = plot_logCPM_byGroup(ptc_vs_normal[['y']],genes = c('TG','FOXJ1'),geneMap = geneMap,group='cancerType')




fig4c_bulkDEG_TCGA_ptc.vs.normal = function(){
  pPTC_deg = read.csv('~/lustre_mt22/Thyroid/Results_v2/10_bulkEdgeR_thyrocytes_foetal_vs_paed/oct24/Sanger_pPTC_vs_normal_TCGA.degs.only_DEGs.csv')
  tcga_deg = read.csv('~/lustre_mt22/Thyroid/Results_v2/10_bulkEdgeR_thyrocytes_foetal_vs_paed/oct24/tcga_retPTC_vs_normal_edgeR_DEG.csv')
  
  ## Do weird column plot
  tcga_deg$deg_category = ifelse(tcga_deg$direction == 'aPTC_up' & tcga_deg$geneSym %in% pPTC_deg$geneSym[pPTC_deg$direction == 'pPTC_up'],'shared',
                                 ifelse(tcga_deg$direction == 'aPTC_down' & tcga_deg$geneSym %in% pPTC_deg$geneSym[pPTC_deg$direction == 'pPTC_down'],'shared','unique'))
  
  dd = tcga_deg %>% group_by(direction,deg_category) %>% summarise(nGene=n())
  dd$nGene[dd$direction == 'aPTC_down'] = -dd$nGene[dd$direction == 'aPTC_down']
  
  plotFun_nDEG_PTC.vs.Normal_inBulkSamples = function(noFrame=FALSE,noPlot=FALSE){
    p1 = ggplot(dd,aes(deg_category,nGene))+
      geom_col(aes(fill=direction),width = 0.6)+
      scale_fill_manual(values = c('#1a4a87','#a4282c','#5878a1','#b36d6e'))+
      geom_hline(yintercept = 0,col='black',lwd=0.5)+
      theme_classic(base_size = 11)+xlab('')+ylab('Number of differentially expressed genes')+
      theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),
            axis.line = element_blank(),#legend.position = 'bottom',
            axis.ticks = element_line(colour = 'black'),
            legend.title = element_text(size=10,colour = 'black'),
            legend.text = element_text(size=8,colour = 'black'),legend.key.size = unit(0.5,'cm'),
            axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,colour = 'black'),
            axis.text = element_text(color='black'))
    print(p1)  
  }
  saveFig(file.path(plotDir,'Fig4c_nDEG_Tum.vs.Normal_bulkSamples'),plotFun_nDEG_PTC.vs.Normal_inBulkSamples,rawData=dd,width = 4,height = 3.5,res = 500,useDingbats = F)
  
}






##--------------------------------------##
##    EnrichR on adult-unique DEGs    ####
##--------------------------------------##

## Perform enrichR on adultPTC-unique DEGs
upDEG_enriched <- enrichr(TCGA_unique$geneSym[TCGA_unique$direction == 'aPTC_up'], dbs)
downDEG_enriched <- enrichr(TCGA_unique$geneSym[TCGA_unique$direction == 'aPTC_down'], dbs)
## UP regulated module
allTerms_up = extract_enrichR.res(upDEG_enriched,module_direction = 'up',pVal_cutoff = 0.2,min_overlapFrac = 0.1,
                                  db_toExtract = names(upDEG_enriched)[c(5,6)])
allTerms_up$db = factor(allTerms_up$db,c("MSigDB_Hallmark_2020","KEGG_2021_Human"))
table(allTerms_up$db)
min(allTerms_up$nGene)
allTerms_down = extract_enrichR.res(downDEG_enriched,module_direction = 'down',pVal_cutoff = 0.2,min_overlapFrac = 0.1,
                                  db_toExtract = names(downDEG_enriched)[c(3,5,6)])
table(allTerms_down$db)
  



## Plot results - UP-regulated
plotFun_enrichR.up_pPTC_DEGs = function(noFrame=FALSE,noPlot=FALSE){
  p3=ggplot(allTerms_up,aes(Combined.Score,Term))+
    geom_point(aes(size = nGene,col=overlap2))+
    facet_grid(db~.,scales = 'free_y',space = 'free_y')+
    geom_segment(aes(x=0,xend=(Combined.Score),y=yStart,yend=yStart),col=grey(0))+
    xlab('Combined Score')+ylab('') + 
    ggtitle('Up-regulated DEGs in aPTC only')+
    scale_color_gradient(low='#eba9a9',high = '#a10505')+
    theme_bw(base_size = 11)+
    theme(panel.border = element_rect(fill=F,colour = 'black'),
          axis.text = element_text(colour = 'black',size=9.5),
          axis.ticks = element_line(colour = 'black'),
          panel.spacing.y = unit(0.3,'cm'),
          axis.line = element_blank())
  print(p3)
}

saveFig(file.path(plotDir,'Fig4_enrichR.up_aPTC.unique.DEGs'),plotFun_enrichR.up_pPTC_DEGs,rawData=allTerms_up,width = 8,height = 4,res = 500)

  
  
## Plot results - DOWN-regulated
allTerms_down = extract_enrichR.res(enrichR_result = downDEG_enriched,module_direction = 'down',pVal_cutoff = 0.05,min_overlapFrac = 0.05,
                                    db_toExtract = names(downDEG_enriched)[c(5,6)])
dim(allTerms_down)
## No significant term for KEGG and MSigDB

allTerms_down = extract_enrichR.res(downDEG_enriched,module_direction = 'down',pVal_cutoff = 0.1,min_overlapFrac = 0.05,
                                    db_toExtract = names(downDEG_enriched)[c(1,2,3)])

plotFun_enrichR.down_pPTC_DEGs = function(noFrame=FALSE,noPlot=FALSE){
  p3=ggplot(allTerms_down,aes(Combined.Score,Term))+
    geom_point(aes(size = nGene,col=overlap2))+
    facet_grid(db~.,scales = 'free_y',space = 'free_y')+
    geom_segment(aes(x=0,xend=(Combined.Score),y=yStart,yend=yStart),col=grey(0))+
    xlab('Combined Score')+ylab('') + 
    ggtitle('Down-regulated DEGs in pPTC')+
    scale_color_gradient(low='#eba9a9',high = '#a10505')+
    theme_bw(base_size = 11)+
    theme(panel.border = element_rect(fill=F,colour = 'black'),
          axis.text = element_text(colour = 'black',size=9.5),
          axis.ticks = element_line(colour = 'black'),
          axis.line = element_blank())
  print(p3)
}
  
saveFig(file.path(plotDir,'Fig4_enrichR.down_aPTC.unique.DEGs'),plotFun_enrichR.down_pPTC_DEGs,rawData=allTerms_down,width = 7,height = 3.7,res = 500)






##-----------------------------------------------------------------------------------##
##    Scoring for the enrichment of adultPTC-unique DEGs in normal fThy scRNAseq   ####
##-----------------------------------------------------------------------------------##

## Define the gene module
TCGA_unique = TCGA_unique[order(abs(TCGA_unique$logFC),decreasing = T),]
geneList = split(TCGA_unique$geneSym,TCGA_unique$direction)
length(geneList[[1]])
geneList[[1]] = geneList[[1]][1:300]
geneList[[2]] = geneList[[2]][1:300]
geneList[['aPTC_combined_1']] = c(paste0(geneList[[1]],'-'),paste0(geneList[[2]],'+'))
geneList[['aPTC_combined_2']] = c(paste0(geneList[[2]],'-'),paste0(geneList[[1]],'+'))

##----- Normal foetal
fThy <- UCell::AddModuleScore_UCell(fThy, features = geneList,ncores = 3)
fThy$dataset = 'fThy'
FeaturePlot(fThy,'aPTC_combined_2_UCell',cols = c(grey(0.2),grey(0.8),'red')) 

dd = fThy@meta.data[fThy$annot != 'Thyrocytes',c('cellID','annot','aPTC_up_UCell','aPTC_down_UCell','aPTC_combined_1_UCell','aPTC_combined_2_UCell','dataset')]
dd$group_fill = dd$annot
dd = dd[dd$annot != 'thy_Cycling',]
dd$group_fill = ifelse(dd$annot %in% c('fTFC1','aTFC1'),'TFC1',
                       ifelse(dd$annot %in% c('fTFC2','aTFC2'),'TFC2','others'))


write.csv(dd,'tcga_retPTC_vs_normal_edgeR_ModuleScore_inscRNAseq_2412.csv')

groupCol = c('PTC'='purple','TFC2'=grey(0.1),'TFC1'='orange','Thyrocytes'=grey(0.7),'others'='red')
plotFun_fTFC2_moduleScore_in.scFThySamples = function(noFrame=FALSE,noPlot=FALSE){

  if(noPlot & !noFrame){
    
    p = ggplot(dd,aes(annot,aPTC_down_UCell))+
      geom_boxplot(aes(fill=group_fill),outlier.size = 0.01,alpha=0.7,width=0.24)+
      scale_fill_manual(values = groupCol)+
      facet_grid(~dataset,scales = 'free_x')+
      theme_classic(base_size = 11)+
      theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 0.8),axis.line = element_blank(),
            axis.text = element_text(colour = 'black',size=9.5),
            axis.ticks = element_line(colour = 'black'),
            axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
            strip.background = element_blank()) +
      xlab('') + ylab('TCGA - aPTC down-regulated genes')
  }else{
    p = ggplot(dd,aes(annot,aPTC_down_UCell))+
      geom_quasirandom(width = 0.2,size=0.25,alpha=0.2,col=grey(0.7))+
      geom_boxplot(aes(fill=group_fill),outlier.size = 0.01,alpha=0.7,width=0.24)+
      scale_fill_manual(values = groupCol)+
      facet_grid(~dataset,scales = 'free_x')+
      theme_classic(base_size = 11)+
      theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 0.8),axis.line = element_blank(),
            axis.text = element_text(colour = 'black',size=9.5),
            axis.ticks = element_line(colour = 'black'),
            axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
            strip.background = element_blank()) +
      xlab('') + ylab('TCGA - aPTC down-regulated genes')
  }
  
  
  print(p)
}
saveFig(file.path(plotDir,'Fig4e_adultPTC.unique_moduleScore_in.scFThy.Samples'),plotFun_fTFC2_moduleScore_in.scFThySamples,rawData=dd,width = 3.7,height = 3.2,res = 500,useDingbats = F)




##--------------------------------##
##    Validation in REBC-THYR   ####
##--------------------------------##
source('~/FetalThyroidAtlas/R/helperFunctions.R')

gene_modules = read.csv('~/FetalThyroidAtlas/SupplementaryTables/TableS8_adult.PTCvsNormal_imprints_in_paed.PTCvsNormal.csv',row.names = 1)
gene_modules$ensID = geneMap$ensID[match(gene_modules$gene_symbol..adult.PTC.vs.normal.,geneMap$gene_name)]
moduleList = list(aPTC_up = gene_modules$ensID[gene_modules$direction..adult.PTC.vs.normal. == 'upregulated_in_adultPTC'],
                   aPTC_down = gene_modules$ensID[gene_modules$direction..adult.PTC.vs.normal. == 'downregulated_in_adultPTC'],
                   pPTC_up = gene_modules$ensembl_ID..paed.PTC.vs.normal.[!is.na(gene_modules$direction..paed.PTC.vs.normal.) & gene_modules$direction..paed.PTC.vs.normal. == 'upregulated_in_paed.PTC'],
                   pPTC_down = gene_modules$ensembl_ID..paed.PTC.vs.normal.[!is.na(gene_modules$direction..paed.PTC.vs.normal.) & gene_modules$direction..paed.PTC.vs.normal. == 'downregulated_in_paed.PTC'])
moduleList[['adult_up_only']] = moduleList[['aPTC_up']][!moduleList[['aPTC_up']] %in% moduleList[['pPTC_up']]]
moduleList[['adult_down_only']] = moduleList[['aPTC_down']][!moduleList[['aPTC_down']] %in% moduleList[['pPTC_down']]]

##--- import bulk counts and calculate cpmCnt in xx01_moduleScoring.R
bulkRNA = import_bulkRNA_thyroid(bulk_sources = c('inhouse'='~/FetalThyroidAtlas/Data/inhouse_bulk/inhouse_bulkRNA_fetalThyroid_paedPTC.RDS',
                                                  'TCGA_Thyroid'='~/FetalThyroidAtlas/Data/published_bulkRNAseq/TCGA_Thyroid/TCGA_Thyroid_bulkRNA_se.RDS',
                                                  'REBC_THYR' = '~/FetalThyroidAtlas/Data/published_bulkRNAseq/REBC-THYR/REBC_THYR_2508.RData',
                                                  'He2021'='~/FetalThyroidAtlas/Data/published_bulkRNAseq/He_etal_21/aPTC_He_2021_se.RDS',
                                                  'Lee2024' = '~/FetalThyroidAtlas/Data/published_bulkRNAseq/Lee_etal_24/aPTC_Lee_2024_se.RDS'),
                                 gene_map = gene_map)
bulk_samples = bulkRNA[['bulk_samples']]
rownames(bulk_samples) = bulk_samples$sampleID
cpmCnt = bulkRNA[['cpmCnt']]
tpmCnt = bulkRNA[['tpm_count']]
rawCnt = bulkRNA[['raw_count']]

##---  Score the modules -----##
rebc_se = readRDS('~/FetalThyroidAtlas/Data/published_bulkRNAseq/REBC-THYR/REBC_THYR_2508.RData')
rebc_mdat = colData(rebc_se)
tpmCnt = assays(rebc_se)[['tpm_unstranded']]
colnames(tpmCnt) = colData(rebc_se)[['File.ID']]
table(rownames(tpmCnt) %in% gene_modules$ensID[gene_modules$direction..adult.PTC.vs.normal. != 'not_DE'])
table(gene_modules$ensID[gene_modules$direction..adult.PTC.vs.normal. != 'not_DE'] %in% rownames(tpmCnt))

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

allScore$cancerType[allScore$source2 %in% c('REBC_THYR_paed','REBC_THYR_adult') & 
                      allScore$cancerType != 'Normal' & 
                      allScore$sampleID %in% rebc_mdat$File.ID[grepl('RET',rebc_mdat$WGS_CandidateDriverFusion)]] = paste0(allScore$cancerType[allScore$source2 %in% c('REBC_THYR_paed','REBC_THYR_adult') & 
                                                                                                                                                 allScore$cancerType != 'Normal' & 
                                                                                                                                                 allScore$sampleID %in% rebc_mdat$File.ID[grepl('RET',rebc_mdat$WGS_CandidateDriverFusion)]],'_RET')
allScore$others = rebc_mdat$DriverGroup_CoMutFig[match(allScore$sampleID,rebc_mdat$File.ID)]

library(ggbeeswarm)
ggplot(allScore,aes(cancerType,TotalScore))+
  geom_boxplot(outlier.shape = NA,aes(fill=cancerType))+
  geom_quasirandom(aes(col=others),width = 0.3,size=0.1,alpha=1)+
  facet_grid(moduleType ~ source2,scales = 'free',space = 'free')+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))

rebc_mdat$WGS_CandidateDriverFusion

allScore_sub = allScore[allScore$cancerType == 'Normal',]
allScore_sub = merge(allScore_sub,rebc_mdat,by.x='Row.names',by.y='File.ID')
allScore_sub = allScore_sub[allScore_sub$moduleType == 'aPTC_up',]
libsize = colSums(assays(rebc_se)[['unstranded']])
names(libsize) = colData(rebc_se)[['File.ID']]
allScore_sub$libSize = libsize[match(allScore_sub$sampleID,names(libsize))]
pdf('~/rebc_test.pdf')
for (i in 2:ncol(allScore_sub)) {
  p <- ggplot(allScore_sub, aes(x = TotalScore, y = libSize)) +
    geom_point(aes(col = allScore_sub[[i]])) +
    scale_y_log10() +
    theme_bw() +
    theme(legend.position = "none") +
    labs(col = colnames(allScore_sub)[i])  # label legend with column name
  
  print(p)   # needed inside a loop to display plots
}
dev.off()


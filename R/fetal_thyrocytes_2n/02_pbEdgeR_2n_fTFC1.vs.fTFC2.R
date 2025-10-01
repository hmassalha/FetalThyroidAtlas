## Perform pseudobulk DE analysis on scRNAseq diploid thyrocytes (from 27 foetuses), using edgeR
## to compare fTFC1 vs fTFC2 ##

outDir = '~/lustre_mt22/Thyroid/Results_v2/2.1_DEG_2n_fTFC1.vs.fTFC2/oct24'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
setwd(outDir)



##-------------##
##  Libraries  ##
##-------------##
library(tidyverse)
library(Seurat)
library(GenomicFeatures)
library(RColorBrewer)
library(edgeR)

source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
#source("~/lustre_mt22/generalScripts/utils/pseudobulk.R")



##-----------------------##
##        Params          #
##-----------------------##
plotDir = outDir
keepCylcingCells=F
tgtChrs = paste0('chr',c(1:22))

#Define genomic coordinates
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
geneMap = read.csv('~/FetalThyroidAtlas/Results/geneMap.csv')
geneMap$gene_biotype = geneMap$gene_type
geneMap$geneSym = geneMap$gene_name
geneMap$ensID = geneMap$gene_id
geneMap$chr = geneMap$seqnames

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
  df = cbind(df,colDat[match(df$pb_groupID,colDat[[pb_groupID]]),!colnames(colDat) %in% colnames(df)])
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






##----------------------------------------------------##
##  Import 2n-T21 age-matched fThy seurat object    ####
##----------------------------------------------------##
tissue = 'thyroid'

## Import foetal 2n thyrocytes only
srat_in_fp = '/lustre/scratch125/casm/team274sb/mt22/Thyroid/Data/fetalThyroid/fThyrocytes_2n_jul23.RDS'
srat_in_fp = '~/FetalThyroidAtlas/Data/fThyrocytes_2n_atlas.RDS'

if(!file.exists(srat_in_fp)){
  stop(sprintf('Cannot find input directory for annotated AK-REF merged sratObj for tissue %s \nThe path is: \n%s',tissue, srat_in_fp))
  # srat = Read10X('~/lustre_mt22/Thyroid/Data/agematched_2nT21/oct24')
  # srat = CreateSeuratObject(srat)
  # srat = standard_clustering(srat)
  # srat$cellID = rownames(srat@meta.data)
  # umap_coord = read.csv('~/lustre_mt22/Thyroid/Data/agematched_2nT21/oct24/UMAP.csv.gz',row.names = 1)
  # colnames(umap_coord) =c('UMAP_1','UMAP_2')
  # rownames(umap_coord) = srat$cellID
  # srat@reductions$umap@cell.embeddings = as.matrix(umap_coord)
  # mdat = read.csv('~/lustre_mt22/Thyroid/Data/agematched_2nT21/oct24/metadata.csv.gz',row.names = 1)
  # srat@meta.data = cbind(srat@meta.data,mdat[match(srat$cellID,rownames(mdat)),!colnames(mdat) %in% colnames(srat@meta.data)])
  # DimPlot(srat,group.by = 'donorID')
  # saveRDS(srat,srat_in_fp)
}else{
  srat = readRDS(srat_in_fp)   
  srat$tissue = tissue
  srat$cellID = rownames(srat@meta.data)
}

## check that only cells in G1 are retained
if(!keepCylcingCells){
  cyclingCells = srat$cellID[grepl('Cycling',srat$celltype)]
  if(length(cyclingCells) > 0){
    message(sprintf('Removing cycling cells from srat for tissue %s',tissue))
    srat = subset(srat, subset = cellID %in% srat$cellID[!srat$cellID %in% cyclingCells])
  }
}

## Adjust some metadata
srat$finalAnn = srat$celltype
srat$finalAnn[srat$finalAnn == 'thy_TH_processing'] = 'fTFC1'
srat$finalAnn[srat$finalAnn == 'thy_Lumen-forming'] = 'fTFC2'
srat$Genotype = srat$karyotype
srat$donorID = srat$donor
srat$age_group = ifelse(srat$pcw %in% c(9,10),'9-10',
                        ifelse(srat$pcw %in% c(11,12,13),'11-13','14-20'))
srat$group = paste0(srat$donor,'_',srat$finalAnn)





##-------------------------------------------##
##  ~ Celltype + donorID (all age_group)   ####
##-------------------------------------------##


outDir = '~/lustre_mt22/Thyroid/Results_v2/2.1_DEG_2n_fTFC1.vs.fTFC2/oct24/celltype_donorID_allAgeGroup'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
setwd(outDir)


result_fp = file.path(outDir,'pb_edgeR_allAgeGroup.RDS')
if(file.exists(result_fp)){
  out = readRDS(result_fp)
}else{
  out = list()
  
  message(sprintf("\n\n------- Consider thyrocytes from all ageGroup"))
  srat.sub = srat
  
  #Check we have a reasonable number of counts from both settings. ie Get number of cells for the relevant cell type and genenotype
  nCellsGroup = table(factor(srat.sub@meta.data$finalAnn))
  
  if((!all(nCellsGroup>=50))){
    message(sprintf('Low number of cells detected'))
    print(nCellsGroup)
    next
  }
  
  #Check how many from individual donors
  nCells = table(paste0(srat.sub@meta.data$group))
  if(sum(nCells>50)<3){
    message(sprintf("Too few effective replicates.  Skipping..."))
    print(nCells)
    next
  }
  message("Found the following number of high confidence cells")
  print(nCells)
  
  
  # remove individuals with < 30 cells
  group_toKeep = names(nCells[nCells >= 30])
  
  # if(ageGroup == '9-10'){
  #   group_toKeep = group_toKeep[!group_toKeep %in% c('Srv39_fTFC1','Srv39_fTFC2')]
  # }
  #OK, we're going ahead, create the needed objects
  toc = srat.sub@assays$RNA@counts[,colnames(srat.sub@assays$RNA@counts) %in% srat.sub$cellID[srat.sub$group %in% group_toKeep]]
  toc = toc[rownames(toc) %in% geneMap$geneSym,]
  m = match(rownames(toc),geneMap$geneSym)
  sum(is.na(m))
  rownames(toc) = geneMap$ensID[m]
  mDat = srat.sub@meta.data[match(colnames(toc),srat.sub$cellID),c('cellID','donorID','pcw','age_group','sex','group','finalAnn')]
  mDat = mDat[match(colnames(toc),mDat$cellID),]
  rownames(mDat) = colnames(toc)
  
  # check that rownames(mDat) is in correct order
  if(!all(rownames(mDat) == mDat$cellID)){
    stop(sprintf('Incorrect mDat cellID order for tissue %s cell type %s genotype %s. Please check!',tissue, tgtCell, geno))
  }
  
  # Only keep genes present in gns
  toc = toc[rownames(toc) %in% geneMap$ensID[geneMap$chr %in% tgtChrs & 
                                               !is.na(geneMap$gene_biotype) & 
                                               geneMap$gene_biotype == 'protein_coding'],]
  #toc = toc[rownames(toc) %in% gns$gene_id,]
  coords = gns[rownames(toc)]
  
  
  donorID = 'group'
  group = 'finalAnn'
  
  ## Drop irrelevant genes
  #Uninteresting genes
  w = which(as.character(seqnames(coords)) %in% tgtChrs)
  coords = coords[w]
  seqlevels(coords) = tgtChrs
  toc = toc[coords$gene_id,]
  #Non-expressed genes
  w = which(rowSums(toc>0)>=3)
  coords = coords[w]
  toc = toc[w,]
  
  ##=========================##
  # Get genomic coordinates
  #Order aa and bb by genomic coords
  #If it's positive stranded or unknown use start, else end
  coords$TSS = ifelse(strand(coords)=='-',end(coords),start(coords))
  o = order(seqnames(coords),coords$TSS)
  coords = coords[o]
  toc = toc[o,]
  
  ##=========================##
  # Create pseudobulk 
  pb = do.call(cbind,lapply(split(colnames(toc),mDat[,donorID]),function(e) rowSums(toc[,e,drop=FALSE])))
  
  
  colDat = mDat[match(colnames(pb),mDat[,donorID]),]
  rownames(colDat) = colDat[,donorID]
  
  
  ##=========================##
  # Fit EDGER model
  #mycontrast = 'finalAnnfTFC2'
  pdf(file.path(outDir,paste0('pbEdgeR_plots.pdf')))
  out[[g]] = fit_model(pb,colDat,formula = '~ %s + donorID',geneMap=geneMap,
                       MDS_groups = c('finalAnn','donorID','sex','age_group'),
                       groupID='finalAnn',pb_groupID = donorID,coef = 2)
  dev.off()
  
  saveRDS(out,result_fp)
  
}



## Extract log2FC and p-val for all genes
coef=2
allGenes = data.frame()
for(i in 1:length(out)){
  y = out[[i]][['y']]
  pb = y$counts
  write.csv(pb,'pb_2n_fTFC1.vs.fTFC2.csv',row.names = T,col.names = T)
  
  fit = out[[i]][['fit']]
  # Extract coefficients
  qlf<-glmQLFTest(fit, coef = coef)
  # get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
  tt <- topTags(qlf, n = Inf)
  tt <- tt$table
  tt = annotateGenes(tt,geneMap = geneMap)
  tt$comp = names(out)[i]
  allGenes = rbind(allGenes,tt)
}

allGenes$isDEG = (abs(allGenes$logFC) >= 0.2 & allGenes$FDR < 0.05)
allGenes$group = ifelse(allGenes$isDEG,'DEG','nonDEG')
allGenes$direction = ifelse(allGenes$logFC > 0,'fTFC2_up','fTFC2_down')
allGenes$comp = '~celltype+donorID (allAgeGroup)'

## Add percentage of cells from each cell type expressing each gene
allGenes$pct_fTFC1 = 100*rowSums(srat@assays$RNA@counts[allGenes$geneSym,srat$cellID[srat$finalAnn == 'fTFC1']] > 0) / length(srat$cellID[srat$finalAnn == 'fTFC1'])
allGenes$pct_fTFC2 = 100*rowSums(srat@assays$RNA@counts[allGenes$geneSym,srat$cellID[srat$finalAnn == 'fTFC2']] > 0) / length(srat$cellID[srat$finalAnn == 'fTFC2'])
allGenes$pct_diff = allGenes$pct_fTFC2 - allGenes$pct_fTFC1

write.csv(allGenes,file.path(outDir,'allGenes_log2FC_allAgeGroup.csv'))
allGenes = read.csv(file.path(outDir,'allGenes_log2FC_allAgeGroup.csv'))

deg = rbind(allGenes[allGenes$FDR < 0.01 & abs(allGenes$logFC) > 1 & allGenes$pct_fTFC2 > 20 & allGenes$direction == 'fTFC2_up',],
            allGenes[allGenes$FDR < 0.01 & abs(allGenes$logFC) > 1 & allGenes$pct_fTFC1 > 20 & allGenes$direction == 'fTFC2_down',])
table(deg$direction)
min(abs(deg$pct_diff))
table(deg$direction,deg$isTF)



##--------------------------------------------------##
##    Plot expression of DEGs in scRNAseq data    ####
##--------------------------------------------------##
genes_toPlot = c(deg$geneSym[abs(deg$pct_diff) > 20 & deg$direction == 'fTFC2_down'][1:100],
                 deg$geneSym[abs(deg$pct_diff) > 20 & deg$direction == 'fTFC2_up'][1:100])
length(genes_toPlot)

Idents(srat) = srat$finalAnn
srat$group = paste0(srat$finalAnn,':',srat$age_group,':',srat$donor)
DotPlot(srat,idents = c('fTFC1','fTFC2'),group.by = 'finalAnn',
        features = genes_toPlot
)+RotatedAxis() + 
  theme(axis.text.x = element_text(size=8,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size=8),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = 'right') + xlab('') + ylab('')






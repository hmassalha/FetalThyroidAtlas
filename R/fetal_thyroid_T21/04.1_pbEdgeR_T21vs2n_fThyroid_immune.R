## Perform pseudobulk DE analysis using edgeR 
## Comparing T21 - vs - 2n fetal immune cell types within the thyroid

setwd('~/FetalThyroidAtlas')

outDir = 'Results/2505/DEG_T21vs2n_immune'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}



##-------------##
##  Libraries  ##
##-------------##
library(Seurat)
library(tidyverse)
library(Matrix)
library(GenomeInfoDb)
library(GenomicFeatures)
library(GenomicRanges)
library(edgeR)

source("R/utils/misc.R")
#source("R/utils/sc_utils.R")
source("R/utils/04_pbEdgeR_helperFunctions.R")

# Create output directory
outDir_fp = file.path(outDir)  
plotDir = outDir_fp
if(!dir.exists(plotDir)){
  message('Making plot directory')
  dir.create(plotDir,recursive = T)
}



##-----------------------##
##        Params          #
##-----------------------##
#Define genomic coordinates
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
txdb = GenomicFeatures::makeTxDbFromGFF(gtf)
gns = GenomicFeatures::genes(txdb)

geneMap = read.csv('~/FetalThyroidAtlas/Results/geneMap.csv',row.names = 1)
geneMap$geneSym = geneMap$gene_name
geneMap$ensID = geneMap$gene_id
geneMap$chr = geneMap$seqnames
geneMap$gene_biotype = geneMap$gene_type
geneMap = geneMap[geneMap$gene_id %in% names(gns),]

keepCylcingCells=F
tissue='thyroid_imm'
tgtChrs = paste0('chr',c(1:22))

##--------------------------------------##
##  Import thyrocyte seurat object    ####
##--------------------------------------##
# Import thyroid scRNA-seq objects
fThyroid_2nT21_fp = 'Data/fThyroid_2nT21_agematched_atlas.RDS'
srat = readRDS(fThyroid_2nT21_fp)   
srat$cellID = rownames(srat@meta.data)
srat$Genotype = srat$karyotype
srat$donorID = srat$donor
srat$age_group = gsub('-','_',srat$age_group)
srat$finalAnn = srat$celltype

# Subset to just immune cell types
srat = subset(srat,subset = cellID %in% srat$cellID[grepl('^imm_',srat$celltype)])
### check that only cells in G1 are retained
if(!keepCylcingCells){
  cyclingCells = srat$cellID[grepl('Cycling',srat$celltype)]
  if(length(cyclingCells) > 0){
    message(sprintf('Removing cycling cells from srat for tissue %s',tissue))
    srat = subset(srat, subset = cellID %in% srat$cellID[!srat$cellID %in% cyclingCells])
  }
}



##-----------------------------##
##    Karyotype + ageGroup   ####
##-----------------------------##

out = list()
pb_byCT = list()
minCell = 20

for(tgtCell in unique(srat$finalAnn)){
  
  message(sprintf("\n\n------- Consider cell type %s",tgtCell))
  srat.sub = subset(srat,subset = cellID %in% srat$cellID[srat$finalAnn == tgtCell])
  
  #Check we have a reasonable number of counts from both settings. ie Get number of cells for the relevant cell type and genenotype
  nCellsGroup = table(factor(srat.sub@meta.data$Genotype))
  
  if((!all(nCellsGroup>=minCell))){
    message(sprintf('Low number of cells detected'))
    print(nCellsGroup)
    next
  }
  
  #Check how many from individual donors
  nCells = table(srat.sub@meta.data$donorID)
  if(sum(nCells>=minCell)<3){
    message(sprintf("Too few effective replicates.  Skipping..."))
    print(nCells)
    next
  }
  message("Found the following number of high confidence cells")
  print(nCells)
  
  
  # remove individuals with < minCell cells
  donorID_toKeep = names(nCells[nCells >= minCell])
  
  # Make sure that we still have at least 2 groups for comparison afer this step
  if(n_distinct(srat$Genotype[srat$donor %in% donorID_toKeep]) < 2){
    message(sprintf('Not enough groups for comparison'))
    print(unique(srat$Genotype[srat$donor %in% donorID_toKeep]))
    next
  }
  
  
  #OK, we're going ahead, create the needed objects
  toc = srat.sub@assays$RNA@counts[,colnames(srat.sub@assays$RNA@counts) %in% srat.sub$cellID[srat.sub$donorID %in% donorID_toKeep]]
  toc = toc[rownames(toc) %in% geneMap$geneSym,]
  m = match(rownames(toc),geneMap$geneSym)
  sum(is.na(m))
  rownames(toc) = geneMap$ensID[m]
  mDat = srat.sub@meta.data[match(colnames(toc),srat.sub$cellID),c('cellID','donorID','Genotype','pcw','age_group','sex')]
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
  coords = gns[rownames(toc)]
  
  formula = '~ %s'
  donorID = 'donorID'
  group = 'Genotype'
  
  ## Drop irrelevant genes
  #Uninteresting genes
  w = which(as.character(GenomicRanges::seqnames(coords)) %in% tgtChrs)
  coords = coords[w]
  GenomeInfoDb::seqlevels(coords) = tgtChrs
  toc = toc[coords$gene_id,]
  
  #Non-expressed genes: these are genes expressed in <10 cells or <10% cells in all groups
  min_cells = 10
  w = which(Matrix::rowSums(toc>0)>=min_cells)
  print(sprintf('Removing %d genes',nrow(toc) - length(w)))
  
  coords = coords[w]
  toc = toc[w,]
  
  
  # Remove genes expressed in < 10% of cells in at least 1 groups
  mDat$group = paste0(mDat[,'Genotype'],':',mDat$age_group)
  percCell_expressed_perGene = do.call(cbind,lapply(split(colnames(toc),mDat[,'group']),function(e){
    tmp = data.frame(percCell_expressed = 100*Matrix::rowSums(toc[,e,drop=FALSE] > 0) / length(e))
    colnames(tmp) = paste0(colnames(tmp),'_',unique(mDat$group[mDat$cellID %in% e]))
    return(tmp)
  } ))
  percCell_expressed_perGene$max_perc = apply(percCell_expressed_perGene,1,max)
  
  # only keep genes expressed >=10% cells in at least 1 group
  min_percCell = 10
  genes_toKeep = rownames(percCell_expressed_perGene[percCell_expressed_perGene$max_perc >=min_percCell,])
  print(sprintf('Removing %d genes',nrow(toc) - length(genes_toKeep)))
  coords = coords[genes_toKeep]
  toc = toc[genes_toKeep,]
  
  
  ##=========================##
  # Get genomic coordinates
  #Order aa and bb by genomic coords
  #If it's positive stranded or unknown use start, else end
  coords$TSS = ifelse(GenomicRanges::strand(coords)=='-',GenomicRanges::end(coords),GenomicRanges::start(coords))
  o = order(GenomicRanges::seqnames(coords),coords$TSS)
  coords = coords[o]
  toc = toc[o,]
  
  ##=========================##
  # Create pseudobulk 
  ## Split each donor into 3 pb
  ## mDat = mDat %>% group_by(donorID) %>% mutate(pb_ID = (seq(1:n())%%2))
  ## mDat$donorID = paste0(mDat$donorID,'_',mDat$pb_ID)
  pb = do.call(cbind,lapply(split(colnames(toc),mDat[,donorID]),function(e) Matrix::rowSums(toc[,e,drop=FALSE])))
  pb_byCT[[tgtCell]] = pb
  colDat = mDat[match(colnames(pb),mDat[[donorID]]),]
  colDat$Genotype[colDat$Genotype == '2n'] = 'diploid'
  colDat$Genotype = factor(colDat$Genotype,c('diploid','T21'))
  rownames(colDat) = colDat[[donorID]]
  colDat$group = colDat[[group]]
  
  
  
  
  
  ##=========================##
  # Fit EDGER model
  pdf(file.path(outDir,paste0(tgtCell,'_pbEdgeR_plots.pdf')))
  out[[tgtCell]] = fit_model(pb,colDat,formula = '~ 0 + %s + age_group',geneMap,groupID='group',MDS_groups = c('Genotype','donorID','sex','age_group'))
  dev.off()
}



saveRDS(out,file.path(outDir,'pb_edgeR_byCT.RDS'))


out = readRDS(file.path(outDir,'pb_edgeR_byCT.RDS'))
## Extract log2FC and p-val for all genes
allGenes = data.frame()
for(i in 1:length(out)){
  y = out[[i]][['y']]
  fit = out[[i]][['fit']]
  # Extract coefficients
  contrast = makeContrasts(groupT21-groupdiploid, levels=y$design)
  qlf<-glmQLFTest(fit, contrast=contrast)
  # get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
  tt <- topTags(qlf, n = Inf)
  tt <- tt$table
  tt = annotateGenes(tt,geneMap = geneMap)
  tt$comp = names(out)[i]
  allGenes = rbind(allGenes,tt)
}

allGenes$isDEG = (abs(allGenes$logFC) >= 0.2 & allGenes$FDR < 0.05)
allGenes$ischr21 = (allGenes$chr == 'chr21')
allGenes$group = ifelse(allGenes$ischr21,'chr21',ifelse(allGenes$isDEG,'DEG','nonDEG'))

colnames(allGenes)[colnames(allGenes) == 'comp'] = 'comparison'
table(allGenes$comparison)
write.csv(allGenes,file.path(outDir,'allGenes_log2FC_splitbyCT.csv'))
#write.csv(allGenes,file.path(outDir,'allGenes_log2FC_splitbyCT.csv'))

allGenes = read.csv(file.path(outDir,'allGenes_log2FC_splitbyCT.csv'))
allGenes$comp = allGenes$comparison
## Do volcano plots
ggplot(allGenes,aes(logFC,-log10(FDR),col=group))+
  geom_point(size=0.1)+
  scale_color_manual(values = c('purple','red',grey(0.7)))+
  geom_hline(yintercept = -log10(0.05),lty=2,linewidth=0.1)+
  geom_vline(xintercept = c(-0.2,0.2),lty=2,linewidth=0.1)+
  facet_wrap(vars(comp),nrow=1)+
  theme_classic(base_size = 13)+xlab('log2FC')+ylab('- log10 (FDR)')+
  theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),
        axis.line = element_blank(),#legend.position = 'bottom',
        axis.ticks = element_line(colour = 'black'),
        legend.title = element_text(size=10,colour = 'black'),
        legend.text = element_text(size=8,colour = 'black'),legend.key.size = unit(0.5,'cm'),
        axis.text = element_text(color='black'),strip.background = element_blank(),strip.text = element_text(size=13,colour = 'black'))

## Plot log2FC by chromosome
allGenes$chr = factor(allGenes$chr,paste0('chr',c(1:22)))
ggplot(allGenes,aes(chr,logFC))+
  geom_jitter(data=allGenes[allGenes$isDEG == F,],size=0.1,width=0.3,col=grey(0.7))+
  geom_jitter(data=allGenes[allGenes$isDEG == T,],size=0.2,width=0.3,col='red')+
  geom_boxplot(alpha=0.3,outlier.shape = NA,width=0.75)+
  #scale_color_manual(values = c(grey(0.7),'red'))+
  #geom_hline(yintercept = -log10(0.05),lty=2,linewidth=0.1)+
  geom_hline(yintercept = c(-0.2,0.2),lty=2,linewidth=0.1)+
  facet_wrap(vars(comp),nrow=2)+
  theme_classic(base_size = 13)+ylab('log2FC')+xlab('')+
  theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),
        axis.line = element_blank(),#legend.position = 'bottom',
        axis.ticks = element_line(colour = 'black'),
        legend.title = element_text(size=10,colour = 'black'),
        legend.text = element_text(size=8,colour = 'black'),legend.key.size = unit(0.5,'cm'),
        axis.text = element_text(color='black'),strip.background = element_blank(),strip.text = element_text(size=13,colour = 'black'))



# Plot number of DEGs

deg = data.frame()
for(i in 1:length(out)){
  tt = out[[i]][['tt']]
  tt$comp = names(out)[i]
  deg = rbind(deg,tt)
}

deg = deg[abs(deg$logFC)>= 0.2 & deg$FDR < 0.05,]
deg$direction = ifelse(deg$logFC > 0,'T21_up','T21_down')
df = deg %>% group_by(comp,direction) %>% summarise(nGene = n())
df$nGene[df$direction == 'T21_down'] = -df$nGene[df$direction == 'T21_down']
ggplot(df,aes('',nGene,fill=direction))+
  geom_col(width = 0.58)+
  geom_text(aes(label = nGene), 
            position = position_stack(vjust = 0.5), 
            vjust = ifelse(df$nGene > 0, -0.5, 1.5), 
            size = 3) +
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = c('#1a4a87','#a4282c','#5878a1','#b36d6e'))+
  facet_wrap(vars(comp),nrow=1)+
  theme_classic(base_size = 13)+xlab('')+ylab('# DEGs')+
  theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),
        axis.line = element_blank(),#legend.position = 'bottom',
        axis.ticks = element_line(colour = 'black'),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size=10,colour = 'black'),
        legend.text = element_text(size=8,colour = 'black'),legend.key.size = unit(0.5,'cm'),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,colour = 'black'),
        axis.text = element_text(color='black'),strip.background = element_blank(),strip.text = element_text(size=10))


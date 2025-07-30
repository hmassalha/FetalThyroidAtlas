## Perform pseudobulk on scRNAseq data using edgeR on T21 and diploid fetal thyrocytes in single-cell RNA-seq dataset ###
## oct24 is the current version used in the Manuscript
## nov24 we tried to add additional stringent gene filters, to reduce number of genes in stat test --> see if FDR values for thyroid hormone genes are any lower so that we might be able to comment on them in the MS.
##       unfortunately no, FDR values range from 0.2xx to 0.99 for these genes across both analyses
## --> we should stick with oct24 version for now

setwd('~/FetalThyroidAtlas')

outDir = 'Results/2505/DEG_T21vs2n_thyrocytes'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}



##-------------##
##  Libraries  ##
##-------------##
library(Seurat)
library(GenomicFeatures)

library(tidyverse)
library(readxl)
library(reshape2)
library(zoo)
library(RColorBrewer)
library(ComplexHeatmap)


library(edgeR)
library(MAST)

source("R/utils/misc.R")
source("R/utils/sc_utils.R")
source("R/utils/04_pbEdgeR_helperFunctions.R")
#source("R/utils/pseudobulk.R")


##-----------------------##
##        Params          #
##-----------------------##
#Define genomic coordinates
geneMap = read.csv('~/FetalThyroidAtlas/Results/geneMap.csv')

kMerge=11
keepCylcingCells=F
ageMatch = T
tgtChrs = paste0('chr',c(1:22))


##--------------------------------------##
##  Import thyrocyte seurat object    ####
##--------------------------------------##

# Create output directory
outDir_fp = file.path(outDir)  
plotDir = outDir_fp
if(!dir.exists(plotDir)){
  message('Making plot directory')
  dir.create(plotDir,recursive = T)
}

# Import thyroid scRNA-seq objects
fThyroid_2nT21_fp = 'Data/fThyroid_2nT21_agematched_atlas.RDS'
srat = readRDS(fThyroid_2nT21_fp)   


### Plot proportion of fTFC1 vs fTFC2 in 2n vs T21
dd = srat@meta.data %>% filter(celltype %in% c('thy_TH_processing','thy_Lumen-forming')) %>% 
  group_by(karyotype,donor,pcw,celltype) %>% summarise(nCell = n()) %>% 
  group_by(karyotype,donor,pcw) %>% mutate(totalCell = sum(nCell), frac = nCell/totalCell)
dd$celltype[dd$celltype == 'thy_Lumen-forming'] = 'fTFC2'
dd$celltype[dd$celltype == 'thy_TH_processing'] = 'fTFC1'
dd$pcw = factor(dd$pcw,c(11,12,14,15,17))
ggplot(dd,aes(karyotype,frac))+
  geom_point(aes(col=karyotype))+
  scale_color_manual(values = c('2n'='darkgreen','T21'='purple'))+
  #geom_boxplot(aes(fill=karyotype))+
  facet_grid(celltype~pcw,scales = 'free_y')+xlab('')+
  theme_classic()+
  theme(panel.border = element_rect(fill=F,colour = 'black'),
        axis.line = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(color='black'))


### check that only cells in G1 are retained
if(!keepCylcingCells){
  cyclingCells = srat$cellID[grepl('Cycling',srat$celltype)]
  if(length(cyclingCells) > 0){
    message(sprintf('Removing cycling cells from srat for tissue %s',tissue))
    srat = subset(srat, subset = cellID %in% srat$cellID[!srat$cellID %in% cyclingCells])
  }
}

srat$finalAnn = srat$celltype
srat$finalAnn[srat$finalAnn == 'thy_TH_processing'] = 'fTFC1'
srat$finalAnn[srat$finalAnn == 'thy_Lumen-forming'] = 'fTFC2'
srat$Genotype = srat$karyotype
srat$donorID = srat$donor
srat$age_group = gsub('-','_',srat$age_group)




## Using EDGER ##


##------------------##
##    Karyotype   ####
##------------------##

out = list()
pb_byCT_byAgeGroup = list()
# pb_outDir = '/lustre/scratch126/cellgen/team292/mt22/ThyroidDev/pb_T21vs2n'
# if(!dir.exists(pb_outDir)){
#   dir.create(pb_outDir,recursive = T)
# }
# 
# for(i in 1:length(pb_byCT_byAgeGroup)){
#   pb = pb_byCT_byAgeGroup[[i]]
#   write.csv(pb,file.path(pb_outDir,paste0(names(pb_byCT_byAgeGroup)[i],'.csv')),row.names = T)
# }
# 
# for(i in 1:length(pb_byCT)){
#   pb = pb_byCT[[i]]
#   write.csv(pb,file.path(pb_outDir,paste0(names(pb_byCT)[i],'.csv')),row.names = T)
# }

for(tgtCell in unique(srat$finalAnn)){
  for(ageGroup in unique(srat$age_group)){
    message(sprintf("\n\n------- Consider cell type %s from ageGroup %s",tgtCell,ageGroup))
    srat.sub = subset(srat,subset = cellID %in% srat$cellID[srat$finalAnn == tgtCell & srat$age_group == ageGroup])
    
    #Check we have a reasonable number of counts from both settings. ie Get number of cells for the relevant cell type and genenotype
    nCellsGroup = table(factor(srat.sub@meta.data$Genotype))
    
    if((!all(nCellsGroup>=50))){
      message(sprintf('Low number of cells detected'))
      print(nCellsGroup)
      next
    }
    
    #Check how many from individual donors
    nCells = table(srat.sub@meta.data$donorID)
    if(sum(nCells>50)<3){
      message(sprintf("Too few effective replicates.  Skipping..."))
      print(nCells)
      next
    }
    message("Found the following number of high confidence cells")
    print(nCells)
    
    
    # remove individuals with < 30 cells
    donorID_toKeep = names(nCells[nCells >= 30])
    if(tgtCell %in% c('fTFC1','fTFC2') & ageGroup =='14_20'){donorID_toKeep = donorID_toKeep[donorID_toKeep!='Hrv231']}
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
    w = which(as.character(seqnames(coords)) %in% tgtChrs)
    coords = coords[w]
    seqlevels(coords) = tgtChrs
    toc = toc[coords$gene_id,]
    #Non-expressed genes: these are genes expressed in <10 cells or <10% cells in all groups
    min_cells = 10
    w = which(rowSums(toc>0)>=min_cells)
    print(sprintf('Removing %d genes',nrow(toc) - length(w)))
    coords = coords[w]
    toc = toc[w,]
    # Remove genes expressed in < 10% of cells in at least 1 groups
    mDat$group = mDat[,'Genotype']
    percCell_expressed_perGene = do.call(cbind,lapply(split(colnames(toc),mDat[,'group']),function(e){
      tmp = data.frame(percCell_expressed = 100*rowSums(toc[,e,drop=FALSE] > 0) / length(e))
      colnames(tmp) = paste0(colnames(tmp),'_',unique(mDat$group[mDat$cellID %in% e]))
      return(tmp)
    } ))
    percCell_expressed_perGene$max_perc = apply(percCell_expressed_perGene,1,max)
    
    # only keep genes expressed >=10% cells in at least 1 group
    min_percCell = 10
    genes_toKeep = rownames(percCell_expressed_perGene[percCell_expressed_perGene$max_perc >=10,])
    print(sprintf('Removing %d genes',nrow(toc) - length(genes_toKeep)))
    coords = coords[genes_toKeep]
    toc = toc[genes_toKeep,]
    
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
    pb_byCT_byAgeGroup[[paste0(tgtCell,'_',ageGroup)]] = pb
    
    colDat = mDat[match(colnames(pb),mDat[,donorID]),]
    colDat$Genotype[colDat$Genotype == '2n'] = 'diploid'
    colDat$Genotype = factor(colDat$Genotype,c('diploid','T21'))
    rownames(colDat) = colDat[,donorID]
   
    
    ##=========================##
    # Fit EDGER model
    pdf(file.path(outDir,paste0(tgtCell,'_',ageGroup,'_pbEdgeR_plots.pdf')))
    out[[paste0(tgtCell,'_',ageGroup)]] = fit_model(pb,colDat,formula = '~ 0 + %s',geneMap,groupID='group')
    dev.off()
    
  }
}


saveRDS(out,file.path(outDir,'pb_edgeR_byCT_byAgeGroup.RDS'))


out = readRDS(file.path(outDir,'pb_edgeR_byCT_byAgeGroup.RDS'))
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
write.csv(allGenes,file.path(outDir,'allGenes_log2FC_splitbyCT.byAgeGroup.csv'))
#write.csv(allGenes,file.path(outDir,'allGenes_log2FC_splitbyCT.csv'))

allGenes = read.csv(file.path(outDir,'allGenes_log2FC_splitbyCT.byAgeGroup.csv'))
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

## Plot number of DEGs as function of number of cells
#-- Calculate number of cells in each pb comparison
nCell_perComp = data.frame()
for(tgtCell in unique(srat$finalAnn)){
  for(ageGroup in unique(srat$age_group)){
    message(sprintf("\n\n------- Consider cell type %s from ageGroup %s",tgtCell,ageGroup))
    srat.sub = subset(srat,subset = cellID %in% srat$cellID[srat$finalAnn == tgtCell & srat$age_group == ageGroup])
    
    #Check we have a reasonable number of counts from both settings. ie Get number of cells for the relevant cell type and genenotype
    nCellsGroup = table(factor(srat.sub@meta.data$Genotype))
    
    if((!all(nCellsGroup>=50))){
      message(sprintf('Low number of cells detected'))
      print(nCellsGroup)
      next
    }
    
    #Check how many from individual donors
    nCells = table(srat.sub@meta.data$donorID)
    if(sum(nCells>50)<3){
      message(sprintf("Too few effective replicates.  Skipping..."))
      print(nCells)
      next
    }
    message("Found the following number of high confidence cells")
    print(nCells)
    
    
    # remove individuals with < 30 cells
    donorID_toKeep = names(nCells[nCells >= 30])
    if(tgtCell %in% c('fTFC1','fTFC2') & ageGroup =='14_20'){donorID_toKeep = donorID_toKeep[donorID_toKeep!='Hrv231']}
    #OK, we're going ahead, create the needed objects
    toc = srat.sub@assays$RNA@counts[,colnames(srat.sub@assays$RNA@counts) %in% srat.sub$cellID[srat.sub$donorID %in% donorID_toKeep]]
    tmp = data.frame(nCell = ncol(toc),
                     ageGroup = ageGroup,
                     celltype = tgtCell)
    nCell_perComp = rbind(nCell_perComp,tmp)
  }
}

write.csv(nCell_perComp,file.path(outDir,'nCell_perComp.csv'))

nCell_perComp = read.csv(file.path(outDir,'nCell_perComp.csv'))


## add to the number of genes
df = rbind(df,data.frame(comp='fTFC2_11_13',direction=c('T21_up','T21_down'),nGene=0))
nCell_perComp$comp = paste0(nCell_perComp$celltype,'_',nCell_perComp$ageGroup)
df = cbind(df,nCell_perComp[match(df$comp,nCell_perComp$comp),!colnames(nCell_perComp) %in% colnames(df)])
df = df %>% group_by(comp,nCell,ageGroup,celltype) %>% summarise(nDEG = sum(abs(nGene)))

ggplot(df,aes(nCell,nDEG))+
  geom_line(aes(col=celltype))+
  geom_point(aes(col=celltype,shape=ageGroup),size=3)+
  scale_color_manual(values = col25)+
  xlab('Number of cells')+ylab('# DEGs')+
  theme_classic(base_size = 10)+xlim(0,16100)+
  theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),
        axis.line = element_blank(),#legend.position = 'bottom',
        axis.ticks = element_line(colour = 'black'),
        legend.title = element_text(size=11,colour = 'black'),
        legend.text = element_text(size=9,colour = 'black'),legend.key.size = unit(0.5,'cm'),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,colour = 'black'),
        axis.text = element_text(color='black'))


## Plot log2FC of age-specific genes
deg$celltype = gsub('_.*$','',deg$comp)
deg$ageGroup = gsub('fTFC1_|fTFC2_','',deg$comp)

dd = deg %>% group_by(celltype,direction,geneSym) %>% mutate(nAgeGroup = n_distinct(ageGroup)) %>% 
  group_by(ageGroup,direction,geneSym) %>% mutate(nCT = n_distinct(celltype))
dd$geneCat = ifelse(dd$nAgeGroup == 1 & dd$nCT == 1,dd$comp,
                    ifelse(dd$nAgeGroup == 1 & dd$nCT == 2,dd$ageGroup,
                           ifelse(dd$nAgeGroup == 2 & dd$nCT == 1,dd$celltype,'allShared')))
dd$geneCat = paste0(dd$direction,':',dd$geneCat)


## Results of Karyotype + ageGroup
out2 = readRDS(file.path(outDir,'pb_edgeR_byCT.RDS'))
deg2 = data.frame()
for(i in 1:length(out2)){
  tt = out2[[i]][['tt']]
  tt$comp = names(out2)[i]
  deg2 = rbind(deg2,tt)
}
deg2$direction = ifelse(deg2$logFC > 0,'T21_up','T21_down')

dd=rbind(dd[dd$direction == 'T21_down' & !dd$geneSym %in% deg2$geneSym[deg2$direction == 'T21_down'],],
         dd[dd$direction == 'T21_up' & !dd$geneSym %in% deg2$geneSym[deg2$direction == 'T21_up'],])

mtx = allGenes[allGenes$geneSym %in% dd$geneSym,]
mtx = pivot_wider(mtx,id_cols = c('geneSym'),names_from = 'comp',values_from = 'logFC')
mtx = column_to_rownames(mtx,'geneSym')
mtx = t(mtx)

col_fun = circlize::colorRamp2(c(-1, 0, 1), c('#1a4a87','white','#a4282c'))
Heatmap(mtx,show_column_dend = F,show_column_names = T,show_row_dend = F,show_row_names = T,
        cluster_rows = F,cluster_columns = T,
        column_title_rot = 90,column_title_gp = gpar(fontsize=7),column_names_gp = gpar(fontsize=7),column_names_rot = 90,
        split = gsub('_.*$','',rownames(mtx)),
        col=col_fun,
        column_split = dd$geneCat[match(colnames(mtx),dd$geneSym)])



df = plot_logCPM_byGroup(out2[[1]][['y']],genes = c(dd$geneSym[dd$geneCat == 'T21_up:fTFC1_14_20']),geneMap = geneMap)


genes = c('MT1X',"MT1H",'MC1R','NQO1','SULT1A2','FGF2','TOX2','RCAN1')
df = plot_logCPM_byGroup(out2[[2]][['y']],genes = c(dd$geneSym[dd$geneSym %in% genes]),geneMap = geneMap)

plotFun_pb_logCPM = function(noFrame=FALSE,noPlot=FALSE){
  df$geneSym = factor(df$geneSym,genes)
  p=ggplot(df,aes(group,logCPM,fill=Genotype))+
    geom_boxplot(outlier.shape = NA,width=0.45)+
    geom_jitter(size=0.7,width = 0.2)+
    scale_fill_manual(values = c('T21'='#81aa84','diploid'='#8062a4'))+
    facet_wrap(vars(geneSym),scales = 'free_y',ncol=5)+
    theme_classic(base_size = 13)+xlab('')+
    theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
          axis.line = element_blank(),axis.text = element_text(color='black'),
          strip.background = element_blank(),
          panel.border = element_rect(fill=F,colour = 'black'))
  print(p)
}
saveFig(file.path(plotDir,'Fig3f_pb_logCPM_fTFC2'),plotFun_pb_logCPM,rawData=df,width = 10,height = 5,res = 500,useDingbats = F)


## Plot overlapps in DEG
library(UpSetR)
deg_list = split(deg$geneSym,paste0(deg$comp,'_',deg$direction))
upset(fromList(deg_list),text.scale = 1.7,nsets = 20)

## 3 genes consistently upregulated in T21 TFC1

# Plot normalised expression in pb - logCPM
plot_logCPM_byGroup = function(y,genes,geneMap){
  df = y$logCPM[geneMap$ensID[geneMap$geneSym %in% genes],]
  rownames(df) = geneMap$geneSym[match(rownames(df),geneMap$ensID)]
  df = as.data.frame(df)
  df$geneSym = rownames(df)
  df = pivot_longer(df,col=colnames(df)[colnames(df) != 'geneSym'],names_to='donorID',values_to='logCPM')
  colDat = y$samples
  df = cbind(df,colDat[match(df$donorID,colDat$donorID),!colnames(colDat) %in% colnames(df)])
  df$group = paste0(df$age_group,':',df$Genotype)
  df$geneSym = factor(df$geneSym,unique(genes[genes %in% df$geneSym]))
  p=ggplot(df,aes(group,logCPM,fill=Genotype))+
    geom_boxplot(outlier.shape = NA,width=0.5)+
    geom_jitter(size=0.7,width = 0.25)+
    scale_fill_manual(values = c('T21'='#81aa84','diploid'='#8062a4'))+
    facet_wrap(vars(geneSym),scales = 'free_y')+
    theme_classic(base_size = 12)+xlab('')+
    theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
          axis.line = element_blank(),axis.text = element_text(color='black'),
          strip.background = element_blank(),
          panel.border = element_rect(fill=F,colour = 'black'))
  print(p)
  return(df)
}

df = plot_logCPM_byGroup(out[[1]][['y']],genes = c(geneMap$geneSym[geneMap$chr=='chr21' &
                                                              !geneMap$geneSym %in% deg$geneSym &
                                                              geneMap$ensID %in% rownames(out[[1]][['y']]$counts)][21:31],'SLC5A5'),geneMap = geneMap)

df = plot_logCPM_byGroup(out[[1]][['y']],genes = c('SLC5A5','COL6A1','DYRK1A','COL18A1','MT1H','MT1X','MT1G','FGD4','BMP8A'),geneMap = geneMap)

srat$group = paste0(srat$finalAnn,'_',srat$karyotype,'_',srat$age_group,'_',srat$donor)
DotPlot(srat,group.by = 'group',
        #features = intersect(intersect(deg_list[['fTFC2_14_20_T21_up']],deg_list[['fTFC1_14_20_T21_up']]),deg_list[['fTFC1_11_13_T21_up']])
        features = unique(c('SLC5A5','COL6A1','DYRK1A','COL18A1','MT1H','MT1X','MT1G','FGD4','BMP8A',
                     'TPO','SLC5A5','DUOXA2','PAX8','TG','TPO','DIO2','GLIS3','SLC26A4','TSHR'))
        )+RotatedAxis()+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        axis.text = element_text(size=10))+xlab('') + ylab('')




##-----------------------------##
##    Karyotype + ageGroup   ####
##-----------------------------##

out = list()
pb_byCT = list()

for(tgtCell in unique(srat$finalAnn)){

  message(sprintf("\n\n------- Consider cell type %s",tgtCell))
  srat.sub = subset(srat,subset = cellID %in% srat$cellID[srat$finalAnn == tgtCell])
  
  #Check we have a reasonable number of counts from both settings. ie Get number of cells for the relevant cell type and genenotype
  nCellsGroup = table(factor(srat.sub@meta.data$Genotype))
  
  if((!all(nCellsGroup>=50))){
    message(sprintf('Low number of cells detected'))
    print(nCellsGroup)
    next
  }
  
  #Check how many from individual donors
  nCells = table(srat.sub@meta.data$donorID)
  if(sum(nCells>50)<3){
    message(sprintf("Too few effective replicates.  Skipping..."))
    print(nCells)
    next
  }
  message("Found the following number of high confidence cells")
  print(nCells)
  
  
  # remove individuals with < 30 cells
  donorID_toKeep = names(nCells[nCells >= 30])
  if(tgtCell %in% c('fTFC1','fTFC2')){donorID_toKeep = donorID_toKeep[donorID_toKeep!='Hrv231']}
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
  w = which(as.character(seqnames(coords)) %in% tgtChrs)
  coords = coords[w]
  seqlevels(coords) = tgtChrs
  toc = toc[coords$gene_id,]
  
  #Non-expressed genes: these are genes expressed in <10 cells or <10% cells in all groups
  min_cells = 10
  w = which(rowSums(toc>0)>=min_cells)
  print(sprintf('Removing %d genes',nrow(toc) - length(w)))
  
  coords = coords[w]
  toc = toc[w,]
  
  
  # Remove genes expressed in < 10% of cells in at least 1 groups
  mDat$group = paste0(mDat[,'Genotype'],':',mDat$age_group)
  percCell_expressed_perGene = do.call(cbind,lapply(split(colnames(toc),mDat[,'group']),function(e){
    tmp = data.frame(percCell_expressed = 100*rowSums(toc[,e,drop=FALSE] > 0) / length(e))
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
  coords$TSS = ifelse(strand(coords)=='-',end(coords),start(coords))
  o = order(seqnames(coords),coords$TSS)
  coords = coords[o]
  toc = toc[o,]
  
  ##=========================##
  # Create pseudobulk 
  ## Split each donor into 3 pb
  ## mDat = mDat %>% group_by(donorID) %>% mutate(pb_ID = (seq(1:n())%%2))
  ## mDat$donorID = paste0(mDat$donorID,'_',mDat$pb_ID)
  pb = do.call(cbind,lapply(split(colnames(toc),mDat[,donorID]),function(e) rowSums(toc[,e,drop=FALSE])))
  pb_byCT[[tgtCell]] = pb
  colDat = mDat[match(colnames(pb),mDat[[donorID]]),]
  colDat$Genotype[colDat$Genotype == '2n'] = 'diploid'
  colDat$Genotype = factor(colDat$Genotype,c('diploid','T21'))
  rownames(colDat) = colDat[[donorID]]
  
  
  
  
  
  ##=========================##
  # Fit EDGER model
  pdf(file.path(outDir,paste0(tgtCell,'_pbEdgeR_plots.pdf')))
  out[[tgtCell]] = fit_model(pb,colDat,formula = '~ 0 + %s + age_group',geneMap,groupID='group',MDS_groups = c('Genotype','donorID','sex','age_group'))
  dev.off()
  
  
}



saveRDS(out,file.path(outDir,'pb_edgeR_byCT.RDS'))


##----  Gene set enrichment analysis with limma::camera() ####
## Obtain the gene sets
C2_genesets <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c2.all.v7.1.entrez.rds"))
HM_genesets <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.h.all.v7.1.entrez.rds")) 

library(org.Hs.eg.db)
entrez_geneid <- mapIds(org.Hs.eg.db, keys=rownames(fit), 
                        column="ENTREZID", 
                        keytype="ENSEMBL")
idx_cam <- ids2indices(c(C2_genesets,HM_genesets),id=entrez_geneid)

## Import the DEG object
out = readRDS(file.path(outDir,'pb_edgeR_byCT.RDS'))
fit = out[[1]][['fit']]
contrast = makeContrasts(groupT21-groupdiploid, levels=out[[1]][['y']]$design)  
fit<-glmQLFTest(fit, contrast=contrast)
summary(decideTests(fit))


## Run camera
cam_fTFC1 = limma::camera(out[[1]][['y']]$logCPM, index=idx_cam, design=out[[1]][['y']]$design, contrast=contrast) 
cam_fTFC1$group = 'fTFC1'
cam_fTFC2 <- limma::camera(out[[2]][['y']]$logCPM, index=idx_cam, design=out[[2]][['y']]$design, contrast=contrast) 
cam_fTFC2$group = 'fTFC2'

library(UpSetR)
upset(fromList(list('fTFC1_up' = rownames(cam_fTFC1)[cam_fTFC1$FDR < 0.05 & cam_fTFC1$Direction == 'Up'],
                    'fTFC1_down' = rownames(cam_fTFC1)[cam_fTFC1$FDR < 0.05 & cam_fTFC1$Direction == 'Down'],
                    'fTFC2_up' = rownames(cam_fTFC2)[cam_fTFC2$FDR < 0.05 & cam_fTFC2$Direction == 'Up'],
                    'fTFC2_down' = rownames(cam_fTFC2)[cam_fTFC2$FDR < 0.05 & cam_fTFC2$Direction == 'Down'])))

cam = cam[order(cam$PValue),]
cam <- cam[(cam$PValue < 0.05),] #  & (cam$FDR < 0.05)
C2HM_enrich <- rownames(cam)
cam$log10Pval <- -log10(cam$PValue)
cam$Pval_dir <- unlist(lapply(1:nrow(cam), function(x){
  if (cam[x,]$Direction == "Up"){
    cam[x,]$log10Pval
  } else {
    -cam[x,]$log10Pval
  }
}))



##----------------------------------------##
##  Assess overlaps between 2 analyses  ####
##----------------------------------------##
deg_analysis1 = deg
deg_analysis2 = deg
deg_combined = rbind(deg_analysis1,deg_analysis2)

write.csv(deg_combined,file.path(outDir,'combined_DEG_results.csv'))
deg_combined = read.csv(file.path(outDir,'combined_DEG_results.csv'),row.names = 1)

thyroid_genes = c('NKX2-1','FOXE1','PAX8','GLIS3','TG','TPO','TSHR','SLC5A5','DIO2','DUOXA2','IYD','SLC26A4')
table(thyroid_genes %in% deg_combined$geneSym)



colnames(deg_combined)[colnames(deg_combined) == 'comp'] = 'comparison'
deg_combined$comparison[deg_combined$comparison == 'fTFC1'] = 'analysis2_fTFC1'
deg_combined$comparison[deg_combined$comparison == 'fTFC2'] = 'analysis2_fTFC2'
write.csv(deg_combined,file.path(outDir,'combined_DEG_results_forHassan.csv'))
deg_combined = read.csv(file.path(outDir,'combined_DEG_results_forHassan.csv'))
library(UpSetR)
deg_list = split(deg_combined$geneSym,paste0(deg_combined$comp,'_',deg_combined$direction))
names(deg_list)[names(deg_list) == 'fTFC1_T21_up'] = 'analysis2_fTFC1_T21up'
names(deg_list)[names(deg_list) == 'fTFC1_11_13_T21_up'] = 'analysis1_fTFC1_11.13_T21up'
names(deg_list)[names(deg_list) == 'fTFC1_14_20_T21_up'] = 'analysis1_fTFC1_14.20_T21up'
upset(fromList(deg_list[grepl('fTFC1',names(deg_list)) & grepl('up',names(deg_list))]),text.scale = 1.9,nsets = 20,keep.order = T)

genes = deg_list[['fTFC1_T21_up']][!deg_list[['fTFC1_T21_up']] %in% c(deg_list[['fTFC1_11_13_T21_up']],
                                                                      deg_list[['fTFC1_14_20_T21_up']])]
genesL = deg_list[['fTFC1_14_20_T21_up']][!deg_list[['fTFC1_14_20_T21_up']] %in% c(deg_list[['fTFC1_11_13_T21_up']],
                                                                      deg_list[['fTFC1_T21_up']])]
genesE = deg_list[['fTFC1_11_13_T21_up']][!deg_list[['fTFC1_11_13_T21_up']] %in% c(deg_list[['fTFC1_14_20_T21_up']],
                                                                                   deg_list[['fTFC1_T21_up']])]
genes = c(genesL,genesE,'SLC5A5')

# Group C
genes = deg_list[['fTFC1_T21_up']][deg_list[['fTFC1_T21_up']] %in% deg_list[['fTFC1_11_13_T21_up']] &
                                     ! deg_list[['fTFC1_T21_up']] %in% deg_list[['fTFC1_14_20_T21_up']]]

genes = deg_list[['fTFC1_T21_up']][deg_list[['fTFC1_T21_up']] %in% intersect(deg_list[['fTFC1_11_13_T21_up']],
                                                                      deg_list[['fTFC1_14_20_T21_up']])]


genes = deg_list[['fTFC2_T21_up']][!deg_list[['fTFC2_T21_up']] %in% c(deg_list[['fTFC2_11_13_T21_up']],
                                                                             deg_list[['fTFC2_14_20_T21_up']])]

length(genes)
df = plot_logCPM_byGroup(out[[1]][['y']],genes = genes[1:16],geneMap = geneMap)





srat$group = paste0(srat$finalAnn,'_',srat$karyotype,'_',srat$age_group)
DotPlot(srat,group.by = 'group',
        #features = intersect(intersect(deg_list[['fTFC2_14_20_T21_up']],deg_list[['fTFC1_14_20_T21_up']]),deg_list[['fTFC1_11_13_T21_up']])
        features = c('SLC5A5','COL6A1','DYRK1A','COL18A1','MT1H','MT1X','MT1G','FGD4','BMP8A')
)




##----------------------------------------------------##
##    Karyotype + age_group + karyotype:age_group   ####
##----------------------------------------------------##
out = list()

for(tgtCell in unique(srat$finalAnn)){
  
  message(sprintf("\n\n------- Consider cell type %s",tgtCell))
  srat.sub = subset(srat,subset = cellID %in% srat$cellID[srat$finalAnn == tgtCell])
  
  #Check we have a reasonable number of counts from both settings. ie Get number of cells for the relevant cell type and genenotype
  nCellsGroup = table(factor(srat.sub@meta.data$Genotype))
  
  if((!all(nCellsGroup>=50))){
    message(sprintf('Low number of cells detected'))
    print(nCellsGroup)
    next
  }
  
  #Check how many from individual donors
  nCells = table(srat.sub@meta.data$donorID)
  if(sum(nCells>50)<3){
    message(sprintf("Too few effective replicates.  Skipping..."))
    print(nCells)
    next
  }
  message("Found the following number of high confidence cells")
  print(nCells)
  
  
  # remove individuals with < 30 cells
  donorID_toKeep = names(nCells[nCells >= 30])
  if(tgtCell %in% c('fTFC1','fTFC2')){donorID_toKeep = donorID_toKeep[donorID_toKeep!='Hrv231']}
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
  colDat$Genotype[colDat$Genotype == '2n'] = 'diploid'
  colDat$Genotype = factor(colDat$Genotype,c('diploid','T21'))
  rownames(colDat) = colDat[,donorID]
  
  
  
  
  
  ##=========================##
  # Fit EDGER model
  pdf(file.path(outDir,paste0(tgtCell,'_pbEdgeR_plots.pdf')))
  out[[tgtCell]] = fit_model(pb,colDat,formula = '~ %s + age_group + group:age_group',geneMap,groupID='group',MDS_groups = c('Genotype','donorID','sex','age_group'),coef = c(2))
  dev.off()
  
  
}



saveRDS(out,file.path(outDir,'pb_edgeR_byCT.RDS'))
out = readRDS(file.path(outDir,'pb_edgeR_byCT.RDS'))










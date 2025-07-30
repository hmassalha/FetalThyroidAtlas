## Statistical test to assess the increase in proportion of mesenchymal and X in mid and late developmental stages compared to early stage


##----------------##
##   Libraries  ####
##----------------##
library(tidyverse)
library(Seurat)
library(GenomicFeatures)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
source("~/lustre_mt22/generalScripts/utils/pseudobulk.R")

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

fThy = readRDS('~/lustre_mt22/Thyroid/Data/fetalThyroid/fThyroid_annotated_fromHassan_jul23.RDS')
table(fThy$cluster,fThy$celltype)
DimPlot(fThy,group.by = 'celltype',cols = c(col25,pal34H),label = T,repel = T,label.box = T,label.size = 3) + NoAxes() + NoLegend()
DimPlot(fThy,cells.highlight = fThy$cellID[fThy$celltype == 'pat_Parathyrocytes'])
df = fThy@meta.data[,c('cellID','donor','pcw','cluster','celltype')]
df$age_group = ifelse(df$pcw %in% c(9,10),'9-10',
                      ifelse(df$pcw %in% c(11,12,13),'11-13','14-20'))
df$celltype_group = '?'
df$celltype_group[df$celltype == 'thy_TH_processing'] = 'fTFC1'
df$celltype_group[df$celltype == 'thy_Lumen-forming'] = 'fTFC2'
df$celltype_group[df$celltype == 'smc_Pericytes'] = 'Pericytes'
df$celltype_group[df$celltype %in% c('mes_CYGB','mes_IGF1R')] = 'fTS1'
df$celltype_group[df$celltype %in% c('mes_KCNB2')] = 'fTS2'
df$celltype_group[df$celltype %in% c('pat_Parathyrocytes')] = 'Parathyroid_cells'
df$celltype_group[df$celltype %in% c('paf_C_cells')] = 'C_cells'

table(df$celltype_group)
df = df[df$celltype_group != '?',]
df2 = df %>% group_by(celltype_group,age_group) %>% summarise(n=n()) %>% 
  group_by(age_group) %>% mutate(totalCell = sum(n),perc = 100*n/totalCell)
df2$age_group = factor(df2$age_group,c('9-10','11-13','14-20'))
df2$celltype_group = factor(df2$celltype_group,rev(c('C_cells','Parathyroid_cells','fTS1','fTS2','Pericytes','fTFC2','fTFC1')))
ggplot(df2,aes(age_group,y=perc,fill=celltype_group))+
  geom_col()+
  theme_classic()+
  scale_fill_manual(values = col25)+
  theme(panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank())



## Thoughts on what test to use:
# if used two-proportion test, this is a fancy version of binomial test (perhaps when n is large enough? need to read up on the differences a little more)
# but it is design to compare the proportion between two groups. 
# With multiple samplings per group (eg. different donorID per group), 
# unless we expect a systematic difference (then would be more appropriate to perhaps do chi-squared), 
# we can treat different samplings as if it has come from the same big group --> n_cells_interest_across_ALL_donors / total_cells_across_ALL_donors
# http://www.sthda.com/english/wiki/two-proportions-z-test-in-r
# Also, a note on the Yate's continuity correction: https://library.virginia.edu/data/articles/continuity-corrections-imperfect-responses-to-slight-problems
# 
# Alternatively, if wanted to test for the "mean proportion" per group --> then can do a standard wilcoxon / t-test

## Do proportion t-test (z-test)
dd = df %>% group_by(celltype_group,age_group) %>% summarise(n=n()) %>% 
  group_by(age_group) %>% mutate(totalCell = sum(n),perc = 100*n/totalCell)

##--- on fTS cells, not grouped by donorID
prop.test(x=c(dd$n[dd$celltype_group == 'fTS1' & dd$age_group == '11-13'],
              dd$n[dd$celltype_group == 'fTS1' & dd$age_group == '9-10']), 
          n=c(unique(dd$totalCell[dd$celltype_group == 'fTS1' & dd$age_group == '11-13']),
              unique(dd$totalCell[dd$celltype_group == 'fTS1' & dd$age_group == '9-10'])), 
          p = NULL, alternative = "greater",
          correct = FALSE)

prop.test(x=c(dd$n[dd$celltype_group == 'fTS1' & dd$age_group == '14-20'],
              dd$n[dd$celltype_group == 'fTS1' & dd$age_group == '9-10']), 
          n=c(unique(dd$totalCell[dd$celltype_group == 'fTS1' & dd$age_group == '14-20']),
              unique(dd$totalCell[dd$celltype_group == 'fTS1' & dd$age_group == '9-10'])), 
          p = NULL, alternative = "greater",
          correct = FALSE)


##--- on pericytes cells, not grouped by donorID
prop.test(x=c(dd$n[dd$celltype_group == 'Pericytes' & dd$age_group == '11-13'],
              dd$n[dd$celltype_group == 'Pericytes' & dd$age_group == '9-10']), 
          n=c(unique(dd$totalCell[dd$celltype_group == 'Pericytes' & dd$age_group == '11-13']),
              unique(dd$totalCell[dd$celltype_group == 'Pericytes' & dd$age_group == '9-10'])), 
          p = NULL, alternative = "greater",
          correct = FALSE)

prop.test(x=c(dd$n[dd$celltype_group == 'Pericytes' & dd$age_group == '14-20'],
              dd$n[dd$celltype_group == 'Pericytes' & dd$age_group == '9-10']), 
          n=c(unique(dd$totalCell[dd$celltype_group == 'Pericytes' & dd$age_group == '14-20']),
              unique(dd$totalCell[dd$celltype_group == 'Pericytes' & dd$age_group == '9-10'])), 
          p = NULL, alternative = "greater",
          correct = FALSE)

## Do normal wilcoxon test
dd = df %>% group_by(donor,celltype_group,age_group) %>% summarise(n=n()) %>% 
  group_by(donor,age_group) %>% mutate(totalCell = sum(n),perc = 100*n/totalCell)
dd$age_group = factor(dd$age_group,c('9-10','11-13','14-20'))
library(ggbeeswarm)
ggplot(dd[dd$celltype_group %in% c('fTS1','fTS2','Pericytes','fTFC1','fTFC2','Parathyroid_cells'),],aes(age_group,perc))+
  geom_boxplot(width=0.4,aes(fill=age_group))+
  facet_wrap(vars(celltype_group),nrow=1,scales = 'free_y')+
  geom_quasirandom(size=1.3,width = 0.3)+
  theme_classic()+
  theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),
        axis.line = element_blank())


wilcox.test(dd$perc[dd$celltype_group == 'fTS1' & dd$age_group == '11-13'],
            dd$perc[dd$celltype_group == 'fTS1' & dd$age_group == '9-10'],alternative = 'greater')
wilcox.test(dd$perc[dd$celltype_group == 'fTS1' & dd$age_group == '14-20'],
            dd$perc[dd$celltype_group == 'fTS1' & dd$age_group == '9-10'],alternative = 'greater')


wilcox.test(dd$perc[dd$celltype_group == 'Pericytes' & dd$age_group == '11-13'],
            dd$perc[dd$celltype_group == 'Pericytes' & dd$age_group == '9-10'],alternative = 'greater')

wilcox.test(dd$perc[dd$celltype_group == 'Pericytes' & dd$age_group == '14-20'],
            dd$perc[dd$celltype_group == 'Pericytes' & dd$age_group == '9-10'],alternative = 'greater')





## Do normal wilcoxon test but group mid and late together
df$age_group2 = ifelse(df$age_group != '9-10','late','early')
dd = df %>% group_by(donor,celltype_group,age_group2) %>% summarise(n=n()) %>% 
  group_by(donor,age_group2) %>% mutate(totalCell = sum(n),perc = 100*n/totalCell)
dd$age_group2 = factor(dd$age_group2,c('early','late'))
library(ggbeeswarm)
ggplot(dd[dd$celltype_group %in% c('fTS1','fTS2','Pericytes','fTFC1','fTFC2','Parathyroid_cells'),],aes(age_group2,perc))+
  geom_boxplot(width=0.4,aes(fill=age_group2))+
  facet_wrap(vars(celltype_group),nrow=1,scales = 'free_y')+
  geom_quasirandom(size=1.3,width = 0.3)+
  theme_classic()+
  theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),
        axis.line = element_blank())

dd$age_group = dd$age_group2
wilcox.test(dd$perc[dd$celltype_group == 'fTS1' & dd$age_group == 'late'],
            dd$perc[dd$celltype_group == 'fTS1' & dd$age_group == 'early'],alternative = 'greater')

wilcox.test(dd$perc[dd$celltype_group == 'Pericytes' & dd$age_group == 'late'],
            dd$perc[dd$celltype_group == 'Pericytes' & dd$age_group == 'early'],alternative = 'greater')




##--------   Ratio of fTFC1 : fTFC2   ---------##
## Across all donors within the age group
dd = fThy@meta.data[fThy$celltype %in% c('thy_Cycling','thy_Lumen-forming','thy_TH_processing'),c('cellID','donor','pcw','cluster','celltype')]
dd$age_group = df$age_group[match(dd$donor,df$donor)]
dd = dd %>% group_by(celltype,age_group) %>% summarise(n=n()) %>% 
  group_by(age_group) %>% mutate(totalCell = sum(n),perc = 100*n/totalCell)
dd = dd[dd$celltype_group %in% c('fTFC1','fTFC2'),]
dd = pivot_wider(dd,id_cols = c('age_group'),values_from = 'n',names_from = 'celltype_group')
dd$ratio = dd$fTFC1 / dd$fTFC2

## Across all donors within the age group
dd = df %>% group_by(donor,celltype_group,age_group) %>% summarise(n=n()) %>% 
  group_by(donor,age_group2) %>% mutate(totalCell = sum(n),perc = 100*n/totalCell)









outDir = "~/lustre_mt22/Thyroid/Results/x_2nT21_edgeR_pseudobulk"
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

setwd(outDir)



##-----------------------------------------------##
##    Checking the pseudobulk normalisation    ####
##-----------------------------------------------##
## Try creating the pseudobulks myself


ggRLE <- function(dat_x, annot, col_str,isLog=TRUE, isLarge=FALSE,
                  ylim = c(-2,2),zero_line=TRUE, zero_col="skyblue", medPoint=FALSE, whisk=1.5){
  
  col_title <- str_to_title(col_str)
  
  if (!(setequal(colnames(dat_x), rownames(annot)))){
    message("Make sure annot rownames match dat_x sample names")
  }
  
  if (!(isLog)){
    dat_x <- log2(dat_x+1)
  }
  
  ## Subset to common sample names
  annot <- annot[order(as.vector(annot[[col_str]])), , drop = F]
  annot$ColourBy <- as.vector(annot[[col_str]])
  annot$Sample <- rownames(annot)
  dat_x <- dat_x[, rownames(annot)]
  
  ## RLE boxplots
  rle <- dat_x - rowMedians(dat_x)
  rleLong <- reshape2::melt(rle, value.name = "RLE", varnames = c("genes", "Sample"))
  rleLong$Sample <- as.character(rleLong$Sample)
  
  rleLong <- merge(data.table::data.table(varhandle::unfactor(rleLong)),  # faster merging
                   data.table::data.table(varhandle::unfactor(annot)),
                   by = "Sample", sort=F)
  
  ## Calculate the median of the RLE boxplots:
  rleLong <- rleLong %>%
    group_by(Sample) %>%
    mutate(MedRLE = median(RLE)) %>%
    ungroup() %>%
    data.frame()
  rleLong$Sample <- factor(rleLong$Sample , levels=unique(rleLong$Sample))
  
  ## Remove Whiskers of the boxplots if the sample size is very large -- Sep's code
  if(isLarge){
    whisk=0
  }
  
  if (medPoint){
    gg <- ggplot(rleLong, aes(x = Sample, y = RLE, fill = ColourBy))+
      stat_boxplot(geom = "errorbar", width = 0.3)+
      geom_boxplot(outlier.shape = NA, coef=whisk)+
      geom_point(data = rleLong[! duplicated(rleLong$Sample), ],
                 aes(x = Sample, y = MedRLE, fill = ColourBy),
                 size = 2, shape = 21, colour = "black", lwd = 2)+
      scale_fill_manual(name=col_title)+
      scale_y_continuous(name = "RLE",limits = ylim)+theme_bw()+
      theme(panel.grid.minor = element_blank(),
            #  panel.grid.major = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank())
  } else {
    gg <- ggplot(rleLong, aes(x = Sample, y = RLE, fill = ColourBy))+
      stat_boxplot(geom = "errorbar", width = 0.3)+
      geom_boxplot(outlier.shape = NA,coef=whisk)+
      #scale_fill_manual(name=col_title)+
      scale_y_continuous(name = "RLE",limits = ylim)+theme_bw()+
      theme(panel.grid.minor = element_blank(),
            #  panel.grid.major = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank())
  }
  
  if (zero_line){
    gg <- gg + geom_hline(yintercept = 0, col = zero_col, lwd = 1)+
      theme(panel.grid.major=element_blank())
  }
  return(gg)
}




## Hassans pseudobulks
pb_hm11 = Matrix::readMM('~/lustre_mt22/Thyroid/Results/x_2nT21_edgeR_pseudobulk/matrix.mtx.gz')
pb_hm11 = as.matrix(pb_hm11)
obs = read.csv('~/lustre_mt22/Thyroid/Results/x_2nT21_edgeR_pseudobulk/metadata.csv')
colnames(pb_hm11) = obs$X
var = read.delim('~/lustre_mt22/Thyroid/Results/x_2nT21_edgeR_pseudobulk/geneMap.csv',header = F)
rownames(pb_hm11) = var$V1


pb_donorMerged = do.call(cbind,lapply(split(obs$X,obs$replicate),function(e){
  pb = pb_hm11[,e]
  pb = rowSums(pb)
}))

pb_donorMerged2 = Matrix::readMM('~/lustre_mt22/Thyroid/Results/x_2nT21_edgeR_pseudobulk/noTechRep/matrix.mtx.gz')
pb_donorMerged2 = as.matrix(pb_donorMerged2)
obs = read.csv('~/lustre_mt22/Thyroid/Results/x_2nT21_edgeR_pseudobulk/noTechRep/metadata.csv')
colnames(pb_donorMerged2) = obs$X
obs$age_group = ifelse(obs$pcw %in% c(11:13),'11_13','14_20')
var = read.delim('~/lustre_mt22/Thyroid/Results/x_2nT21_edgeR_pseudobulk/noTechRep/geneMap.csv',header = F)
rownames(pb_donorMerged2) = var$V1

pb_donorMerged2 = pb
obs = sampDat
obs$karyotype = obs$group
obs$age_group = obs$ageGroup
obs$celltype = 'fTFC1'
#obs = obs[match(colnames(pb_donorMerged),obs$replicate),]
# create an edgeR object with counts and grouping factor
y <- DGEList(pb_donorMerged2, group = obs$karyotype)
y <- DGEList(pb_hm11, group = obs$karyotype)




edgeR_fitModel = function(pb,obs,group,design){
  library(edgeR)  
  # create an edgeR object with counts and grouping factor
  y <- DGEList(pb, group = obs[[group]])
  
  # filter out genes with low counts
  print("Dimensions before subsetting:")
  print(dim(y))
  print("")
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  print("Dimensions after subsetting:")
  print(dim(y))
  print("")
  
  # normalize
  y <- calcNormFactors(y)
  
  # estimate dispersion
  y <- estimateDisp(y, design = design)
  
  # fit the model
  fit <- glmQLFit(y, design,robust = T)
  
  return(list("fit"=fit, "design"=design, "y"=y))
}


sampDat$karyotype = sampDat$group
sampDat$age_group = sampDat$ageGroup

pb_fTFC1 = pb
obs_fTFC1 = sampDat[match(colnames(pb_fTFC1),sampDat$donor),]

obs = obs_fTFC1



# pb_fTFC2 = pb
# obs_fTFC1 = sampDat[match(colnames(pb_fTFC1),sampDat$donor),]
# 
# obs = obs_fTFC1
obs$karyotype[obs$karyotype == '2n'] = 'diploid'

##--------------------------##
##    ~ 0 + karyotype     ####
##--------------------------##

# build design matrix      
karyotype <- as.factor(obs[['karyotype']])
print(levels(karyotype))
age_group <- as.factor(gsub('-','_',obs$age_group))
print(levels(age_group))
donor <- as.factor(obs$donor)
print(levels(donor))

design <- model.matrix(~ 0 + karyotype + age_group)

## Fit the model
out = edgeR_fitModel(pb=pb_fTFC1,obs=obs,group='karyotype',design=design)
y = out[['y']]
fit = out[['fit']]
plotMDS(y, col=ifelse(y$samples$group == "diploid", "red", "blue"))
plotBCV(y)
colnames(y$design)

myContrast <- makeContrasts('karyotypeT21-karyotypediploid', levels = y$design)
qlf <- glmQLFTest(fit, contrast=myContrast)

# get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
tt <- topTags(qlf, n = Inf)
tt <- tt$table

plotSmear(qlf, de.tags = rownames(tt)[which(tt$FDR<0.05)])

table(tt$FDR < 0.05)
# tt = tt[!is.na(geneMap$ensID[match(rownames(tt),geneMap$ensID)]),]
# rownames(tt) = geneMap$ensID[match(rownames(tt),geneMap$geneSym)]
tt = annotateGenes(tt,geneMap = geneMap)
tt$DE = (tt$FDR < 0.05 & abs(tt$logFC) > 0.5 )
table(tt$DE)
tt$chr = factor(tt$chr,c(paste0('chr',c(1:22)),'chrM','chrX','chrY',unique(tt$chr[grepl('GL\\d+|KI\\d+',tt$chr)])))
library(ggbeeswarm)
ggplot(tt,aes(chr,logFC))+
  geom_quasirandom(aes(col=DE),size=0.1)+
  geom_boxplot(outlier.shape = NA,alpha=0.5)+
  scale_color_manual(values = c('black','red'))+
  geom_hline(yintercept = c(0,0.5,-0.5))+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank())+
  xlab('') + ggtitle('T21 vs 2n',subtitle = '1 replicate per donor / ~ 0 + karyotype')




##-----------------------------------##
##    ~ 0 + karyotype_ageGroup     ####
##-----------------------------------##

# build design matrix     
obs$group = paste0(obs$karyotype,'_',gsub('-','_',obs$age_group))
group <- as.factor(obs[['group']])
print(levels(group))
# karyotype <- as.factor(obs[['karyotype']])
# print(levels(karyotype))
# age_group <- as.factor(gsub('-','_',obs$age_group))
# print(levels(age_group))
# donor <- as.factor(obs$donor)
# print(levels(donor))

design <- model.matrix(~ 0 + group)

## Fit the model
out = edgeR_fitModel(pb=pb_fTFC1,obs=obs,group='group',design=design)
y = out[['y']]
fit = out[['fit']]
plotMDS(y, col=ifelse(y$samples$group == "diploid", "red", "blue"))
plotBCV(y)
colnames(y$design)

myContrast <- makeContrasts('groupT21_14_20-group2n_14_20', levels = y$design)
qlf <- glmQLFTest(fit, contrast=myContrast)

# get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
tt <- topTags(qlf, n = Inf)
tt <- tt$table

plotSmear(qlf, de.tags = rownames(tt)[which(tt$FDR<0.05)])

table(tt$FDR < 0.05 & tt$PValue < 0.05)
# tt = tt[!is.na(geneMap$ensID[match(rownames(tt),geneMap$ensID)]),]
# rownames(tt) = geneMap$ensID[match(rownames(tt),geneMap$geneSym)]
tt = annotateGenes(tt,geneMap = geneMap)
tt$DE = (tt$FDR < 0.05 & abs(tt$logFC) > 0.5 )
table(tt$DE)
tt$chr = factor(tt$chr,c(paste0('chr',c(1:22)),'chrM','chrX','chrY',unique(tt$chr[grepl('GL\\d+|KI\\d+',tt$chr)])))
library(ggbeeswarm)
ggplot(tt,aes(chr,logFC))+
  geom_quasirandom(aes(col=DE),size=0.1)+
  geom_boxplot(outlier.shape = NA,alpha=0.5)+
  scale_color_manual(values = c('black','red'))+
  geom_hline(yintercept = c(0,0.5,-0.5))+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank())+
  xlab('') + ggtitle('T21 vs 2n',subtitle = '1 replicate per donor / ~ 0 + karyotype')



##-------------------------------------##
##    ~ 0 + karyotype + ageGroup     ####
##-------------------------------------##

# build design matrix     
karyotype <- as.factor(obs[['karyotype']])
print(levels(karyotype))
age_group <- as.factor(gsub('-','_',obs$age_group))
print(levels(age_group))
# donor <- as.factor(obs$donor)
# print(levels(donor))

design <- model.matrix(~ 0 + age_group + karyotype)

## Fit the model
out = edgeR_fitModel(pb=pb_fTFC1,obs=obs,group='karyotype',design=design)
y = out[['y']]
fit = out[['fit']]
plotMDS(y, col=ifelse(y$samples$group == "diploid", "red", "blue"))
plotBCV(y)
colnames(y$design)

myContrast <- makeContrasts('age_group14_20 - age_group11_13', levels = y$design)
myContrast <- makeContrasts('karyotypeT21', levels = y$design)
qlf <- glmQLFTest(fit, contrast=myContrast)

# get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
tt <- topTags(qlf, n = Inf)
tt <- tt$table

plotSmear(qlf, de.tags = rownames(tt)[which(tt$FDR<0.05)])

table(tt$FDR < 0.05 & tt$PValue < 0.05)
# tt = tt[!is.na(geneMap$ensID[match(rownames(tt),geneMap$ensID)]),]
# rownames(tt) = geneMap$ensID[match(rownames(tt),geneMap$geneSym)]
tt = annotateGenes(tt,geneMap = geneMap)
tt$DE = (tt$FDR < 0.05 & abs(tt$logFC) > 0.5 )
tt$chr = factor(tt$chr,c(paste0('chr',c(1:22)),'chrM','chrX','chrY',unique(tt$chr[grepl('GL\\d+|KI\\d+',tt$chr)])))
library(ggbeeswarm)
ggplot(tt,aes(chr,logFC))+
  geom_quasirandom(aes(col=DE),size=0.1)+
  geom_boxplot(outlier.shape = NA,alpha=0.5)+
  scale_color_manual(values = c('black','red'))+
  geom_hline(yintercept = c(0,0.5,-0.5))+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank())+
  xlab('') + ggtitle('T21 vs 2n',subtitle = '1 replicate per donor / ~ 0 + karyotype')



##----------------------------------------------------------##
##    ~ 0 + karyotype + ageGroup + ageGroup:karyotype     ####
##----------------------------------------------------------##

# build design matrix     
karyotype <- as.factor(obs[['karyotype']])
print(levels(karyotype))
age_group <- as.factor(gsub('-','_',obs$age_group))
print(levels(age_group))
# donor <- as.factor(obs$donor)
# print(levels(donor))

design <- model.matrix(~ 0 + age_group + karyotype + age_group*karyotype)

## Fit the model
out = edgeR_fitModel(pb=pb_fTFC1,obs=obs,group='karyotype',design=design)
y = out[['y']]
fit = out[['fit']]
plotMDS(y, col=ifelse(y$samples$group == "diploid", "red", "blue"))
plotBCV(y)
colnames(y$design)

myContrast <- makeContrasts("age_group14_20:karyotypeT21", levels = y$design)

qlf <- glmQLFTest(fit, contrast=myContrast)
qlf <- glmQLFTest(fit, coef=1)

# get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
tt <- topTags(qlf, n = Inf)
tt <- tt$table

plotSmear(qlf, de.tags = rownames(tt)[which(tt$FDR<0.05)])

table(tt$FDR < 0.05 & tt$PValue < 0.05)
# tt = tt[!is.na(geneMap$ensID[match(rownames(tt),geneMap$ensID)]),]
# rownames(tt) = geneMap$ensID[match(rownames(tt),geneMap$geneSym)]
tt = annotateGenes(tt,geneMap = geneMap)
tt$DE = (tt$FDR < 0.05 & abs(tt$logFC) > 0.5 )
tt$chr = factor(tt$chr,c(paste0('chr',c(1:22)),'chrM','chrX','chrY',unique(tt$chr[grepl('GL\\d+|KI\\d+',tt$chr)])))
library(ggbeeswarm)
ggplot(tt,aes(chr,logFC))+
  geom_quasirandom(aes(col=DE),size=0.1)+
  geom_boxplot(outlier.shape = NA,alpha=0.5)+
  scale_color_manual(values = c('black','red'))+
  geom_hline(yintercept = c(0,0.5,-0.5))+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank())+
  xlab('') + ggtitle('T21 vs 2n',subtitle = '1 replicate per donor / ~ 0 + karyotype')






##-----------------------------------------------##
##    ~ 0 + karyotype (separate age groups)    ####
##-----------------------------------------------##

##--- young age group
obs = obs_fTFC1[obs_fTFC1$age_group == '11-13',]
pb_fTFC1_young = pb_fTFC1[,obs$donor]
dim(pb_fTFC1_young)

# build design matrix     
karyotype <- as.factor(obs[['karyotype']])
print(levels(karyotype))
# age_group <- as.factor(gsub('-','_',obs$age_group))
# print(levels(age_group))
# donor <- as.factor(obs$donor)
# print(levels(donor))

design <- model.matrix(~ 0 + karyotype)
colnames(design)

## Fit the model
out = edgeR_fitModel(pb=pb_fTFC1_young,obs=obs,group='karyotype',design=design)
y = out[['y']]
fit = out[['fit']]
plotMDS(y, col=ifelse(y$samples$group == "diploid", "red", "blue"))
plotBCV(y)
colnames(y$design)

myContrast <- makeContrasts("karyotypeT21-karyotypediploid", levels = y$design)
qlf <- glmQLFTest(fit, contrast=myContrast)
# qlf <- glmQLFTest(fit, coef=4)

# get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
tt <- topTags(qlf, n = Inf)
tt <- tt$table

plotSmear(qlf, de.tags = rownames(tt)[which(tt$FDR<0.05)])

table(tt$FDR < 0.05)
# tt = tt[!is.na(geneMap$ensID[match(rownames(tt),geneMap$ensID)]),]
# rownames(tt) = geneMap$ensID[match(rownames(tt),geneMap$geneSym)]
tt = annotateGenes(tt,geneMap = geneMap)
tt$DE = (tt$FDR < 0.05 & abs(tt$logFC) > 0.2 )
table(tt$DE)
tt$chr = factor(tt$chr,c(paste0('chr',c(1:22)),'chrM','chrX','chrY',unique(tt$chr[grepl('GL\\d+|KI\\d+',tt$chr)])))
library(ggbeeswarm)
ggplot(tt,aes(chr,logFC))+
  geom_quasirandom(aes(col=DE),size=0.1)+
  geom_boxplot(outlier.shape = NA,alpha=0.5)+
  scale_color_manual(values = c('black','red'))+
  geom_hline(yintercept = c(0,0.5,-0.5))+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank())+
  xlab('') + ggtitle('T21 vs 2n',subtitle = '1 replicate per donor / ~ 0 + karyotype')




##--- OLD age group
obs = obs_fTFC1[obs_fTFC1$age_group == '14-20',]
pb_fTFC1_young = pb_fTFC1[,obs$donor]
dim(pb_fTFC1_young)

# build design matrix     
karyotype <- as.factor(obs[['karyotype']])
print(levels(karyotype))
# age_group <- as.factor(gsub('-','_',obs$age_group))
# print(levels(age_group))
# donor <- as.factor(obs$donor)
# print(levels(donor))

design <- model.matrix(~0+ karyotype)
colnames(design)

## Fit the model
out = edgeR_fitModel(pb=pb_fTFC1_young,obs=obs,group='karyotype',design=design)
y = out[['y']]
fit = out[['fit']]
plotMDS(y, col=ifelse(y$samples$group == "diploid", "red", "blue"))
plotBCV(y)
colnames(y$design)

myContrast <- makeContrasts("karyotypeT21-karyotypediploid", levels = y$design)
qlf <- glmQLFTest(fit, contrast=myContrast)
#qlf <- glmQLFTest(fit, coef=2)

# get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
tt <- topTags(qlf, n = Inf)
tt <- tt$table

plotSmear(qlf, de.tags = rownames(tt)[which(tt$FDR<0.05)])

table(tt$FDR < 0.05)
# tt = tt[!is.na(geneMap$ensID[match(rownames(tt),geneMap$ensID)]),]
# rownames(tt) = geneMap$ensID[match(rownames(tt),geneMap$geneSym)]
tt = annotateGenes(tt,geneMap = geneMap)
tt$DE = (tt$FDR < 0.05 & abs(tt$logFC) > 0.2 )
table(tt$DE)
tt$chr = factor(tt$chr,c(paste0('chr',c(1:22)),'chrM','chrX','chrY',unique(tt$chr[grepl('GL\\d+|KI\\d+',tt$chr)])))
library(ggbeeswarm)
ggplot(tt,aes(chr,logFC))+
  geom_quasirandom(aes(col=DE),size=0.1)+
  geom_boxplot(outlier.shape = NA,alpha=0.5)+
  scale_color_manual(values = c('black','red'))+
  geom_hline(yintercept = c(0,0.5,-0.5))+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank())+
  xlab('') + ggtitle('T21 vs 2n',subtitle = '1 replicate per donor / ~ 0 + karyotype')




##------------------------------------------------------##
##    ~ 0 + karyotype (separate age groups) - fTFC2   ####
##------------------------------------------------------##
#obs_fTFC2 = obs
obs_fTFC2$karyotype[obs_fTFC2$karyotype == '2n'] = 'diploid'
pb_fTFC2 = pb_donorMerged2
colnames(pb_fTFC2) = gsub('donor_|_0','',colnames(pb_fTFC2))
#rownames(pb_fTFC2) = geneMap$ensID[match(rownames(pb_fTFC2),geneMap$geneSym)]
##--- young age group
obs = obs_fTFC2[obs_fTFC2$age_group == '11_13',]
pb_fTFC2_young = pb_fTFC2[,obs$donor]
dim(pb_fTFC2_young)

# build design matrix     
karyotype <- as.factor(obs[['karyotype']])
print(levels(karyotype))
# age_group <- as.factor(gsub('-','_',obs$age_group))
# print(levels(age_group))
# donor <- as.factor(obs$donor)
# print(levels(donor))

design <- model.matrix(~ 0 + karyotype)
colnames(design)

## Fit the model
out = edgeR_fitModel(pb=pb_fTFC2_young,obs=obs,group='karyotype',design=design)
y = out[['y']]
fit = out[['fit']]
plotMDS(y, col=ifelse(y$samples$group == "diploid", "red", "blue"))
plotBCV(y)
colnames(y$design)

myContrast <- makeContrasts("karyotypeT21-karyotypediploid", levels = y$design)
qlf <- glmQLFTest(fit, contrast=myContrast)
# qlf <- glmQLFTest(fit, coef=4)

# get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
tt <- topTags(qlf, n = Inf)
tt <- tt$table

plotSmear(qlf, de.tags = rownames(tt)[which(tt$FDR<0.05)])

table(tt$FDR < 0.05)
# tt = tt[!is.na(geneMap$ensID[match(rownames(tt),geneMap$ensID)]),]
# rownames(tt) = geneMap$ensID[match(rownames(tt),geneMap$geneSym)]
tt = annotateGenes(tt,geneMap = geneMap)
tt$DE = (tt$FDR < 0.05 & abs(tt$logFC) > 0.2 )
table(tt$DE)
tt$chr = factor(tt$chr,c(paste0('chr',c(1:22)),'chrM','chrX','chrY',unique(tt$chr[grepl('GL\\d+|KI\\d+',tt$chr)])))
library(ggbeeswarm)
ggplot(tt,aes(chr,logFC))+
  geom_quasirandom(aes(col=DE),size=0.1)+
  geom_boxplot(outlier.shape = NA,alpha=0.5)+
  scale_color_manual(values = c('black','red'))+
  geom_hline(yintercept = c(0,0.5,-0.5))+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank())+
  xlab('') + ggtitle('T21 vs 2n',subtitle = '1 replicate per donor / ~ 0 + karyotype')




##--- OLD age group
obs = obs_fTFC2[obs_fTFC2$age_group == '14_20',]
pb_fTFC2_young = pb_fTFC2[,obs$donor]
dim(pb_fTFC2_young)
# build design matrix     
karyotype <- as.factor(obs[['karyotype']])
print(levels(karyotype))
# age_group <- as.factor(gsub('-','_',obs$age_group))
# print(levels(age_group))
# donor <- as.factor(obs$donor)
# print(levels(donor))

design <- model.matrix(~0+ karyotype)
colnames(design)

## Fit the model
out = edgeR_fitModel(pb=pb_fTFC2_young,obs=obs,group='karyotype',design=design)
y = out[['y']]
fit = out[['fit']]
plotMDS(y, col=ifelse(y$samples$group == "diploid", "red", "blue"))
plotBCV(y)
colnames(y$design)

myContrast <- makeContrasts("karyotypeT21-karyotypediploid", levels = y$design)
qlf <- glmQLFTest(fit, contrast=myContrast)
#qlf <- glmQLFTest(fit, coef=2)

# get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
tt <- topTags(qlf, n = Inf)
tt <- tt$table

plotSmear(qlf, de.tags = rownames(tt)[which(tt$FDR<0.05)])

table(tt$FDR < 0.05)
# tt = tt[!is.na(geneMap$ensID[match(rownames(tt),geneMap$ensID)]),]
# rownames(tt) = geneMap$ensID[match(rownames(tt),geneMap$geneSym)]
tt = annotateGenes(tt,geneMap = geneMap)
tt$DE = (tt$FDR < 0.05 & abs(tt$logFC) > 0.2 )
table(tt$DE)
tt$chr = factor(tt$chr,c(paste0('chr',c(1:22)),'chrM','chrX','chrY',unique(tt$chr[grepl('GL\\d+|KI\\d+',tt$chr)])))
library(ggbeeswarm)
ggplot(tt,aes(chr,logFC))+
  geom_quasirandom(aes(col=DE),size=0.1)+
  geom_boxplot(outlier.shape = NA,alpha=0.5)+
  scale_color_manual(values = c('black','red'))+
  geom_hline(yintercept = c(0,0.5,-0.5))+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank())+
  xlab('') + ggtitle('T21 vs 2n',subtitle = '1 replicate per donor / ~ 0 + karyotype')













## RLE before normalisation
rownames(obs) = obs$X
normCnt = sweep(y$counts,2,STATS = y$samples$norm.factors,FUN = '/')
cpms <- cpm(y, log=FALSE)
ggRLE(dat_x=cpms, annot=obs, col_str='karyotype',isLog=T, isLarge=FALSE,
      ylim = c(-2,2),zero_line=TRUE, zero_col="skyblue", medPoint=FALSE, whisk=1.5)

rownames(cpms) = geneMap$geneSym[match(rownames(cpms),geneMap$ensID)]

g = cpms['SLC5A5',]
g = data.frame(donor = names(g),logCPM = g)
g$donor = gsub('donor_|_0$','',g$donor)
g = merge(g,obs,by='donor')
ggplot(g,aes(karyotype,logCPM))+
  geom_boxplot()+
  facet_wrap(vars(age_group))+
  geom_point()

contr <- makeContrasts(Group60vs40 = Group60um - Group40um, levels=design) #specify comparison between 2 groups
> lrt <- glmLRT(fit, contrast=contr) 








pb_hm11 = Matrix::readMM('~/lustre_mt22/Thyroid/Results/x_2nT21_edgeR_pseudobulk/matrix.mtx.gz')
pb_hm11 = as.matrix(pb_hm11)
obs = read.csv('~/lustre_mt22/Thyroid/Results/x_2nT21_edgeR_pseudobulk/metadata.csv')
colnames(pb_hm11) = obs$X
var = read.delim('~/lustre_mt22/Thyroid/Results/x_2nT21_edgeR_pseudobulk/geneMap.csv',header = F)
rownames(pb_hm11) = var$V1
pb = pb_hm11
group = 'karyotype'
library(edgeR)  
# create an edgeR object with counts and grouping factor
y <- DGEList(pb, group = obs[[group]])

# filter out genes with low counts
print("Dimensions before subsetting:")
print(dim(y))
print("")
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
print("Dimensions after subsetting:")
print(dim(y))
print("")

# normalize
y <- calcNormFactors(y)

PooledCounts <- sumTechReps(y, ID=obs$replicate[match(colnames(y$counts),obs$X)])
y=PooledCounts


# build design matrix      
obs = obs[match(colnames(y$counts),obs$replicate),]
karyotype <- as.factor(obs[['karyotype']])
print(levels(karyotype))
age_group <- as.factor(gsub('-','_',obs$age_group))
print(levels(age_group))
donor <- as.factor(obs$donor)
print(levels(donor))

design <- model.matrix(~ 0 + karyotype + age_group)
# estimate dispersion
y <- estimateDisp(y, design = design)

# fit the model
fit <- glmQLFit(y, design,robust = T)

## Fit the model
out = edgeR_fitModel(pb=pb_fTFC1,obs=obs,group='karyotype',design=design)
y = out[['y']]
fit = out[['fit']]
plotMDS(y, col=ifelse(y$samples$group == "2n", "red", "blue"))
plotBCV(y)
colnames(y$design)

myContrast <- makeContrasts('karyotypeT21-karyotype2n', levels = y$design)
qlf <- glmQLFTest(fit, contrast=myContrast)

# get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
tt <- topTags(qlf, n = Inf)
tt <- tt$table

plotSmear(qlf, de.tags = rownames(tt)[which(tt$FDR<0.05)])

table(tt$FDR < 0.05)
# tt = tt[!is.na(geneMap$ensID[match(rownames(tt),geneMap$ensID)]),]
# rownames(tt) = geneMap$ensID[match(rownames(tt),geneMap$geneSym)]
tt = annotateGenes(tt,geneMap = geneMap)
tt$DE = (tt$FDR < 0.05 & abs(tt$logFC) > 0.5 )
table(tt$DE)
tt$chr = factor(tt$chr,c(paste0('chr',c(1:22)),'chrM','chrX','chrY',unique(tt$chr[grepl('GL\\d+|KI\\d+',tt$chr)])))
library(ggbeeswarm)
ggplot(tt,aes(chr,logFC))+
  geom_quasirandom(aes(col=DE),size=0.1)+
  geom_boxplot(outlier.shape = NA,alpha=0.5)+
  scale_color_manual(values = c('black','red'))+
  geom_hline(yintercept = c(0,0.5,-0.5))+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
        panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank())+
  xlab('') + ggtitle('T21 vs 2n',subtitle = '1 replicate per donor / ~ 0 + karyotype')








### Perform MILO to assess differential abundance on fetal Thyroid Immune cells ####

## Install development version
#devtools::install_github("MarioniLab/miloR") 

library(miloR)
library(SingleCellExperiment)
library(tidyverse)
library(scater)
library(dplyr)
library(patchwork)
library(Seurat)
library(RColorBrewer)
library(ggbeeswarm)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")



import_seurat_object = function(object_type,data_dir){
  
  srat_fp = file.path(data_dir,paste0(object_type,'_sratObject.RDS'))
  if(!file.exists(srat_fp)){
    ## Import seurat object
    srat = Read10X(data_dir)
    srat = CreateSeuratObject(srat)
    srat = standard_clustering(srat)
    ## Add metadata
    mdat = read.csv(file.path(data_dir,'metadata.csv'))
    mdat$cellID = mdat$X
    srat@meta.data = cbind(srat@meta.data,mdat[match(rownames(srat@meta.data),mdat$cellID),])
    
    # Add scVI latent space
    scVI = read.csv(file.path(data_dir,'obsm_scVI40.csv'))
    scVI = column_to_rownames(scVI,var = 'X')
    table(rownames(scVI) %in% srat$cellID)
    srat@reductions$scVI = srat@reductions$pca
    srat@reductions$scVI@cell.embeddings = as.matrix(scVI)
    # Add umap latent space
    umap = read.csv(file.path(data_dir,'obsm_umap.csv'))
    umap = column_to_rownames(umap,var = 'X')
    colnames(umap) = c('UMAP_1','UMAP_2')
    srat@reductions$umap@cell.embeddings = as.matrix(umap)
    
    saveRDS(srat,srat_fp)
  }else{
    srat = readRDS(srat_fp)
  }
  
  return(srat)  
}




#'SCVI'
differential_abundance_testing = function(srat,reduced_dims = 'PCA',k = 50, d = 30,prop=0.1,da_design,sample_col='donor',annot_col="celltype",design_formula='',outDir,...){
  ## Convert seurat object to SCE using function from scater package
  srat.sce = as.SingleCellExperiment(srat)
  ## Visualize data
  #plotReducedDim(srat.sce, colour_by="karyotype", dimred = reduced_dims)
  
  # Create a Milo object
  srat.milo <- Milo(srat.sce) # This extends the SingleCellExperiment class to store information about neighbourhoods on the KNN graph.
  srat.milo
  
  # Construct KNN graph
  srat.milo <- buildGraph(srat.milo, k = k, d = d, reduced.dim = reduced_dims)
  
  # Defining representative neighbourhoods on the KNN graph
  srat.milo <- makeNhoods(srat.milo, prop = prop, k = k, d=d, refined = TRUE, reduced_dims = reduced_dims)
  plotNhoodSizeHist(srat.milo)
  
  # Counting cells in neighbourhoods
  srat.milo <- countCells(srat.milo, meta.data = data.frame(colData(srat.milo)), sample=sample_col)
  head(nhoodCounts(srat.milo))
  
  # Computing neighbourhood connectivity
  srat.milo <- calcNhoodDistance(srat.milo, d=d, reduced.dim = reduced_dims)
  
  # Testing
  da_results <- testNhoods(srat.milo, design = as.formula(design_formula), design.df = da_design,reduced.dim = reduced_dims,...)
  
  # Annotate the neighbourhoods
  da_results <- annotateNhoods(srat.milo, da_results, coldata_col = annot_col)
  
  if(!is.null(outDir)){
    if(dir.exists(outDir)){
      write.csv(da_results,file.path(outDir,'DA_results.csv'))    
    }
  }
  
  return(list(da_results=da_results,srat.milo=srat.milo))
}






###---------------------------------------------------------------------------------------------##
###      Differential abundance in 2n + T21 fetal thy_immune compartment, by karyotype       #####
###--------------------------------------------------------------------------------------------##

object_type = 'thy_imm_2n_T21' # or 'thy_imm_2n'
data_dir = '~/lustre_mt22/Thyroid/Data/adata_imm_all_T21_2n_harm_anno'
design_formula = '~ 0 + karyotype + age_group'


outDir = file.path('/lustre/scratch126/cellgen/team292/mt22/ThyroidDev/milo_analysis',object_type)
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}



## 1. Import thyroid seurat object  ##
srat = import_seurat_object(object_type,data_dir)

## Add age group
srat$age_group = ''
srat$age_group[srat$pcw %in% c(9,10)] = '9.10'
srat$age_group[srat$pcw %in% c(11,12,13)] = '11.13'
srat$age_group[srat$pcw >= 14 ] = '14.20'

## Add immune-lineage
srat$lineage = ifelse(srat$celltype %in% c('imm_PreB_cells','imm_B_cells'),'B cells',
                      ifelse(srat$celltype %in% c('imm_ILC','imm_Cycling_T_cells','imm_DN/DP_T_cells','imm_T_cells','imm_Type1_Innate_T','imm_Cycling_T',
                                                  'imm_Cycling_NK_cells','imm_NK_cells'),'NK/T cells',
                             ifelse(srat$celltype %in% c('imm_DC_prec','imm_DC2','imm_PDC/DC1','imm_Macrophages','imm_Monocytes'),'Myeloid cells',
                                    ifelse(srat$celltype %in% c('imm_MEMP','imm_Mast_cells'),'MK/Mast cells','others'))))

#"celltype","karyotype", "leiden","phase"
DimPlot(srat,group.by = 'karyotype',label = T,label.box = T,repel = T,cols = col25) + theme(legend.position = 'none')
DimPlot(srat,group.by = 'celltype',label = T,label.box = T,repel = T,cols = col25) + theme(legend.position = 'none')
DimPlot(srat,group.by = 'gender',label = T,label.box = T,repel = T,cols = col25) + theme(legend.position = 'none')
DimPlot(srat,group.by = 'age_group',label = T,label.box = T,repel = T,cols = col25) + theme(legend.position = 'none')
DimPlot(srat,group.by = 'lineage',label = T,label.box = T,repel = T,cols = col25) + theme(legend.position = 'none')

## Plot fraction of cells per group
dd = as.data.frame(table(srat$donor,srat$karyotype,srat$celltype,srat$gender,srat$age_group,srat$lineage))
colnames(dd) = c('donor','karyotype','celltype','gender','age_group','lineage','nCell')
dd = dd %>% group_by(donor) %>% mutate(totalCell = sum(nCell))
dd$frac = dd$nCell/dd$totalCell
dd = dd[dd$nCell > 0,]
dd$age_group = factor(dd$age_group,c('9.10','11.13','14.20'))
dd$lineage = factor(dd$lineage,c('B cells','NK/T cells', 'Myeloid cells','MK/Mast cells'))
dd$celltype = gsub('^imm_','',dd$celltype)
dd$celltype = factor(dd$celltype,c('PreB_cells','B_cells',
                                   'ILC','T_cells','Cycling_T','Type1_Innate_T','NK_cells','Cycling_NK_cells',
                                   'DC2','PDC/DC1','Monocytes','Macrophages','Mast_cells'))
p = ggplot(dd,aes(karyotype,frac))+
  geom_boxplot(aes(fill=age_group),outlier.shape = NA)+
  geom_quasirandom(size=0.5,width = 0.5)+
  scale_fill_manual(values = c('9.10'='lightblue','11.13'='gold','14.20'='red'),name='Age group')+
  facet_grid(.~lineage+celltype)+
  theme_classic()+xlab('Karyotype') + ylab('Fraction of all immune cells (per donor)')+
  theme(panel.border = element_rect(fill=NA,colour = 'black',linewidth = 1),
        axis.line = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))



pdf(file.path(outDir,'imm_celltype_fraction_perDonor.pdf'),width = 12,height = 4)
print(p)
dev.off()




# Defining experimental design

da_design = distinct(srat@meta.data[,c("donor", "karyotype","gender",'age_group')])
rownames(da_design) <- paste0(da_design$donor)
head(da_design)


## Differential abundance testing

milo_results = differential_abundance_testing(srat,reduced_dims = 'SCVI',k = 30, d = 30,prop=0.1,
                                              da_design=da_design,sample_col='donor',annot_col="celltype",design_formula=design_formula,outDir = outDir)
  
srat.milo = milo_results[['srat.milo']]
da_results = milo_results[['da_results']]
da_results %>%
  arrange(SpatialFDR) %>%
  head()

plotNhoodMA(da_results)


## Inspecting DA testing results
## Look at distribution of uncorrected P-value to help verify that the test was balanced.
ggplot(da_results, aes(PValue)) + geom_histogram(bins=100)+
  theme_classic()

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = 1)+ ## Mark significance threshold (10% FDR)
  theme_classic()


## Visualize DA results relating them to the sc embedding
srat.milo <- buildNhoodGraph(srat.milo)
# Plot single-cell UMAP
umap_pl <- plotReducedDim(srat.milo, dimred = "UMAP", colour_by="celltype", text_by = "celltype", text_size = 3) +
  guides(fill="none")
# Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(srat.milo, da_results, layout="UMAP",alpha=0.05)

p1 = umap_pl + nh_graph_pl +
  plot_layout(guides="collect")

pdf(file.path(outDir,'milo_umap.pdf'),width = 14,height = 8)
print(p1)
dev.off()


## Visualizing whether DA is particularly evident in certain cell types
ggplot(da_results, aes(celltype_fraction)) + geom_histogram(bins=100)+
  theme_classic()+
  geom_vline(xintercept = 0.7,col='red')
da_results$celltype <- ifelse(da_results$celltype_fraction < 0.7, "Mixed", da_results$celltype)
# visualize the distribution of DA Fold Changes in different cell types
p = plotDAbeeswarm(da_results, group.by = "celltype")
pdf(file.path(outDir,'milo_DAbeeswarm_byCellType.pdf'),width = 9,height = 7)
print(p)
dev.off()





## Identifying signatures of DA subpopulations

# ## Add log normalized count to Milo object
# srat.milo <- logNormCounts(srat.milo)
# 
# da_results$NhoodGroup <- as.numeric(da_results$SpatialFDR < 0.1 & da_results$logFC < 0)
# da_nhood_markers <- findNhoodGroupMarkers(srat.milo, da_results, subset.row = rownames(srat.milo),
#                                           aggregate.samples = TRUE, sample_col = sample_col)
# markers <- dge_smp[which(dge_smp$adj.P.Val_1 < 0.1 ), "GeneID"]
# logcounts(srat.milo) <- log1p(counts(srat.milo))
# srat.milo <- calcNhoodExpression(srat.milo, subset.row=markers)


##----- Automatic grouping of neighbourhoods -----##
## Run buildNhoodGraph to store nhood adjacency matrix
srat.milo <- buildNhoodGraph(srat.milo)

## Find groups
da_results <- groupNhoods(srat.milo, da_results, max.lfc.delta = 10)
head(da_results)
plotNhoodGroups(srat.milo, da_results, layout="UMAP") 
plotDAbeeswarm(da_results, "NhoodGroup")

# check how changing the grouping parameters changes the groups we obtain, starting with the LFC delta by plotting with different values of max.lfc.delta 
plotDAbeeswarm(groupNhoods(srat.milo, da_results, max.lfc.delta = 1) , group.by = "NhoodGroup") + ggtitle("max LFC delta=1")
plotDAbeeswarm(groupNhoods(srat.milo, da_results, max.lfc.delta = 2)   , group.by = "NhoodGroup") + ggtitle("max LFC delta=2")
plotDAbeeswarm(groupNhoods(srat.milo, da_results, max.lfc.delta = 3)   , group.by = "NhoodGroup") + ggtitle("max LFC delta=3")

set.seed(42)
da_results <- groupNhoods(srat.milo, da_results, max.lfc.delta = 10, overlap=1)
plotNhoodGroups(srat.milo, da_results, layout="UMAP")

# Finding gene signatures for neighbourhoods
library(scran)
## Exclude zero counts genes
keep.rows <- rowSums(logcounts(srat.milo)) != 0
srat.milo <- srat.milo[keep.rows, ]

## Find HVGs
set.seed(101)
dec <- modelGeneVar(srat.milo)
hvgs <- getTopHVGs(dec, n=2000)

# this vignette randomly fails to identify HVGs for some reason
if(!length(hvgs)){
  set.seed(42)
  dec <- modelGeneVar(srat.milo)
  hvgs <- getTopHVGs(dec, n=2000)
}

head(hvgs)
set.seed(42)
nhood_markers <- findNhoodGroupMarkers(srat.milo, da_results, subset.row = hvgs, 
                                       aggregate.samples = TRUE, sample_col = "donor")

head(nhood_markers)

# check out the markers for group 6
gr6_markers <- nhood_markers[c("logFC_6", "adj.P.Val_6")] 
colnames(gr6_markers) <- c("logFC", "adj.P.Val")

head(gr6_markers[order(gr6_markers$adj.P.Val), ])


## Visualize detected markers
ggplot(nhood_markers, aes(logFC_6, -log10(adj.P.Val_6 ))) + 
  geom_point(alpha=0.5, size=0.5) +
  geom_hline(yintercept = 3)
markers <- nhood_markers$GeneID[nhood_markers$adj.P.Val_6 < 0.01 & nhood_markers$logFC_6 > 0]

set.seed(42)
plotNhoodExpressionGroups(srat.milo, da_results, features=intersect(rownames(srat.milo), markers[1:10]),
                          #subset.nhoods = da_results$NhoodGroup %in% c('6','5'), 
                          scale=TRUE,
                          grid.space = "fixed")


plotNhoodExpressionDA(srat.milo, da_results, features = markers,
                      #subset.nhoods = da_results$celltype %in% c("Anterior Primitive Streak", "Def. endoderm", "Gut", "Visceral endoderm"),
                      assay="logcounts",
                      scale_to_1 = TRUE, cluster_features = TRUE
)










###---------------------------------------------------------------------------------------------##
###       Differential abundance in Diploid fetal thy_immune compartment, by age group       #####
###--------------------------------------------------------------------------------------------##

# SUMMARY: Due to large imbalances in biological replicates (number of samples per age group), its difficult to estimate within group variance 
#          --> Milo is NOT appropriate for this particular datasets
#          I tried downsampling per donor, which reduced the imbalances in cell number. However, this still didn't resolve the issue with group imbalances
# Found another method called scCODA, but depends on how much emphasis we want to put on this analysis - might be ok with just the boxplot

object_type = 'thy_imm_2n'
data_dir = '~/lustre_mt22/Thyroid/Data/adata_imm_all_2n_harm_anno'
design_formula = '~ 0 + age_group + donor'


outDir = file.path('/lustre/scratch126/cellgen/team292/mt22/ThyroidDev/milo_analysis',object_type)
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

## 1. Import thyroid seurat object  ##
srat = import_seurat_object(object_type,data_dir)

## Add age group
srat$age_group = ''
srat$age_group[srat$pcw %in% c(9,10)] = '9.10'
srat$age_group[srat$pcw %in% c(11,12,13)] = '11.13'
srat$age_group[srat$pcw >= 14 ] = '14.20'

## Add immune-lineage
srat$lineage = ifelse(srat$celltype %in% c('imm_Pro-B_cells','imm_B_cells'),'B cells',
                      ifelse(srat$celltype %in% c('imm_ILCs','imm_Cycling_T_cells','imm_DN/DP_T_cells','imm_T_cells','imm_Cycling_NK','imm_NK_cells'),'NK/T cells',
                             ifelse(srat$celltype %in% c('imm_DC_prec','imm_DC2','imm_Macrophages','imm_Monocytes'),'Myeloid cells',
                                    ifelse(srat$celltype %in% c('imm_MEMP','imm_Mast_cells'),'MK/Mast cells','others'))))

#"celltype","karyotype", "leiden","phase"
DimPlot(srat,group.by = 'karyotype',label = T,label.box = T,repel = T,cols = col25) + theme(legend.position = 'none')
DimPlot(srat,group.by = 'celltype',label = T,label.box = T,repel = T,cols = col25) + theme(legend.position = 'none')
DimPlot(srat,group.by = 'gender',label = T,label.box = T,repel = T,cols = col25) + theme(legend.position = 'none')
DimPlot(srat,group.by = 'age_group',label = T,label.box = T,repel = T,cols = col25) + theme(legend.position = 'none')
DimPlot(srat,group.by = 'lineage',label = T,label.box = T,repel = T,cols = col25) + theme(legend.position = 'none')

## Plot fraction of cells per group
dd = as.data.frame(table(srat$donor,srat$karyotype,srat$celltype,srat$gender,srat$age_group,srat$lineage))
colnames(dd) = c('donor','karyotype','celltype','gender','age_group','lineage','nCell')
dd = dd %>% group_by(donor) %>% mutate(totalCell = sum(nCell))
dd$frac = dd$nCell/dd$totalCell
dd = dd[dd$nCell > 0,]
dd$age_group = factor(dd$age_group,c('9.10','11.13','14.20'))
dd$lineage = factor(dd$lineage,c('B cells','NK/T cells', 'Myeloid cells','MK/Mast cells'))
dd$celltype = gsub('^imm_','',dd$celltype)
if(object_type = "thy_imm_2n"){
  p = ggplot(dd,aes(age_group,frac))+
    geom_boxplot(aes(fill=age_group),outlier.shape = NA)+
    geom_quasirandom(size=0.5,width = 0.5)+
    scale_fill_manual(values = c('9.10'='lightblue','11.13'='gold','14.20'='red'),name='Age group')+
    facet_grid(.~lineage+celltype)+
    theme_classic()+xlab('Age group (pcw)') + ylab('Fraction of all immune cells (per donor)')+
    theme(panel.border = element_rect(fill=NA,colour = 'black',linewidth = 1),
          axis.line = element_blank(),
          strip.background = element_blank(),
          axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))
  
}

pdf(file.path(outDir,'imm_celltype_fraction_perDonor.pdf'),width = 12,height = 4)
print(p)
dev.off()


##----- DOWN-SAMPLING ------####
# I first did the analysis without down-sampling --> the results always show signficant enrichment in the older age group, potentially as there are a lot more cells / samples from older age group.
# Checking MA plot of the results, not centralised around 0 on the y-axis --> cell number imblance is too strong confounding factor
# Milo does normalised for cell number but this is not enough - see discussion here: https://github.com/MarioniLab/miloR/issues/208
# --> Due to huge imblanaces in number of cells / samples per age group, authors advised to down-sample the object instead

# I think downsampling should be done at the donor level, such that each donor have the same total number of cells but randomly sampling across cell types
# --> to ensure that the relative proportion of different cell types are preserved

# There's 1 sample from pcw=10 with 64 cells... maybe we just look at the 11-13 vs 14-20 age group?
minCell_perDonor = table(srat$donor)
minCell_perDonor = minCell_perDonor[minCell_perDonor == min(minCell_perDonor)]
minCell_perDonor = 58

srat$toKeep = F

for(d in unique(srat$donor)){
  if(sum(srat$donor == d) > minCell_perDonor){
    set.seed(1245)
    srat$toKeep[srat$donor == d] = (srat$cellID[srat$donor==d] %in% sample(srat$cellID[srat$donor == d],minCell_perDonor))
  }else if(sum(srat$donor == d) < minCell_perDonor){
    srat$toKeep[srat$donor == d] = F
  }else{
    srat$toKeep[srat$donor == d] = T
  }
}

table(srat$toKeep,srat$donor)
srat = subset(srat,subset = toKeep == TRUE)
srat = standard_clustering(srat)

# Defining experimental design

da_design = distinct(srat@meta.data[,c("donor","gender",'age_group')])  
rownames(da_design) <- paste0(da_design$donor)


## Differential abundance testing
contrast.all <- c("9.10 - 11.13", "11.13 - 14.20", "9.10 - 14.20")
result = list()
for(contrast in contrast.all){
  
  result[[contrast]] = differential_abundance_testing(srat,reduced_dims = 'PCA',k = 30, d = 30,prop=0.1,
                                                      da_design=da_design,sample_col='donor',annot_col="celltype",design_formula=design_formula,outDir = NULL,
                                                      fdr.weighting="graph-overlap", model.contrasts = contrast)
}

# # this is the edgeR code called by `testNhoods`
# model <- model.matrix(~ 0 + age_group, data=da_design)
# mod.constrast <- makeContrasts(contrasts=contrast.all, levels=model)
# contrast1.res <- testNhoods(srat.milo, design= ~ 0+ age_group, design.df=da_design, fdr.weighting="graph-overlap", model.contrasts = contrast.all)
# head(contrast1.res)

## Group results into 1 table
combined_results = data.frame()
for(i in 1:length(result)){
  result[[i]][[1]]$contrast = names(result)[i]
  combined_results = rbind(combined_results,result[[i]][[1]])
}

write.csv(combined_results,file.path(outDir,'DA_allContrasts_results.csv'))


head(da_design)



# ## Differential abundance testing

srat.milo = result[[2]][['srat.milo']]
da_results = result[[2]][['da_results']]
da_results %>%
  arrange(SpatialFDR) %>%
  head()

plotNhoodMA(da_results)

nh_counts <- countCells(srat.milo, meta.data = data.frame(colData(srat.milo)), samples = "donor")
boxplot(cell_count ~ age_group, data = nh_counts)

## Inspecting DA testing results
## Look at distribution of uncorrected P-value to help verify that the test was balanced.
ggplot(da_results, aes(PValue)) + geom_histogram(bins=100)+
  theme_classic()

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = 1)+ ## Mark significance threshold (10% FDR)
  theme_classic()


## Visualize DA results relating them to the sc embedding
srat.milo <- buildNhoodGraph(srat.milo)
# Plot single-cell UMAP
umap_pl <- plotReducedDim(srat.milo, dimred = "UMAP", colour_by="celltype", text_by = "celltype", text_size = 3) +
  guides(fill="none")
# Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(srat.milo, da_results, layout="UMAP",alpha=0.05)

p1 = umap_pl + nh_graph_pl +
  plot_layout(guides="collect")

pdf(file.path(outDir,'milo_umap.pdf'),width = 14,height = 8)
print(p1)
dev.off()


## Visualizing whether DA is particularly evident in certain cell types
ggplot(da_results, aes(celltype_fraction)) + geom_histogram(bins=100)+
  theme_classic()+
  geom_vline(xintercept = 0.7,col='red')
da_results$celltype <- ifelse(da_results$celltype_fraction < 0.7, "Mixed", da_results$celltype)
# visualize the distribution of DA Fold Changes in different cell types
p = plotDAbeeswarm(da_results, group.by = "celltype")
pdf(file.path(outDir,'milo_DAbeeswarm_byCellType.pdf'),width = 9,height = 7)
print(p)
dev.off()






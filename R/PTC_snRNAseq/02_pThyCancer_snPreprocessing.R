##--- Preprocess the 10X single-nuclei paediatric thyroid cancer dataset ---##


#outDir = "/nfs/team274/mt22/Thyroid/Results_2505/PTC_scRNAseq/02_pThyCancer_snPreprocessing"
outDir = "Results/2505/PTC_snRNAseq/02_pThyCancer_snPreprocessing"
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

#setwd(outDir)


##----------------##
##   Libraries  ####
##----------------##
library(Seurat)
library(tidyverse)
#source("R/utils/misc.R")
source("R/utils/sc_basicQC.R")


##----------------------------##
##   Set Global parameters  ####
##----------------------------##
maxMT = 30
minGenes = 300
minUMIs = 500  
maxBadFrac = 0.5
numPCs = 75
clusteringRes = 10
skipScrub = F
skipSoup = F
scrubScoreMax = 0.5
scrubPath='./scrubletScores.tsv'
scPath="./strainedCounts_highRho"
doPlot=T
verbose = T
skipIfExists=F
keepMTCells=T
rho_max_limit = NULL


##----------------------------------##
##   Preprocessing snRNAseq data  ####
##----------------------------------##
### 1. Import cellranger output data
###    Run SoupX
###    Subset to keep only cells present in the original publications
###    Add cell labels (as published)

outDir_sub = file.path(outDir,'may25') 
#cleanCountDir = file.path(outDir_sub,'cleanCount')

if(!is.null(rho_max_limit)){
  plotDir = paste0(file.path(outDir_sub,'scCancerThyroid_10X_rhoLim'),rho_max_limit,'_')
  outPath = paste0(file.path(outDir_sub,'scCancerThyroid_10X_rhoLim'),rho_max_limit,'_')
}else{
  plotDir = paste0(file.path(outDir_sub,'scCancerThyroid_10X_rhoLimNone_'))
  outPath = paste0(file.path(outDir_sub,'scCancerThyroid_10X_rhoLimNone_'))
}

cleanSrat_fp = ifelse(keepMTCells,paste0(outPath,'_clean_withMTCells.RDS'),paste0(outPath,'_clean_noMTCells.RDS'))
if(file.exists(cleanSrat_fp) & skipIfExists){
  cleanSrat = readRDS(cleanSrat_fp)  
}else{
  if(!dir.exists(outDir_sub)){
    message(sprintf('Creating output directory'))
    dir.create(outDir_sub,recursive = T)
  }
  
  #setwd(outDir_sub)
    
  dataDirs = c('Data/thyroid_10X/cellranger710_count_46320_SB_Thy_R13236839_GRCh38-2020-A/filtered_feature_bc_matrix/',
               'Data/thyroid_10X/cellranger710_count_46320_SB_Thy_R13236840_GRCh38-2020-A/filtered_feature_bc_matrix/',
               'Data/thyroid_10X/cellranger710_count_46320_SB_Thy_R13236841_GRCh38-2020-A/filtered_feature_bc_matrix/',
               'Data/thyroid_10X/cellranger710_count_46320_SB_Thy_R13236842_GRCh38-2020-A/filtered_feature_bc_matrix/',
               'Data/thyroid_10X/cellranger710_count_46320_SB_Thy_R13236843_GRCh38-2020-A/filtered_feature_bc_matrix/',
               'Data/thyroid_10X/cellranger710_count_46320_SB_Thy_R13236844_GRCh38-2020-A/filtered_feature_bc_matrix/',
               # Y46
               'Data/thyroid_10X/cellranger710_count_48713_CG_SB_NB14695466_GRCh38-2020-A/filtered_feature_bc_matrix/',
               'Data/thyroid_10X/cellranger710_count_48713_CG_SB_NB14695467_GRCh38-2020-A/filtered_feature_bc_matrix/',
               'Data/thyroid_10X/cellranger710_count_48853_CG_SB_NB14664105_GRCh38-2020-A/filtered_feature_bc_matrix/',
               'Data/thyroid_10X/cellranger710_count_48853_CG_SB_NB14664106_GRCh38-2020-A/filtered_feature_bc_matrix/',
               'Data/thyroid_10X/cellranger710_count_48853_CG_SB_NB14664107_GRCh38-2020-A/filtered_feature_bc_matrix/',
               'Data/thyroid_10X/cellranger710_count_48853_CG_SB_NB14664108_GRCh38-2020-A/filtered_feature_bc_matrix/',
               'Data/thyroid_10X/cellranger710_count_48853_CG_SB_NB14664109_GRCh38-2020-A/filtered_feature_bc_matrix/',
               'Data/thyroid_10X/cellranger710_count_48853_CG_SB_NB14664110_GRCh38-2020-A/filtered_feature_bc_matrix/'
               )
  

  cleanCountDir = dataDirs
  names(dataDirs) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_SB_','',dataDirs))
  names(dataDirs) = gsub('_','.',names(dataDirs))
  
  dataDirs=dataDirs[file.exists(dataDirs)]
  dataDirs=dataDirs[sapply(dataDirs, function(x){length(list.files(x))>0})]
  print(n_distinct(dataDirs))
  
  
  metadata = NULL
  matchBy = NULL
  is10X=TRUE
  
  
  # Run basicQC
  message('\nPerforming scRNAseq QC...')
  QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,
                      clusteringRes=clusteringRes,cleanCountDir=cleanCountDir,
                      skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
                      metadata=metadata,matchBy=matchBy,scPath=scPath,rho_max_limit=rho_max_limit,
                      outPath=outPath,skipIfExists=skipIfExists,
                      doPlot=doPlot,plotDir=plotDir,verbose=verbose,is10X=is10X)
  
  cleanSrat = QC.output[[1]]
  df.out = QC.output[[2]]
  
  write.csv(df.out,paste0(outPath,'qc_summary.csv'))
}
##--- Genotyping the WGS + 10X single-nuclei paediatric PTC data ---##

outDir = "Results/2505/PTC_snRNAseq/01_pThyCancer_genotyping"
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

setwd(outDir)

##----------------##
##   Libraries  ####
##----------------##
library(alleleIntegrator)

##----------------------------##
##   Set Global parameters  ####
##----------------------------##
refGenome = '/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa'
refGenome10X = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/fasta/genome.fa'
liftChain = 'hg19ToHg38_noChr.over.chain'
nParallel=48
skipIfExists=T

##----------------------------------------##
##   Get list of DNA and RNA bam files  ####
##----------------------------------------##

#----- RNA BAMs
bams10X = c('Data/thyroid_10X/cellranger710_count_46320_SB_Thy_R13236839_GRCh38-2020-A/possorted_genome_bam.bam',
            'Data/thyroid_10X/cellranger710_count_46320_SB_Thy_R13236840_GRCh38-2020-A/possorted_genome_bam.bam',
            'Data/thyroid_10X/cellranger710_count_46320_SB_Thy_R13236841_GRCh38-2020-A/possorted_genome_bam.bam',
            'Data/thyroid_10X/cellranger710_count_46320_SB_Thy_R13236842_GRCh38-2020-A/possorted_genome_bam.bam',
            'Data/thyroid_10X/cellranger710_count_46320_SB_Thy_R13236843_GRCh38-2020-A/possorted_genome_bam.bam',
            'Data/thyroid_10X/cellranger710_count_46320_SB_Thy_R13236844_GRCh38-2020-A/possorted_genome_bam.bam',
            # Y46
            'Data/thyroid_10X/cellranger710_count_48713_CG_SB_NB14695466_GRCh38-2020-A/possorted_genome_bam.bam',
            'Data/thyroid_10X/cellranger710_count_48713_CG_SB_NB14695467_GRCh38-2020-A/possorted_genome_bam.bam',
            'Data/thyroid_10X/cellranger710_count_48853_CG_SB_NB14664105_GRCh38-2020-A/possorted_genome_bam.bam',
            'Data/thyroid_10X/cellranger710_count_48853_CG_SB_NB14664106_GRCh38-2020-A/possorted_genome_bam.bam',
            'Data/thyroid_10X/cellranger710_count_48853_CG_SB_NB14664107_GRCh38-2020-A/possorted_genome_bam.bam',
            'Data/thyroid_10X/cellranger710_count_48853_CG_SB_NB14664108_GRCh38-2020-A/possorted_genome_bam.bam',
            'Data/thyroid_10X/cellranger710_count_48853_CG_SB_NB14664109_GRCh38-2020-A/possorted_genome_bam.bam',
            'Data/thyroid_10X/cellranger710_count_48853_CG_SB_NB14664110_GRCh38-2020-A/possorted_genome_bam.bam')

names(bams10X) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_SB_','',bams10X))

if(length(bams10X) > 0){
  bams10X = bams10X[file.exists(bams10X)]  
}
# Check that each bams10X file has a unique name
if(length(unique(names(bams10X))) != length(bams10X)){
  stop(sprintf('Duplicated bams10X names detected: %s',names(bams10X)[duplicated(names(bams10X))]))
}


#----- bulk RNA BAMs
rnaBAMs = list.files('Thyroid/Data/inhouse_bulkRNA_thyroid/tic-3426/results',recursive = T,full.names = T,
                     pattern = 'Aligned.sortedByCoord.out.bam$')
names(rnaBAMs) = basename(dirname(rnaBAMs))
rnaBAMs = rnaBAMs[file.exists(rnaBAMs)]


#----- DNA BAMs
dnaBAMs = c('Data/thyroid_DNA/PD62328b/PD62328b.sample.dupmarked.bam',
            'Data/thyroid_DNA/PD62328a/PD62328a.sample.dupmarked.bam',
            'Data/thyroid_DNA/PD64651a/PD64651a.sample.dupmarked.bam',
            'Data/thyroid_DNA/PD64651d/PD64651d.sample.dupmarked.bam',
            'Data/thyroid_DNA/PD64654b/PD64654b.sample.dupmarked.bam')
names(dnaBAMs) = gsub('\\..*$','',basename(dnaBAMs))

if(length(dnaBAMs) > 0){
  dnaBAMs = dnaBAMs[file.exists(dnaBAMs)]
}

# Check that each dnaBAM file has a unique name
if(length(unique(names(dnaBAMs))) != length(dnaBAMs)){
  stop(sprintf('Duplicated bams10X names detected: %s',names(dnaBAMs)[duplicated(names(dnaBAMs))]))
}




##-------------------------------------##
##    Check genotype consistency     ####
##-------------------------------------##

#Are all the BAMs you're going to use from the same individual?  Check before you start

genoCheck = matchBAMs(BAMs = c(dnaBAMs,rnaBAMs,bams10X),
                      refGenomes = rep(c(refGenome,refGenome10X,refGenome10X),c(length(dnaBAMs),length(rnaBAMs),length(bams10X))),
                      outputs = file.path(outDir,paste0(c(names(dnaBAMs),names(rnaBAMs),names(bams10X)),'_genotypeCheck.tsv')),
                      liftOvers=rep(c(liftChain,liftChain,liftChain),c(length(dnaBAMs),length(rnaBAMs),length(bams10X))),
                      is10X=rep(c(FALSE,FALSE,TRUE),c(length(dnaBAMs),length(rnaBAMs),length(bams10X))),
                      nParallel=nParallel,nMaxSim=2,nChunks=24,skipIfExists=skipIfExists)
#If anything is less than 0.8 and you should be concerned...
message(sprintf("The minimum similarity found was %g",min(genoCheck$ibs$ibs)))

plotDir = outDir
if(!is.null(plotDir)){
  library(RColorBrewer)
  library(circlize)
  library(ComplexHeatmap)
  bamGrouping = NULL
  colPal = 'Greens'
  colFun = suppressWarnings(brewer.pal(100, colPal))
  colFun = colorRamp2(seq(0.5, 1, length.out = length(colFun)),colFun)
  hm = Heatmap(genoCheck$ibs$ibs, col = colFun, name = "IBS", show_row_names = TRUE,
               show_column_names = TRUE, show_row_dend = FALSE,
               show_column_dend = FALSE, row_title_rot = 0, column_split = bamGrouping,
               row_split = bamGrouping)
  pdf(paste0(plotDir,'/genotypeCheck.pdf'),width = 20,height = 17)
  draw(hm)
  dev.off()
  
}


## Generate figures for paediatric PTC single-nuclei RNA-seq dataset
setwd('~/FetalThyroidAtlas/')
plotDir = 'Figures/2505'
if(!dir.exists(plotDir)){
  dir.create(plotDir,recursive = T)
}

##----------------##
##   Libraries  ####
##----------------##
library(Seurat)
library(tidyverse)
library(plotrix)



##-------------------------------##
##       Generating plots      ####
##-------------------------------##
source('R/utils/misc.R')


## Import paediatric PTC snRNA-seq dataset
srat_fp = 'Results/2505/PTC_snRNAseq/03_pThyCancer_annotation/pPTC_clean_soupedXrhoLimNone_annotated_2505.RDS'
mdat_fp = 'Results/2505/PTC_snRNAseq/03_pThyCancer_annotation/pPTC_clean_soupedXrhoLimNone_annotated_2505_mdat.csv'

srat = readRDS('Results/2505/PTC_snRNAseq/03_pThyCancer_annotation/pPTC_clean_soupedXrhoLimNone_annotated_2505.RDS')

# SuppFigure A: UMAP of snRNAseq dataset
suppfig_a_pPTC_snRNAseq_UMAP = function(){
  
  if(file.exists(file.path(plotDir,'suppfig_a_pPTC_snRNAseq_UMAP_rawData.tsv'))){
    dd = read.delim(file.path(plotDir,'suppfig_a_pPTC_snRNAseq_UMAP_rawData.tsv'),header = T,sep = '\t')
    
  }else{
    mdat = read.csv(mdat_fp,row.names = 1)
    
    dd = mdat[!mdat$annot %in% c('doublets','unknown','lowQual'),]
    dd$annot[dd$annot == 'Tumour' & dd$etiology == 'Met'] = 'Tumour cells (Met.)'
    dd$annot[dd$annot == 'Tumour' & dd$etiology != 'Met'] = 'Tumour cells'
    dd$annot = gsub('_',' ',dd$annot)
    
    set.seed(2397)
    
    dd=dd[sample(1:nrow(dd),100000),]
    ggplot(dd,aes(UMAP_1,UMAP_2,col=annot))+
      geom_point(size=0.01)+
      scale_color_manual(values = col25)+
      theme_classic()
  }
  
  
  plotFun_celltype = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0,0,0.8,0))
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','Paediatric PTC'),
         frame.plot=F)
    
    if(!noPlot){
      celltype_cols = col25[1:n_distinct(dd$annot)]
      celltype_cols = colAlpha(celltype_cols,alphas = 0.7)
      names(celltype_cols) = unique(dd$annot)
      celltype_cols[names(celltype_cols) == 'Tumour cells'] = colAlpha(grey(0.4),0.7)
      celltype_cols[names(celltype_cols) == 'Tumour cells (Met.)'] = colAlpha(grey(0.7),0.7)
      
      points(dd$UMAP_1,dd$UMAP_2,
             col = celltype_cols[as.character(dd$annot)],
             pch = 19,
             cex=0.07)
      
    }
    
    if(!noFrame){
      #Add coloured labels
      
      mids = aggregate(cbind(UMAP_1,UMAP_2) ~ annot,data=dd,FUN=mean)
      mids$label = mids$annot
      # mids$legend = factor(paste0(mids$label,' - ',mids$broadLineage),levels = c('1 - HSC & prog.',
      #                                                                            '2 - Megakaryocytes',
      #                                                                            '3 - Mast.cell',
      #                                                                            '4 - Erythroblasts',
      #                                                                            '5 - Monocyte/Macrophage',
      #                                                                            '6 - Dendritic cells',
      #                                                                            '7 - Myelocytes',
      #                                                                            '8 - B lineage',
      #                                                                            '9 - T/NK lineage',
      #                                                                            '10 - Stromal'))#,
      #                                                                            # '11 - others'))
      
      #Position tweaks
      # mids[mids$finalAnn=='Leukocyte','UMAP_2'] = mids[mids$finalAnn=='Leukocyte','UMAP_2'] + 4.9
      # mids[mids$finalAnn=='Leukocyte','UMAP_1'] = mids[mids$finalAnn=='Leukocyte','UMAP_1'] - 0.3
      # mids[mids$finalAnn=='PTC','UMAP_2'] = mids[mids$finalAnn=='PTC','UMAP_2'] - 4.9
      # mids[mids$finalAnn=='PTC','UMAP_1'] = mids[mids$finalAnn=='PTC','UMAP_1'] -0.8
      # mids[mids$finalAnn=='Tumour','UMAP_2'] = mids[mids$finalAnn=='Tumour','UMAP_2'] + 3
      # mids[mids$finalAnn=='Tumour','UMAP_1'] = mids[mids$finalAnn=='Tumour','UMAP_1'] + 0.9
      # 
      boxed.labels(mids$UMAP_1,mids$UMAP_2,
                   labels=mids$label,cex = 0.6,xpad = 1.3,ypad = 2.3,border = T,
                   bg=celltype_cols[mids$label],
                   col='black')
      
      legend(x=-5, y=8,legend=unique(mids$label),fill = celltype_cols[unique(mids$label)],lwd = 0,cex = 0.75,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
    }
    
  }
  
  saveFig(file.path(plotDir,'SuppFig_a_pPTC_snRNAseq_UMAP'),plotFun_celltype,rawData=mdat,width = 4.1,height = 4,res = 500)
  
  
  
  
  plotFun_donorID = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0,0,0.8,0))
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','Paediatric PTC'),
         frame.plot=F)
    
    if(!noPlot){
      donor_cols = c('Y24' = '#ecbdc4','Y46' = '#b4d3b2')
      
      points(dd$UMAP_1,dd$UMAP_2,
             col = donor_cols[as.character(dd$donor)],
             pch = 19,
             cex=0.07)
      
    }
  }
  
  saveFig(file.path(plotDir,'SuppFig_a_pPTC_snRNAseq_UMAP_donor'),plotFun_donorID,rawData=mdat,width = 4.1,height = 4,res = 500)
  
  
  plotFun_sampling_site = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0,0,0.8,0))
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','Paediatric PTC'),
         frame.plot=F)
    
    if(!noPlot){
      sampling_site_cols = c('Normal' = '#afd5f0','Tumour'=colAlpha(grey(0.1),0.2),'Met' = grey(0.6))
      
      points(dd$UMAP_1,dd$UMAP_2,
             col = sampling_site_cols[as.character(dd$etiology)],
             pch = 19,
             cex=0.07)
      
    }
  }
  
  saveFig(file.path(plotDir,'SuppFig_a_pPTC_snRNAseq_UMAP_samplingSite'),plotFun_sampling_site,rawData=mdat,width = 4.1,height = 4,res = 500)

}


# SuppFigure X: Immune compartment composition
suppfig_x_pPTC_immune_composition = function(){
  if(file.exists(file.path(plotDir,'suppfig_a_pPTC_snRNAseq_UMAP_rawData.tsv'))){
    dd = read.delim(file.path(plotDir,'suppfig_a_pPTC_snRNAseq_UMAP_rawData.tsv'),header = T,sep = '\t')
    
  }else{
    mdat = read.csv(mdat_fp,row.names = 1)
    
    dd = mdat[mdat$annot %in% c('B_cells','DC1','Mast_cells','Monocytes','Plasma_cells','T_cells'),]
    
    dd = dd %>% group_by(donor,etiology,annot) %>% summarise(n_cell = n()) %>% 
      group_by(donor,etiology) %>% mutate(total_imm_cells = sum(n_cell),frac = n_cell / total_imm_cells)
    dd$annot = gsub('_',' ',dd$annot)
    dd$group = paste0(dd$donor,'_',dd$etiology)
    dd$group = factor(dd$group,c('Y24_Normal','Y46_Normal','Y24_Tumour','Y46_Tumour','Y46_Met'))
    table(is.na(dd$group))
  }
  plotFun_immune_composition = function(noFrame=FALSE,noPlot=FALSE){
    etiology_cols = c('Normal' = '#afd5f0','Tumour'=grey(0.1),'Met' = grey(0.6))
    dd$annot = factor(dd$annot,c('T cells','B cells','Plasma cells','Monocytes','DC1','Mast cells'))
    p= ggplot(dd,aes(group,frac))+
      geom_col(aes(fill = etiology))+
      geom_text(aes(label = n_cell), vjust = 0.5,hjust=-0.2,angle=90,size=3)+
      scale_y_continuous(breaks = seq(0,1,0.25),labels = c(0,0.25,0.5,0.75,1),limits = c(0,1))+
      scale_fill_manual(values = etiology_cols,name='Etiology')+
      facet_wrap(vars(annot),nrow=1)+
      xlab('')+ylab('Fraction of immune population per sample')+
      theme_classic()+
      theme(panel.border = element_rect(fill=F,color='black'),
            axis.line = element_blank(),axis.text.x =element_text(angle=90,vjust = 0.5,hjust = 1),
            strip.background.x = element_blank(),
            strip.text = element_text(color='black'),
            axis.text = element_text(color='black'),axis.ticks = element_line(color='black'))
    
    print(p)
    
  }
  
  saveFig(file.path(plotDir,'SuppFig_X_pPTC_immune_composition'),plotFun_immune_composition,rawData=dd,width = 9,height = 4,res = 500)
  
}




#####------------------------------------------------------------------------------------------
markers = c(
  "TSHR", "NKX2-1", "PAX8", "GLIS3", "TG","SLC26A4","IYD", "HHEX", "FOXE1", "DUOXA1", "DUOXA2", "DUOX1","DUOX2", "SLC5A5", "ZNF804B", 
  # Down-reg in tumour
  "TPO", "COL23A1", "PPARGC1A", "SLC5A8", "DIO2", "TFF3", 
  # up-reg in tumour
  "LRRK2", "HMGA2", "LMO3", "MET", "VAV3", "FN1", "LGALS3", "SERPINA1",  
  # Endothelium
  "FLT1", "PLVAP", "MECOM", "VWF","ANO2", "SLCO2A1", "CDH5", "ECSCR", "TBX1", "FLT4", "PROX1", 
  # Mesenchyme
  "CDH11","LAMA2", "COL6A3", "PDGFRA", "FBLN1", "BICC1", "MGP", 
  # Smooth muscle cells
  "GJC1", "PDGFRB","CLMN", 
  # Immune cells
  "MS4A1", "CD3E", "CD8A", "CD14", "TPSAB1", "CPA3", "KIT", "PTPRC"
)

# SuppFigure X: Dot Plot of canonical markers
suppfig_x_pPTC_celltype_dotplot = function(){
  srat = readRDS(srat_fp)
  srat$donor[srat$donor == 'Y46'] = 'P2'
  srat$donor[srat$donor == 'Y24'] = 'P1'
  srat$annot[srat$annot == 'Tumour' & srat$etiology == 'Met'] = 'Tumour (Met.)'
  srat$annot[srat$annot == 'Thyrocytes'] = paste0('Thyrocytes_',srat$donor[srat$annot == 'Thyrocytes'])
  srat$annot[grepl('Tumour',srat$annot)] = paste0(srat$annot[grepl('Tumour',srat$annot)],'_',srat$donor[grepl('Tumour',srat$annot)])
  
  srat$annot = factor(srat$annot,
                      c('Tumour_P1','Tumour_P2','Tumour (Met.)_P2','Thyrocytes_P1','Thyrocytes_P2',
                        'VECs','LECs','Mesenchymal','Fibroblasts','SMCs',
                        'B_cells','Plasma_cells','T_cells','DC1','Monocytes','Mast_cells',
                        'doublets',  'lowQual',  'unknown'))
  Idents(srat) = srat$annot
  plotFun_pPTC_dotplot = function(noFrame=FALSE,noPlot=FALSE){
    p = DotPlot(srat,idents = unique(Idents(srat)[!Idents(srat) %in% c('doublets',  'lowQual',  'unknown')]),
                group.by = 'annot',
                cols = c(colAlpha(grey(0.95),0.8),'black'),
                features = markers) + 
      coord_flip() + RotatedAxis()+
      theme(axis.text.x = element_text(size=10,vjust = 0.5,hjust = 1,angle = 90),
            axis.text.y = element_text(size=8,face = "italic"),
            legend.position = 'left',
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 9),
            legend.key.size = unit(0.4, "cm")
            )+xlab('')+ylab('')
      
    print(p)
  }
  saveFig(file.path(plotDir,'FigXX_pPTC_snRNAseq_dotplot'),plotFun_pPTC_dotplot,width = 5,height = 10,res = 500)  
}

  



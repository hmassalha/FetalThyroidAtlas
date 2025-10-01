
##------------------------------------------------------------------##
##  ADULT - bulkRNAseq DEG between REBC adult and paediatric PTC  ####
##------------------------------------------------------------------##

##---- Import bulk RNAseq data
# subset to just normal and RET-fusion PTC samples
rebc_se = readRDS('~/FetalThyroidAtlas/Data/published_bulkRNAseq/REBC-THYR/REBC_THYR_2508.RData')
rebc_rawCnt = assays(rebc_se)[['unstranded']]
rebc_mdat = as.data.frame(colData(rebc_se))
# Keep only samples of the relevant cancer types only
# According to this paper here: https://pmc.ncbi.nlm.nih.gov/articles/PMC10298340/
samples_toKeep = rebc_mdat$File.ID

rebc_mdat = as.data.frame(colData(rebc_se))
rebc_mdat$sampleID= rebc_mdat$File.ID
rebc_mdat$source = 'REBC_THYR'
rebc_mdat$sampleName = rownames(rebc_mdat)
rebc_mdat$cancerType = ifelse(rebc_mdat$Tissue.Type == 'Normal','Normal',
                              paste0('PTC_',rebc_mdat$Tumor.Descriptor))
rebc_mdat$age = rebc_mdat$AGE_SURGERY
rebc_mdat$sex = rebc_mdat$SEX

rebc_mdat$source[rebc_mdat$source == 'REBC_THYR' & as.numeric(rebc_mdat$age) <= 16] = 'REBC_THYR_paed'
rebc_mdat$source[rebc_mdat$source == 'REBC_THYR' & as.numeric(rebc_mdat$age) > 16] = 'REBC_THYR_adult'

rebc_mdat = rebc_mdat[rebc_mdat$Tissue.Type !='Normal',]
rebc_rawCnt = rebc_rawCnt[,rebc_mdat$sampleID]
rebc_mdat$group = rebc_mdat$source

##---- Perform edgeR
rebc_ptc_paed_vs_adult = fit_model(pb=rebc_rawCnt,
                               colDat=rebc_mdat,
                               formula='~ %s + sex',
                               geneMap=geneMap,groupID='group',
                               MDS_groups = c('cancerType','group','sex','source'),
                               pb_groupID='sampleID',coef=2,mycontrast=NULL)

saveRDS(rebc_ptc_paed_vs_adult,'rebc_PTC_paed_vs_adult_edgeR.RDS')
rebc_ptc_paed_vs_adult = readRDS('rebc_PTC_paed_vs_adult_edgeR.RDS')

# Extract DEGs
deg = rebc_ptc_paed_vs_adult[['tt']]
deg = deg[abs(deg$logFC) >= 1.5 & deg$FDR < 0.05 & 
            #deg$logCPM > -0.5 & 
            !grepl('^RPL|^RPS|^MT-|^MALAT1|^NEAT1|^GAPDH|^GABARAP|^MRPL|^MRPS|AC\\d+|LINC\\d+',deg$geneSym),]
deg$direction = ifelse(deg$logFC > 0,'aPTC_down','aPTC_up')
table(deg$direction)
deg = deg[order(abs(deg$logFC),decreasing = T),]


moduleList[['rebc_aPTC_up']] = deg$ensID[deg$direction == 'aPTC_up']
moduleList[['rebc_aPTC_down']] = deg$ensID[deg$direction == 'aPTC_down']


plotFun_fTFC_moduleScore = function(noFrame=FALSE,noPlot=FALSE){
  
  allScore$group_facet_ver = allScore$moduleType
  
  dd = allScore[grepl('FFPE|Thyrocytes|Tumour|fTFC|fThy|metastatic|Normal|aTFC|C\\d|PTC',allScore$cancerType) & 
                  #!grepl('follicular|tallCell|metastatic|primary',allScore$cancerType) &
                  !grepl('follicular|tallCell',allScore$cancerType) &
                  #!grepl('primary',allScore$cancerType_details) &
                  allScore$source %in% c('Sanger','TCGA_Thyroid',#'Yoo_2016',
                                         'He_2021','Lee_2024','REBC_THYR_paed','REBC_THYR_adult') &
                  allScore$moduleType %in% moduleType_toUse,]
  
  dd$ageGroup = ifelse(dd$source == 'Sanger',dd$ageCat,'adult')
  dd$source[dd$age]
  dd = dd[dd$ageGroup != 'foetus',]
  dd$cancerNormal = ifelse(dd$cancerType %in% c('Normal'),'Normal',
                           ifelse(dd$cancerType == 'Normal.adj','Normal.adj','Tumour'))
  dd$cancerNormal = factor(dd$cancerNormal,c('Normal','Normal.adj','Tumour'))
  dd$med_normal = NA
  for(dataset in unique(dd$source)){
    for(mod in unique(dd$moduleType)){
      tmp = dd[dd$source == dataset & dd$moduleType == mod,]
      med_normal = median(tmp$TotalScore[tmp$cancerNormal == 'Normal'])
      dd$med_normal[dd$source == dataset & dd$moduleType == mod] = med_normal
    }
  }
  
  dd$normalised_score = dd$TotalScore - dd$med_normal
  dd$source = factor(dd$source,c('Sanger','TCGA_Thyroid','He_2021','Lee_2024','REBC_THYR_paed', 'REBC_THYR_adult'))
  
  table(dd$cancerNormal,dd$cancerNormal,dd$source)
  
  dd$cancerType[dd$source %in% c('REBC_THYR_paed','REBC_THYR_adult') & 
                  dd$cancerNormal != 'Normal' & 
                  dd$sampleID %in% rebc_mdat$File.ID[grepl('RET',rebc_mdat$WGS_CandidateDriverFusion)]] = paste0(dd$cancerType[dd$source %in% c('REBC_THYR_paed','REBC_THYR_adult') & 
                                                                                                                                 dd$cancerNormal != 'Normal' & 
                                                                                                                                 dd$sampleID %in% rebc_mdat$File.ID[grepl('RET',rebc_mdat$WGS_CandidateDriverFusion)]],'_RET')
  
  dd$cancerType[grepl('RET',dd$cancerType)] = 'RET'
  dd$cancerType[!grepl('RET',dd$cancerType)] = dd$cancerNormal[!grepl('RET',dd$cancerType)]
  
  # ## Remove REBC-THYR BRAF samples
  # samples_to_remove = rebc_mdat$File.ID[grepl('BRAF',rebc_mdat$WGS_CandidateDriverFusion) | 
  #                                         grepl('BRAF',rebc_mdat$WGS_CandidateDriverMutation)]
  # dd = dd[!dd$sampleID %in% samples_to_remove,]
  
  # dd$cancerType_details[grepl('REBC',dd$source)] = paste0(dd$cancerType_details[grepl('REBC',dd$source)],'::',rebc_mdat$WGS_CandidateDriverFusion[match(dd$sampleID[grepl('REBC',dd$source)],rebc_mdat$File.ID)])
  # dd$cancerType_details[grepl('REBC',dd$source)] = paste0(dd$cancerType_details[grepl('REBC',dd$source)],'::',rebc_mdat$WGS_CandidateDriverMutation[match(dd$sampleID[grepl('REBC',dd$source)],rebc_mdat$File.ID)])
  # dd$cancerType_details[grepl('REBC',dd$source) & !grepl('Normal::',dd$cancerType_details) & grepl('NCOA4-RET \\(RET-PTC3\\)',dd$cancerType_details)] = 'RET_PTC3'
  # dd$cancerType_details[grepl('REBC',dd$source) & !grepl('Normal::',dd$cancerType_details) & grepl('\\(RET-PTC\\d+\\)',dd$cancerType_details)] = 'RET_PTC_others'
  # dd$cancerType_details[grepl('REBC',dd$source) & !grepl('Normal::',dd$cancerType_details) & grepl('::::.*BRAF',dd$cancerType_details)] = 'BRAF_PTC'
  # dd$cancerType_details[grepl('REBC',dd$source) & grepl('Normal::',dd$cancerType_details)] = 'Normal'
  # dd$cancerType_details[grepl('REBC',dd$source) & grepl('PTC_Metastatic::',dd$cancerType_details)] = 'PTC_Metastatic::others'
  # dd$cancerType_details[grepl('REBC',dd$source) & grepl('PTC_Primary::',dd$cancerType_details)] = 'PTC_Primary::others'
  
  
  #[!dd$sampleName %in% sample_metadata$geo_accession[grepl('_P|-P',sample_metadata$title)],]
  p1 = ggplot(dd, aes(cancerNormal, normalised_score)) +
    geom_hline(yintercept = 0,linetype=2,linewidth=0.3)+
    geom_quasirandom(size=0.4,width = 0.15,alpha=0.6,aes(col=cancerType))+
    geom_boxplot(aes(fill=cancerNormal),outlier.shape = NA,position = 'dodge', alpha = 0.8,width=0.5,linewidth=0.3,colour='black') +
    #geom_point(data=dd[dd$cancerType == 'RET',],size=4)+
    #geom_point(data=dd[dd$sampleName %in% sample_metadata$geo_accession[grepl('_P|-P',sample_metadata$title)],],col='red',size=2)+
    scale_fill_manual(values =c('Normal'=grey(0.8),'Normal.adj'=grey(0.4),'Tumour'='#511378'))+
    #scale_fill_manual(values =c(col25,pal34H))+
    #scale_color_manual(values =c(col25,pal34H))+
    
    #scale_fill_manual(values = c(rep(col25[4],2),'#c7065a',col25[4],rep(colAlpha(col25[1],0.4),3),rep(grey(0.7),5))) +
    #scale_fill_manual(values = c(rep(col25[4],2),rep(grey(0.7),7))) +
    facet_grid(group_facet_ver~source,scales = 'free',space = 'free_x')+
    theme_classic()+
    #ylim(-0.1,0.1)+
    #ggtitle(title)+
    xlab('')+ylab('Centralised fTFC signature score')+
    theme(panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
          strip.background=element_rect(linewidth=0),
          axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1,colour = 'black'))
  
  print(p1)
}


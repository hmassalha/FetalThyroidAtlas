# Supplementary table S12 - list of PTC / normal thyroid bulk RNA-seq samples used in the study

setwd('~/FetalThyroidAtlas/')

#----------------##
#   Libraries  ####
#----------------##
library(tidyverse)
library(GenomicFeatures)
source("R/utils/misc.R")

##----------------------------##
##   Set Global parameters  ####
##----------------------------##
outDir = "SupplementaryTables"
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}


##--- Import bulk counts --------
bulk_sources = c('inhouse'='Data/inhouse_bulk/inhouse_bulkRNA_fetalThyroid_paedPTC.RDS',
                 'TCGA_Thyroid'='Data/published_bulkRNAseq/TCGA_Thyroid/TCGA_Thyroid_bulkRNA_se.RDS',
                 'REBC_THYR' = 'Data/published_bulkRNAseq/REBC-THYR/REBC_THYR_2508.RData',
                 'He2021'='Data/published_bulkRNAseq/He_etal_21/aPTC_He_2021_se.RDS',
                 'Lee2024' = 'Data/published_bulkRNAseq/Lee_etal_24/aPTC_Lee_2024_se.RDS')

columns_to_keep = c('source', 'sample_id', 'donor_id', 'age_category', 'sex', 'cancer_normal', 'cancer_type', 'reference_genome', 'accession_code', 'WGS_sample_id', 'WGS_coverage', 'WGS_accession_code')

mdat_list = list()

for(i in names(bulk_sources)){
  se_obj = readRDS(bulk_sources[[i]])
  mdat = as.data.frame(SummarizedExperiment::colData(se_obj))
  if(i == 'inhouse'){
    mdat = mdat %>% 
      dplyr::mutate(sample_id = sangerSampleID,
                    donor_id = dplyr::case_when(DonorID == 'PR62328' ~ 'PD62328 (P1)',
                                                DonorID == 'PR64651' ~ 'PD64651 (P2)',
                                                .default = gsub('PR','PD',DonorID)),
                    age_category = dplyr::case_when(Sample_owner == 'Hassan' ~ paste0('fetus (',gsub(' .*$','',other_sampleID),')'),
                                                    .default = 'paediatric'),
                    sex = ifelse(is.na(sex),'-',sex),
                    cancer_normal = cancerType,
                    cancer_type = dplyr::case_when(SampleID == 'PR62328a' ~ 'PTC - Primary',
                                                   SampleID == 'PR64651c' ~ 'PTC - Primary (left)',
                                                   SampleID == 'PR64651d' ~ 'PTC - Primary (right)',
                                                   SampleID == 'PR64651e' ~ 'PTC - Metastatic (lymph_node)',
                                                   .default = cancer_normal),
                    reference_genome = "GRCh38 2020-A",
                    accession_code = "EGAD00001015448",
                    WGS_sample_id = dplyr::case_when(SampleID %in% c('PR62328a','PR62328b','PR64651c','PR64651d','PR64651e') ~ gsub('PR','PD',SampleID),
                                                     .default = '-'),
                    WGS_coverage = dplyr::case_when(SampleID == 'PR62328a' ~ '40X',
                                                    SampleID == 'PR62328b' ~ '38X',
                                                    SampleID == 'PR64651d' ~ '91X',
                                                    SampleID == 'PR64651e' ~ '211X',
                                                    .default = '-'),
                    WGS_accession_code = dplyr::case_when(WGS_sample_id != '-' ~ 'EGAD00001015447',
                                                          .default = '-')) %>% 
      dplyr::select(all_of(columns_to_keep))
    
    extra_row_for_P2_blood_WGS <- tibble(
      source = 'Sanger',
      sample_id = "-",
      donor_id = "PD64651 (P2)",
      age_category = "paediatric",
      sex = 'Female',
      cancer_normal = "Normal",
      cancer_type = "Normal (blood)",
      reference_genome = "GRCh38 2020-A",
      accession_code = "-",
      WGS_sample_id = "PD64651a",
      WGS_coverage = "107X",
      WGS_accession_code = "EGAD00001015447"
    )
    
    mdat <- dplyr::bind_rows(mdat, extra_row_for_P2_blood_WGS)
  }
  
  if(i == 'He2021'){
    mdat = mdat %>% 
      dplyr::mutate(sample_id = paste0(sampleID,' (',title,')'),
                    donor_id = gsub('_.*$','',description),
                    age_category = 'adult',
                    sex = dplyr::case_when(sex == 'F' ~ 'Female',
                                           sex == 'M' ~ 'Male',
                                           .default = 'others'),
                    cancer_normal = cancerType,
                    cancer_type = cancerType,
                    reference_genome = "GRCh38p7",
                    accession_code = "GSE165724",
                    WGS_sample_id = '-',
                    WGS_coverage = '-',
                    WGS_accession_code = '-') %>% 
      dplyr::select(all_of(columns_to_keep))
  }
  
  if(i == 'Lee2024'){
    mdat = mdat %>% 
      dplyr::filter(cell.type.ch1 %in% c('Normal','PTC')) %>% 
      dplyr::mutate(sample_id = paste0(sampleID,' (',title,')'),
                    donor_id = gsub('-N.*|N.*|-T.*|T.*','',title),
                    age_category = 'adult',
                    sex = '-',
                    cancer_normal = ifelse(cell.type.ch1 == 'Tumor','PTC','Normal'),
                    cancer_type = dplyr::case_when(cancer_normal == 'Normal' ~ 'Normal',
                                                   .default = cell.subtype.ch1),
                    reference_genome = "GRCh38",
                    accession_code = "GSE213647",
                    WGS_sample_id = '-',
                    WGS_coverage = '-',
                    WGS_accession_code = '-') %>% 
      dplyr::select(all_of(columns_to_keep))
  }
  
  if(i == 'TCGA_Thyroid'){
    mdat = mdat %>% 
      dplyr::filter(grepl('Papillary carcinoma|Papillary adenocarcinoma',primary_diagnosis) &
                      tissue_or_organ_of_origin %in% c('Thyroid gland',"Lymph node, NOS")) %>% 
      dplyr::mutate(source='TCGA_Thyroid',
                    sample_id = barcode,
                    donor_id = patient,
                    age_category = 'adult',
                    sex = dplyr::case_when(gender == 'male' ~ 'Male',
                                           gender == 'female' ~ 'Female',
                                           .default = 'others'),
                    cancer_normal = ifelse(tissue_type == 'Tumor','PTC',tissue_type),
                    cancer_type = dplyr::case_when(sample_type == 'Metastatic' ~ 'PTC - Metastatic',
                                                   sample_type == 'Primary Tumor' ~ 'PTC - Primary',
                                                   sample_type == 'Solid Tissue Normal' ~ 'Normal',
                                                   .default = 'others'),
                    reference_genome = "GRCh38.d1.vd1",
                    accession_code = "TCGA GDC data portal",
                    WGS_sample_id = '-',
                    WGS_coverage = '-',
                    WGS_accession_code = '-') %>% 
      dplyr::select(all_of(columns_to_keep))
  }
  
  if(i == 'REBC_THYR'){
    mdat = mdat %>% 
      dplyr::mutate(source=Project.ID,
                    sample_id = Sample.ID,
                    donor_id = Case.ID,
                    age_category = ifelse(AGE_SURGERY <= 16, 'paediatric','adult'),
                    sex = dplyr::case_when(SEX == 'male' ~ 'Male',
                                           SEX == 'female' ~ 'Female',
                                           .default = 'others'),
                    cancer_normal = ifelse(Tissue.Type == 'Tumor','PTC',Tissue.Type),
                    cancer_type = ifelse(Tissue.Type == 'Tumor',paste0('PTC - ',Tumor.Descriptor),Tissue.Type),
                    reference_genome = "GRCh38.d1.vd1",
                    accession_code = "GDC data portal",
                    WGS_sample_id = '-',
                    WGS_coverage = '-',
                    WGS_accession_code = '-') %>% 
      dplyr::select(all_of(columns_to_keep))
    
  }

  mdat_list[[i]] = mdat
}

bulk_samples = do.call(rbind,mdat_list)
checkmate::assert_false('others' %in% c(bulk_samples$sex,bulk_samples$cancer_type))
checkmate::assert_true(all(colSums(is.na(bulk_samples)) == 0))

## reorder the table
bulk_samples$source = factor(bulk_samples$source,c('Sanger','He_2021','Lee_2024','TCGA_Thyroid','REBC-THYR'))
bulk_samples$donor_id = as.factor(bulk_samples$donor_id)
bulk_samples$cancer_type = factor(bulk_samples$cancer_type,c("Normal", "Normal.adj", "Normal (blood)", "PTC", "PTC - Primary", "PTC - Primary (left)", "PTC - Primary (right)", 
                                                             "PTC - Metastatic", "PTC - Metastatic (lymph_node)"))

bulk_samples = bulk_samples[order(bulk_samples$cancer_type),]
bulk_samples = bulk_samples[order(bulk_samples$donor_id),]
bulk_samples = bulk_samples[order(bulk_samples$source),]

write.csv(bulk_samples,file.path(outDir,'Supplementary_Table_S12_bulkRNAseq_samples.csv'),row.names = F)

## Check that it matches the raw data for figure 5C
fig5C_rawdata = read.delim('Figures/2508/Fig4c_fTFC1.2_moduleScore_bulkSamples_sub_rawData.tsv')
fig5C_rawdata = fig5C_rawdata[fig5C_rawdata$moduleType == 'fTFC2',]
table(fig5C_rawdata$source,fig5C_rawdata$moduleType)
a = table(bulk_samples$source)
a['Sanger'] = a['Sanger'] - 1
b = table(fig5C_rawdata$source[fig5C_rawdata$moduleType == 'fTFC2'])
b['REBC-THYR'] = b['REBC_THYR_adult'] + b['REBC_THYR_paed']
b = b[match(names(a),names(b))]

checkmate::assert_true(all(a == b))

table(fig5C_rawdata$source,fig5C_rawdata$cancerType)
bulk_samples$source2 = ifelse(bulk_samples$source == 'REBC-THYR',paste0('REBC-THYR (',bulk_samples$age_category,')'),bulk_samples$source)
table(bulk_samples$cancer_type,bulk_samples$source2)

table(bulk_samples$age_category,bulk_samples$cancer_normal)

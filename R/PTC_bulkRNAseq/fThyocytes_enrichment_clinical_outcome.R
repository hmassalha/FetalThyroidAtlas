# Correlate fTFC1/2 score with TCGA clinical output 

# Libraries --------------------------------------------------------------------
library(tidyverse)

# Set up -----------------------------------------------------------------------
setwd('~/FetalThyroidAtlas/')

outDir = "Results/2505/PTC_bulkRNAseq"
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

plotDir = 'Figures/2505'

# Import datasets --------------------------------------------------------------
## fTFC1/fTFC2 signature in bulk RNAseq data -----------------------------------
allScore = read.csv(file.path(outDir,'fTFC1.2_top100_geneSignatures_SingScore_bulkRNAseq.csv'))
## TCGA metadata ---------------------------------------------------------------
tcga_se_path <- 'Data/published_bulkRNAseq/TCGA_Thyroid/TCGA_Thyroid_bulkRNA_se.RDS'
tcga_se = readRDS(tcga_se_path)
tcga_mdat = as.data.frame(colData(tcga_se))
tcga_mdat$sampleID= rownames(tcga_mdat)
tcga_mdat$source = 'TCGA_Thyroid'
tcga_mdat$sampleName = rownames(tcga_mdat)
tcga_mdat$cancerType = ifelse(tcga_mdat$tissue_type == 'Normal','Normal',
                              paste0('PTC_',tcga_mdat$classification_of_tumor))
tcga_mdat$age = tcga_mdat$age_at_diagnosis
tcga_mdat$sex = tcga_mdat$gender


table(tcga_mdat$sampleID %in% allScore$sampleID)


# Correlation analysis ---------------------------------------------------------
samples_to_keep = intersect(allScore$sampleID,tcga_mdat$sampleID)
tcga_data = merge(tcga_mdat[match(samples_to_keep, tcga_mdat$sampleID),],
                  allScore[allScore$sampleID %in% samples_to_keep,
                           c("sampleID",colnames(allScore)[!colnames(allScore) %in% colnames(tcga_mdat)])],by='sampleID',all=TRUE)

checkmate::assert_true(nrow(tcga_data) == 526*n_distinct(allScore$moduleType))
tcga_data$definition
tcga_data$tumor_descriptor
table(tcga_data$ajcc_pathologic_stage)
table(tcga_data$prior_treatment)
table(tcga_data$diagnosis_is_primary_disease)
table(tcga_data$primary_diagnosis)
table(tcga_data$vital_status)
table(tcga_data$paper_Risk)
table(tcga_data$paper_medical_history_thyroid)
tcga_data$days_to_death

table(tcga_data$paper_armDriver)
# Tumour sub-types
"definition", "tumor_descriptor","sample_type",
"ajcc_pathologic_stage","synchronous_malignancy","ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m"
# survival metrics 
"vital_status", "days_to_diagnosis", "days_to_last_follow_up"   
"age_at_diagnosis",
"year_of_diagnosis","days_to_death"
# clinical treatment
"treatments","primary_diagnosis", "prior_malignancy","prior_treatment","diagnosis_is_primary_disease",
"residual_disease"                            
"classification_of_tumor", "tumor_focality"                              
# patient metadata
"sex", "age"
# Driver mutation
"paper_BRAF","paper_RAS","paper_Driver"                                
[72] "paper_mutDriver"                             
[74] "paper_mutRareDriver"                         
[76] "paper_armDriver"                             
[77] "paper_armDriver_CN"                          
[78] "paper_focalDriver"                           
[79] "paper_focalDriver_CN"                        
[80] "paper_fusionDriver"                          
[81] "paper_fusionDriverGenes"                     
[88] "paper_fusionDriverVote"                      
[89] "paper_purity"                                
[90] "paper_ploidy"                                
[91] "paper_Genome_doublings"                      
[93] "paper_nonsil_density"                        
[94] "paper_nnon_tot"                              
[95] "paper_nsil_tot"                              
[96] "paper_ndbsnp_tot"                            
[97] "paper_Mb_cov"                                
[98] "paper_nondriver_density"                     
[99] "paper_nSNP"                                  
[100] "paper_nindel"                                
[101] "paper_mut_density"                           
[102] "paper_nmut"                                  
[103] "paper_nmut_driver"            
[109] "paper_pair_type"                             
[110] "paper_Risk"                                  
[111] "paper_MACIS"                                 
[112] "paper_histological_type"                     
[113] "paper_histological_typeFG"                   
[114] "paper_Follicular_fraction"                   
[115] "paper_age"                                   
[116] "paper_braf_genotype_lab"                     
[117] "paper_gender"                                
[118] "paper_radiation_exposure_indicator"          
[132] "paper_BRAF_GROUP"                            
[133] "paper_BRAF_CLUSTER"                          
[134] "paper_BRAF_RAF_score"                        
[135] "paper_BRAF_RAF_class"                        
[136] "paper_Arm_SCNA_Cluster"                      
[137] "paper_Arm_SCNA_Cluster_number"               
[138] "paper_mRNA_Cluster_number"                   
[139] "paper_miRNA_Cluster_number"                  
[140] "paper_meth_Cluster"                          
[141] "paper_meth_Cluster_number"                   
[142] "paper_RPPA_Cluster_number"                   
[143] "paper_SomaticRearrangments"                  
[144] "paper_germline_mutGenes"                     
[145] "paper_germline_mutGenes_variant"
[151] "paper_pathologic_T"                          
[152] "paper_pathologic_N"                          
[153] "paper_pathologic_M"                          
[154] "paper_Neoplasm_Disease_Stage"                
[155] "paper_Extrathyroidal_extension"              
[156] "paper_N0vsN1b"                               
[191] "paper_differentiation_score"                 
[192] "paper_ERK_score"
[210] "paper_nonsil_density_age_fit_residual"       
[211] "paper_Follow_up"                             
[212] "paper_Most_Recent_Days_to_Follow_up"         
[213] "paper_Follow_up_New_Tumor_Event"             
[222] "paper_BRAFV600E_RAS"
[223] "paper_Clinical_Info_Freeze"                  
[230] "paper_SomaticRearrangment_story"             
[231] "paper_medical_history_thyroid"               

ggplot(tcga_data[tcga_data$moduleType %in% c('fTFC1','fTFC2'),],
       aes(reorder(paper_BRAF,TotalScore,median),TotalScore))+
  geom_quasirandom(
    width = 0.15,
    alpha = 1,
    size = 0.7
  ) +
    geom_boxplot(
      width = 0.65,
      outlier.shape = NA,
      alpha = 0.7,
      colour = "black",
      linewidth = 0.4
    ) +
    facet_grid(
       ~ moduleType + cancerNormal,
      scales = "free_x",
      space = "free"
    ) +
    scale_fill_manual(values = age_cols) +
    labs(
      x = NULL,
      y = "TDS enrichment score"
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      ),
      axis.line = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(
        colour = "black",
        size = 9
      ),
      axis.text = element_text(
        colour = "black",
        size = 10
      ),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      # axis.text.x = element_blank(),
      # panel.spacing.x = unit(0.7, "cm"),
      legend.position = "none"
    )


library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(checkmate)
auto_numeric_convert <- function(x, threshold = 0.9){
  
  ## preserve original class
  if(is.numeric(x)) return(x)
  
  if(is.factor(x))
    x <- as.character(x)
  
  if(!is.character(x))
    return(x)
  
  ## decimal comma fix
  x_clean <- gsub(",", ".", x[!is.na(x) & x!=''], fixed = TRUE)
  
  ## attempt conversion
  x_num <- suppressWarnings(as.numeric(x_clean))
  
  ## proportion successfully converted
  prop_numeric <- mean(
    is.na(x_clean) | !is.na(x_num)
  )
  
  ## convert only if mostly numeric
  if(prop_numeric >= threshold){
    x_num_full <- as.numeric(gsub(",", ".", x, fixed = TRUE))
    return(x_num_full)
  }
  
  return(x)
}
auto_detect_plot_type <- function(x,
                                  max_levels_numeric = 5,
                                  numeric_threshold = 0.9){
  
  x <- auto_numeric_convert(x)
  
  ## character / factor → categorical
  if(is.character(x) || is.factor(x)){
    return(list(
      type = "categorical",
      data = x
    ))
  }
  
  ## numeric handling
  if(is.numeric(x)){
    
    n_unique <- data.table::uniqueN(
      stats::na.omit(x)
    )
    
    ## binary / small discrete numeric
    if(n_unique <= max_levels_numeric){
      
      return(list(
        type = "categorical",
        data = factor(x)
      ))
    }
    
    ## continuous numeric
    return(list(
      type = "continuous",
      data = x
    ))
  }
  
  return(list(
    type = "categorical",
    data = factor(x)
  ))
}

plot_tcga_associations <- function(
  data,
  vars,
  score_col = "TotalScore",
  module_keep = c("fTFC1","fTFC2"),
  outfile = "TCGA_clinical_associations.pdf"
){
  
  checkmate::assert_data_frame(data)
  checkmate::assert_character(vars,min.len=1)
  checkmate::assert_string(score_col)
  
  dt <- data.table::as.data.table(data)
  
  checkmate::assert_subset(
    c(score_col,"moduleType"),
    colnames(dt)
  )
  
  vars <- vars[vars %in% colnames(dt)]
  
  dt <- dt[moduleType %in% module_keep]
  
  pdf(
    outfile,
    width = 12,
    height = 8,
    onefile = TRUE
  )
  
  on.exit(dev.off())
  
  for(vv in vars){
    
    message("Processing: ", vv)
    
    dd <- copy(
      dt[
        !is.na(get(vv)) &
          !is.na(get(score_col))
        ]
    )
    
    if(nrow(dd) < 10){
      next
    }
    det <- auto_detect_plot_type(x = dd[[vv]])
    plot_type <- det$type
    x <- det$data
    
    
    try({
      
      ## -------------------------
      ## categorical variables
      ## -------------------------
      
      if(
        !is.numeric(x)
      ){
        
        dd[, plot_var := factor(
          get(vv)
        )]
        
        if(nlevels(dd$plot_var) < 2){
          next
        }
        
        p <- ggplot(
          dd,
          aes(
            reorder(plot_var,
                    get(score_col),
                    median,
                    na.rm=TRUE),
            get(score_col)
          )
        ) +
          
          ggbeeswarm::geom_quasirandom(
            width = 0.15,
            alpha = 0.45,
            size = 0.5
          ) +
          
          geom_boxplot(
            width = 0.65,
            outlier.shape = NA,
            alpha = 0.75,
            colour = "black",
            linewidth = 0.35
          ) +
          
          facet_grid(
            ~ moduleType + cancerNormal,
            scales = "free_x",
            space = "free"
          ) +
          
          labs(
            title = vv,
            x = NULL,
            y = "Enrichment score"
          ) +
          
          theme_classic(base_size = 12) +
          
          theme(
            panel.border = element_rect(
              fill = NA,
              colour = "black",
              linewidth = 0.5
            ),
            strip.background = element_blank(),
            axis.text.x = element_text(
              angle = 90,
              hjust = 1,
              vjust = 0.5,
              size = 8
            ),
            legend.position = "none"
          )
        
      } else {
        
        ## -------------------------
        ## numeric variables
        ## -------------------------
        
        suppressWarnings({
          
          cor_res <- cor.test(
            x,
            dd[[score_col]],
            method = "spearman"
          )
          
        })
        
        rho <- round(
          unname(cor_res$estimate),
          3
        )
        
        pval <- signif(
          cor_res$p.value,
          3
        )
        
        label_txt <- paste0(
          "Spearman ρ = ",
          rho,
          "\nP = ",
          pval
        )
        dd[[vv]] <- x
        p <- ggplot(
          dd,
          aes(
            get(vv),
            get(score_col)
          )
        ) +
          
          geom_point(
            alpha = 0.45,
            size = 0.8
          ) +
          
          geom_smooth(
            method = "lm",
            colour = "#D55E00",
            linewidth = 0.8
          ) +
          
          facet_grid(
            ~ moduleType + cancerNormal
          ) +
          
          annotate(
            "text",
            x = Inf,
            y = Inf,
            label = label_txt,
            hjust = 1.1,
            vjust = 1.2,
            size = 4
          ) +
          
          labs(
            title = vv,
            x = vv,
            y = "Enrichment score"
          ) +
          
          theme_classic(base_size = 12) +
          
          theme(
            panel.border = element_rect(
              fill = NA,
              colour = "black",
              linewidth = 0.5
            ),
            strip.background = element_blank()
          )
        
      }
      
      print(p)
      
    }, silent=TRUE)
  }
}

vars_of_interest <- c(
  
  "definition","tumor_descriptor","sample_type",
  "ajcc_pathologic_stage",
  "ajcc_pathologic_t",
  "ajcc_pathologic_n",
  "ajcc_pathologic_m",
  
  "vital_status",
  "days_to_last_follow_up",
  "days_to_death",
  "age_at_diagnosis",
  
  "sex",
  "age",
  
  "paper_BRAF",
  "paper_RAS",
  "paper_Driver",
  
  "paper_purity",
  "paper_ploidy",
  "paper_mut_density",
  
  "paper_MACIS",
  "paper_histological_type",
  "paper_BRAF_GROUP",
  "paper_BRAF_RAF_score",
  
  "paper_differentiation_score",
  "paper_ERK_score",
  
  "paper_Follow_up_New_Tumor_Event"
)
plot_tcga_associations(
  data = tcga_data,
  vars = vars_of_interest,
  score_col = "TotalScore",
  outfile = file.path(
    outDir,
    "TCGA_reviewer4_clinical_associations.pdf"
  )
)

##-----------------------------------------------------------------------##
##   TDS CORRELATION ANALYSIS - corrected standalone version             ##
##   Assumes in environment:                                             ##
##     bulk_ranked  : ranked TPM matrix (from singscore::rankGenes)     ##
##     bulk_samples : sample metadata dataframe                         ##
##     gene_map     : ensID <-> geneSymbol mapping                      ##
##     allScore     : raw singscore output (before local normalisation)  ##
##     outDir / plotDir : output directories                             ##
##-----------------------------------------------------------------------##

# Mi - temp files --------------------------------------------------------------

# Libraries --------------------------------------------------------------------
library(tidyverse)
library(singscore)
library(ggbeeswarm)
# source('R/helperFunctions.R')

# Set up -----------------------------------------------------------------------
setwd('~/FetalThyroidAtlas/')

outDir = "Results/2505/PTC_bulkRNAseq"
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = T)
}

plotDir = 'Figures/2505'

PATHS <- list(
  # fTFC12score_fp = file.path(
  #   outDir,
  #   'fTFC1.2_top100_geneSignatures_SingScore_bulkRNAseq.csv'
  # ),
  fTFC12score_fp = "../tmp_thyroid/fTFC1.2_top100_geneSignatures_SingScore_bulkRNAseq.csv",
  bulk_ranked_fp = "../tmp_thyroid/tmp_bulkranked_mtx.csv",
  bulk_sample_fp = "../tmp_thyroid/tmp_bulksamples.csv",
  gene_map_fp = "../tmp_thyroid/geneMap.csv",
  tdsscore_fp = "../tmp_thyroid/TDS_SingScore_bulkRNAseq.csv"
)

# Import relevant data ---------------------------------------------------------
## fTFC1/fTFC2 signature in bulk RNAseq data -----------------------------------
allScore = read.csv(PATHS$fTFC12score_fp, row.names = 1)

## bulk RNA-seq data -----------------------------------------------------------
# import bulk counts and calculate cpmCnt in xx01_moduleScoring.R
if (is.null(PATHS$bulk_ranked_fp) || is.na(PATHS$bulk_ranked_fp)) {
  bulkRNA = import_bulkRNA_thyroid(
    bulk_sources = c(
      'inhouse' = 'Data/inhouse_bulk/inhouse_bulkRNA_fetalThyroid_paedPTC.RDS',
      'TCGA_Thyroid' = 'Data/published_bulkRNAseq/TCGA_Thyroid/TCGA_Thyroid_bulkRNA_se.RDS',
      'REBC_THYR' = 'Data/published_bulkRNAseq/REBC-THYR/REBC_THYR_2508.RData',
      'He2021' = 'Data/published_bulkRNAseq/He_etal_21/aPTC_He_2021_se.RDS',
      'Lee2024' = 'Data/published_bulkRNAseq/Lee_etal_24/aPTC_Lee_2024_se.RDS'
    ),
    gene_map = gene_map
  )
  bulk_samples = bulkRNA[['bulk_samples']]
  tpmCnt = bulkRNA[['tpm_count']]
  # cpmCnt = bulkRNA[['cpmCnt']]
  # rawCnt = bulkRNA[['raw_count']]
  mtx = tpmCnt[, !colnames(tpmCnt) %in% c('ensID', 'geneLength')]
  # apply the rankGenes method
  bulk_ranked = rankGenes(mtx)

  # write.csv(bulk_ranked,file.path(outDir,'tmp_bulkranked_mtx.csv'),row.names = T)
  # write.csv(bulk_samples,file.path(outDir,'tmp_bulksamples.csv'),row.names = T)
} else {
  bulk_ranked = read.csv(PATHS$bulk_ranked_fp, row.names = 1) |>
    as.matrix()
  bulk_samples = read.csv(PATHS$bulk_sample_fp, row.names = 1)
}

## gene map --------------------------------------------------------------------
gene_map = read.csv(PATHS$gene_map_fp)

## TDS gene set ----------------------------------------------------------------
## TCGA core panel: https://www.cell.com/cell/fulltext/S0092-8674(14)01238-0
## These are the canonical thyroid differentiation genes
tds_genes_symbols <- c(
  "SLC5A5", # NIS
  "TPO",
  "TG",
  "TSHR",
  "SLC26A4", # pendrin
  "DIO1",
  "DIO2",
  "DUOX1",
  "DUOX2",
  "DUOXA1",
  "DUOXA2",
  "IYD",
  "NKX2-1",
  "FOXE1",
  "PAX8",
  "HHEX"
)

## Convert to ensIDs
tds_genes_ens <- gene_map$gene_id[gene_map$gene_name %in% tds_genes_symbols]
tds_genes_ens <- tds_genes_ens[!is.na(tds_genes_ens)]

## Check coverage
found <- tds_genes_symbols[tds_genes_symbols %in% gene_map$gene_name]
notfound <- tds_genes_symbols[!tds_genes_symbols %in% gene_map$gene_name]
cat(sprintf(
  "\nTDS genes found in gene_map: %d/%d\n",
  length(found),
  length(tds_genes_symbols)
))
if (length(notfound) > 0) {
  cat("Not found:", paste(notfound, collapse = ', '), "\n")
}

## Check these ensIDs are actually in bulk_ranked
in_matrix <- sum(tds_genes_ens %in% rownames(bulk_ranked))
cat(sprintf(
  "TDS ensIDs present in bulk_ranked: %d/%d\n",
  in_matrix,
  length(tds_genes_ens)
))


# Step 0: Safety checks before running -----------------------------------------
## Run these first and verify output looks sensible

cat("allScore columns:\n")
print(colnames(allScore))

cat("\nallScore moduleType values:\n")
print(table(allScore$moduleType))

cat("\nallScore cancerType values (first 20):\n")
print(head(table(allScore$cancerType), 20))

## Check whether normalised_score already exists globally
## (it probably doesn't - it's computed locally in your plot functions)
cat(
  "\n'normalised_score' exists in allScore:",
  'normalised_score' %in% colnames(allScore),
  "\n"
)

# compute normalised fTFC1/2 score ---------------------------------------------
## Add cancerNormal and normalised_score to allScore globally

allScore <- allScore %>%
  mutate(
    cancerNormal = case_when(
      cancerType %in% c('Normal') ~ 'Normal',
      cancerType == 'Normal.adj' ~ 'Normal.adj',
      TRUE ~ 'Tumour'
    )
  )

## Normalise within each dataset + module
## (tumour score - median normal score for that dataset/module)
allScore <- allScore %>%
  group_by(source, moduleType) %>%
  mutate(
    normalised_score = TotalScore -
      median(TotalScore[cancerNormal == 'Normal'], na.rm = TRUE)
  ) %>%
  ungroup()


# Score TDS in bulk ranked matrix ----------------------------------------------
if (is.null(PATHS$tdsscore_fp) || is.na(PATHS$tdsscore_fp)) {
  tds_scores <- singscore::simpleScore(
    rankData = bulk_ranked,
    upSet = tds_genes_ens
  )
  tds_scores$sampleID <- rownames(tds_scores)
  tds_scores <- tds_scores %>%
    dplyr::rename(
      TDS_TotalScore = TotalScore,
      TDS_TotalDispersion = TotalDispersion
    )
} else {
  tds_scores <- read.csv(PATHS$tdsscore_fp, row.names = 1)
}


# Gene overlap analysis --------------------------------------------------------
## Load your fTFC gene module
deg <- read.csv(file.path(outDir, 'fTFC1.2_top100_geneSignatures.csv'))
colnames(deg) <- c(
  'ensID',
  'geneSym',
  'chr',
  'logFC',
  'logCPM',
  'F',
  'PValue',
  'FDR',
  'pct_fTFC1',
  'pct_fTFC2',
  'direction',
  'module'
)

fTFC1_genes <- deg$geneSym[deg$direction == 'fTFC2_down']
fTFC2_genes <- deg$geneSym[deg$direction == 'fTFC2_up']

overlap_fTFC1 <- intersect(fTFC1_genes, tds_genes_symbols)
overlap_fTFC2 <- intersect(fTFC2_genes, tds_genes_symbols)

cat("\n========== Gene overlap with TDS ==========\n")
cat(sprintf(
  "fTFC1 (%d genes): %d overlap with TDS (%.1f%%)\n",
  length(fTFC1_genes),
  length(overlap_fTFC1),
  100 * length(overlap_fTFC1) / length(fTFC1_genes)
))
cat(
  "  Overlapping genes:",
  ifelse(
    length(overlap_fTFC1) == 0,
    "none",
    paste(overlap_fTFC1, collapse = ', ')
  ),
  "\n"
)

cat(sprintf(
  "fTFC2 (%d genes): %d overlap with TDS (%.1f%%)\n",
  length(fTFC2_genes),
  length(overlap_fTFC2),
  100 * length(overlap_fTFC2) / length(fTFC2_genes)
))
cat(
  "  Overlapping genes:",
  ifelse(
    length(overlap_fTFC2) == 0,
    "none",
    paste(overlap_fTFC2, collapse = ', ')
  ),
  "\n"
)

overlap_summary <- data.frame(
  signature = c('fTFC1', 'fTFC2'),
  n_total_genes = c(length(fTFC1_genes), length(fTFC2_genes)),
  n_TDS_overlap = c(length(overlap_fTFC1), length(overlap_fTFC2)),
  pct_overlap = round(
    c(
      100 * length(overlap_fTFC1) / length(fTFC1_genes),
      100 * length(overlap_fTFC2) / length(fTFC2_genes)
    ),
    1
  ),
  overlapping_genes = c(
    paste(overlap_fTFC1, collapse = ';'),
    paste(overlap_fTFC2, collapse = ';')
  )
)
write.csv(
  overlap_summary,
  file.path(outDir, 'TDS_fTFC_gene_overlap.csv'),
  row.names = FALSE
)


# Merge TDS scores with fTFC scores --------------------------------------------
## First check what the sample ID column is called in allScore
## Common options: 'Row.names', 'sampleID', 'Sample'
## Adjust 'sampleID_col' below accordingly

sampleID_col <- 'sampleID'

allScore_merged <- allScore %>%
  dplyr::filter(moduleType %in% c('fTFC1', 'fTFC2')) %>%
  dplyr::left_join(
    tds_scores %>% dplyr::select(sampleID, TDS_TotalScore),
    by = setNames('sampleID', sampleID_col)
  )

## Check merge worked
cat(sprintf(
  "\nMerge check: %d/%d rows have TDS score\n",
  sum(!is.na(allScore_merged$TDS_TotalScore)),
  nrow(allScore_merged)
))

## Normalise TDS within each dataset (same approach as fTFC)
allScore_merged <- allScore_merged %>%
  dplyr::group_by(source) %>%
  dplyr::mutate(
    TDS_normalised = TDS_TotalScore -
      median(TDS_TotalScore[cancerNormal == 'Normal'], na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

## Pivot wide: one row per sample, columns for fTFC1/fTFC2/TDS
allScore_wide <- allScore_merged %>%
  dplyr::select(
    all_of(sampleID_col),
    source,
    cancerType,
    cancerNormal,
    moduleType,
    TotalScore,
    normalised_score,
    TDS_TotalScore,
    TDS_normalised
  ) %>%
  tidyr::pivot_wider(
    names_from = moduleType,
    values_from = c(TotalScore, normalised_score),
    names_glue = "{moduleType}_{.value}"
  )

## Verify
cat("\nallScore_wide columns:\n")
print(colnames(allScore_wide))
cat("Rows:", nrow(allScore_wide), "\n")

write.csv(
  allScore_wide,
  file.path(outDir, 'TDS_fTFC_scores_wide.csv'),
  row.names = FALSE
)


# correlation analysis ---------------------------------------------------------
# conclusion: while there is a strong correlation between fTFC2_score and TDS_score,
#
## Paediatric dataset ----------------------------------------------------------
paed <- allScore_wide[
  allScore_wide$source %in% c('Sanger', 'REBC_THYR_paed', 'REBC_THYR_adult'),
]
paed$age_group = ifelse(paed$source == "REBC_THYR_adult", "adult", "paed")

ggplot(
  paed[paed$cancerNormal == 'Tumour', ],
  aes(TDS_normalised, fTFC2_normalised_score)
) +
  geom_point(aes(col = source)) +
  facet_wrap(vars(age_group)) +
  theme_classic()

## Focus on TCGA (largest N, most reliable) ------------------------------------

tcga <- allScore_wide %>%
  filter(source == 'TCGA_Thyroid', !is.na(TDS_TotalScore))

cat(sprintf(
  "\nTCGA samples for correlation: %d total (%d Normal, %d Tumour)\n",
  nrow(tcga),
  sum(tcga$cancerNormal == 'Normal'),
  sum(tcga$cancerNormal == 'Tumour')
))

## Spearman correlation: fTFC1/fTFC2 vs TDS
## All samples
cor_fTFC1_all <- cor.test(
  tcga$fTFC1_TotalScore,
  tcga$TDS_TotalScore,
  method = 'spearman'
)
cor_fTFC2_all <- cor.test(
  tcga$fTFC2_TotalScore,
  tcga$TDS_TotalScore,
  method = 'spearman'
)

## Tumour only
tcga_t <- tcga %>% filter(cancerNormal == 'Tumour')
tcga_n <- tcga %>% filter(cancerNormal == 'Normal')

cor_fTFC1_tumour <- cor.test(
  tcga_t$fTFC1_TotalScore,
  tcga_t$TDS_TotalScore,
  method = 'spearman'
)
cor_fTFC2_tumour <- cor.test(
  tcga_t$fTFC2_TotalScore,
  tcga_t$TDS_TotalScore,
  method = 'spearman'
)

cat("\n========== Spearman correlations with TDS (TCGA) ==========\n")
cat(sprintf(
  "fTFC1 vs TDS (all):    rho=%.3f, p=%.2e\n",
  cor_fTFC1_all$estimate,
  cor_fTFC1_all$p.value
))
cat(sprintf(
  "fTFC2 vs TDS (all):    rho=%.3f, p=%.2e\n",
  cor_fTFC2_all$estimate,
  cor_fTFC2_all$p.value
))
cat(sprintf(
  "fTFC1 vs TDS (tumour): rho=%.3f, p=%.2e\n",
  cor_fTFC1_tumour$estimate,
  cor_fTFC1_tumour$p.value
))
cat(sprintf(
  "fTFC2 vs TDS (tumour): rho=%.3f, p=%.2e\n",
  cor_fTFC2_tumour$estimate,
  cor_fTFC2_tumour$p.value
))

## Save correlation results
cor_results <- data.frame(
  comparison = c(
    'fTFC1_vs_TDS_all',
    'fTFC2_vs_TDS_all',
    'fTFC1_vs_TDS_tumour',
    'fTFC2_vs_TDS_tumour'
  ),
  rho = c(
    cor_fTFC1_all$estimate,
    cor_fTFC2_all$estimate,
    cor_fTFC1_tumour$estimate,
    cor_fTFC2_tumour$estimate
  ),
  p_value = c(
    cor_fTFC1_all$p.value,
    cor_fTFC2_all$p.value,
    cor_fTFC1_tumour$p.value,
    cor_fTFC2_tumour$p.value
  ),
  n = c(nrow(tcga), nrow(tcga), nrow(tcga_t), nrow(tcga_t))
)
write.csv(
  cor_results,
  file.path(outDir, 'TDS_fTFC_correlation_results.csv'),
  row.names = FALSE
)

# GLM modelling to see if fTFC2 score adds predictive power --------------------
# does fTFC2 add information beyond TDS?
response = allScore_wide$cancerNormal[
  allScore_wide$cancerNormal %in% c('Normal', 'Tumour')
]
## TCGA tumour + normal samples only
library(pROC)
glm_data <- allScore_wide %>%
  filter(
    source == 'TCGA_Thyroid',
    cancerNormal %in% c('Normal', 'Tumour'),
    !is.na(TDS_TotalScore),
    !is.na(fTFC2_TotalScore)
  ) %>%
  mutate(is_tumour = as.integer(cancerNormal == 'Tumour'))

## Model 1: TDS only
m1 <- glm(is_tumour ~ TDS_TotalScore, data = glm_data, family = binomial)

## Model 2: TDS + fTFC2
m2 <- glm(
  is_tumour ~ TDS_TotalScore + fTFC2_TotalScore,
  data = glm_data,
  family = binomial
)

## Likelihood ratio test
lrt <- anova(m1, m2, test = 'LRT')
cat("\n--- Likelihood ratio test: TDS vs TDS + fTFC2 ---\n")
print(lrt)

## AUC comparison
auc_m1 <- pROC::auc(glm_data$is_tumour, fitted(m1))
auc_m2 <- pROC::auc(glm_data$is_tumour, fitted(m2))
auc_test <- pROC::roc.test(
  pROC::roc(glm_data$is_tumour, fitted(m1)),
  pROC::roc(glm_data$is_tumour, fitted(m2))
)

cat(sprintf("\nAUC - TDS only:        %.3f\n", auc_m1))
cat(sprintf("AUC - TDS + fTFC2:     %.3f\n", auc_m2))
cat(sprintf("AUC improvement:       %.3f\n", auc_m2 - auc_m1))
cat(sprintf("DeLong test p-value:   %.2e\n", auc_test$p.value))

## Also test fTFC1 for comparison
m3 <- glm(
  is_tumour ~ TDS_TotalScore + fTFC1_TotalScore,
  data = glm_data,
  family = binomial
)

auc_m3 <- pROC::auc(glm_data$is_tumour, fitted(m3))
lrt_m3 <- anova(m1, m3, test = 'LRT')

cat(sprintf("\nAUC - TDS + fTFC1:     %.3f\n", auc_m3))
cat("\nLRT TDS vs TDS + fTFC1:\n")
print(lrt_m3)

# Partial correlation - fTFC2 after regressing out TDS -------------------------
## Fit regression on NORMAL samples: fTFC2 ~ TDS
## Then ask: do tumour samples deviate from this relationship?
## This is the cleanest test of whether fTFC2 carries information beyond TDS

## Check normal N is adequate
cat(sprintf("\nNormal N for regression: %d\n", nrow(tcga_n)))
if (nrow(tcga_n) < 20) {
  warning("Normal N is low - regression may be unstable")
}

lm_normal <- lm(fTFC2_TotalScore ~ TDS_TotalScore, data = tcga_n)
cat("\nRegression fit (fTFC2 ~ TDS in normal samples):\n")
print(summary(lm_normal))

## Predict expected fTFC2 for all TCGA samples based on normal regression
tcga$fTFC2_expected <- predict(lm_normal, newdata = tcga)
tcga$fTFC2_residual <- tcga$fTFC2_TotalScore - tcga$fTFC2_expected

## Test: do tumour residuals differ from normal residuals?
resid_test <- wilcox.test(
  tcga$fTFC2_residual[tcga$cancerNormal == 'Tumour'],
  tcga$fTFC2_residual[tcga$cancerNormal == 'Normal']
)

cat("\n========== Partial correlation test ==========\n")
cat(
  "Q: After accounting for TDS, does fTFC2 still differ between tumour and normal?\n"
)
cat(sprintf(
  "Wilcoxon W=%.0f, p=%.2e\n",
  resid_test$statistic,
  resid_test$p.value
))
cat(sprintf(
  "Median residual - Normal: %.4f\n",
  median(tcga$fTFC2_residual[tcga$cancerNormal == 'Normal'])
))
cat(sprintf(
  "Median residual - Tumour: %.4f\n",
  median(tcga$fTFC2_residual[tcga$cancerNormal == 'Tumour'])
))

if (resid_test$p.value < 0.05) {
  cat(">> fTFC2 carries information BEYOND TDS in TCGA\n")
} else {
  cat(">> fTFC2 residual NOT significantly different after TDS regression\n")
  cat("   Consider framing fTFC2 as substantially overlapping with TDS\n")
}

write.csv(
  tcga,
  file.path(outDir, 'TCGA_TDS_fTFC_residuals.csv'),
  row.names = FALSE
)


##--- Step 8: Plots ---##

## Plot A: Scatterplot fTFC2 vs TDS coloured by tumour/normal
plotFun_TDS_scatter <- function(noFrame = FALSE, noPlot = FALSE) {
  dd <- tcga %>%
    mutate(cancerNormal = factor(cancerNormal, c('Normal', 'Tumour')))

  ## Annotation text
  ann <- sprintf(
    "rho = %.2f\np = %.2e\n(all samples)",
    cor_fTFC2_all$estimate,
    cor_fTFC2_all$p.value
  )

  p <- ggplot(dd, aes(TDS_TotalScore, fTFC2_TotalScore, col = cancerNormal)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(method = 'lm', se = TRUE, linewidth = 0.8) +
    annotate(
      'text',
      x = -Inf,
      y = Inf,
      hjust = -0.1,
      vjust = 1.3,
      label = ann,
      size = 2.8
    ) +
    scale_color_manual(
      values = c('Normal' = grey(0.5), 'Tumour' = '#511378'),
      name = ''
    ) +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = F, colour = 'black'),
      axis.line = element_blank(),
      axis.text = element_text(colour = 'black'),
      legend.position = 'bottom'
    ) +
    xlab('TDS score') +
    ylab('fTFC2 score')

  print(p)
}
saveFig(
  file.path(plotDir, 'Fig_TDS_fTFC2_scatter'),
  plotFun_TDS_scatter,
  rawData = tcga,
  width = 3.5,
  height = 3.5,
  res = 500
)


## Plot B: fTFC2 residual after TDS regression
plotFun_TDS_residual <- function(noFrame = FALSE, noPlot = FALSE) {
  dd <- tcga %>%
    mutate(cancerNormal = factor(cancerNormal, c('Normal', 'Tumour')))

  plab <- ifelse(
    resid_test$p.value < 0.001,
    sprintf("p = %.2e", resid_test$p.value),
    sprintf("p = %.3f", resid_test$p.value)
  )

  p <- ggplot(dd, aes(cancerNormal, fTFC2_residual)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
    geom_quasirandom(
      size = 0.4,
      width = 0.15,
      alpha = 0.4,
      aes(col = cancerNormal)
    ) +
    geom_boxplot(
      aes(fill = cancerNormal),
      outlier.shape = NA,
      alpha = 0.8,
      width = 0.4,
      linewidth = 0.3,
      colour = 'black'
    ) +
    annotate(
      'text',
      x = 1.5,
      y = max(dd$fTFC2_residual, na.rm = T) * 0.95,
      label = plab,
      size = 3
    ) +
    scale_fill_manual(
      values = c('Normal' = grey(0.8), 'Tumour' = '#511378')
    ) +
    scale_color_manual(
      values = c('Normal' = grey(0.5), 'Tumour' = '#511378')
    ) +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = F, colour = 'black'),
      axis.line = element_blank(),
      axis.text = element_text(colour = 'black'),
      legend.position = 'none'
    ) +
    xlab('') +
    ylab('fTFC2 residual\n(after TDS regression)')

  print(p)
}
saveFig(
  file.path(plotDir, 'Fig_TDS_fTFC2_residual'),
  plotFun_TDS_residual,
  rawData = tcga,
  width = 2.2,
  height = 3.2,
  res = 500
)


## Plot C: Gene overlap bar chart
plotFun_overlap <- function(noFrame = FALSE, noPlot = FALSE) {
  dd <- data.frame(
    signature = rep(c('fTFC1', 'fTFC2'), each = 2),
    category = rep(c('TDS overlap', 'Signature-specific'), 2),
    n = c(
      length(overlap_fTFC1),
      length(fTFC1_genes) - length(overlap_fTFC1),
      length(overlap_fTFC2),
      length(fTFC2_genes) - length(overlap_fTFC2)
    )
  ) %>%
    mutate(category = factor(category, c('TDS overlap', 'Signature-specific')))

  p <- ggplot(dd, aes(signature, n, fill = category)) +
    geom_bar(stat = 'identity', colour = 'black', linewidth = 0.3) +
    scale_fill_manual(
      values = c('TDS overlap' = '#d4a017', 'Signature-specific' = '#511378'),
      name = ''
    ) +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = F, colour = 'black'),
      axis.line = element_blank(),
      axis.text = element_text(colour = 'black'),
      legend.position = 'bottom'
    ) +
    xlab('') +
    ylab('Number of genes')

  print(p)
}
saveFig(
  file.path(plotDir, 'Fig_TDS_gene_overlap'),
  plotFun_overlap,
  rawData = overlap_summary,
  width = 2.5,
  height = 3,
  res = 500
)


##--- Step 9: Summary printout ---##
cat("\n\n========== TDS ANALYSIS SUMMARY ==========\n")
cat("\nGene overlap:\n")
print(overlap_summary[, c(
  'signature',
  'n_total_genes',
  'n_TDS_overlap',
  'pct_overlap'
)])
cat("\nSpearman correlations with TDS (TCGA):\n")
print(cor_results)
cat("\nPartial correlation (fTFC2 residual after TDS):\n")
cat(sprintf(
  "  p = %.2e  -->  %s\n",
  resid_test$p.value,
  ifelse(
    resid_test$p.value < 0.05,
    "fTFC2 carries signal BEYOND TDS",
    "fTFC2 not significantly different from TDS"
  )
))
cat("===========================================\n")

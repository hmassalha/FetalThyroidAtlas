# Correlate fTFC1/2 score with TCGA clinical output

# Libraries --------------------------------------------------------------------
library(tidyverse)

# Set up -----------------------------------------------------------------------
setwd("~/FetalThyroidAtlas/")

outDir <- "Results/2505/PTC_bulkRNAseq"
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = T)
}

plotDir <- "Figures/2505"

# Helpers ----------------------------------------------------------------------
library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(checkmate)
auto_numeric_convert <- function(x, threshold = 0.9) {
  ## preserve original class
  if (is.numeric(x)) {
    return(x)
  }

  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (!is.character(x)) {
    return(x)
  }

  ## decimal comma fix
  x_clean <- gsub(",", ".", x[!is.na(x) & x != ""], fixed = TRUE)

  ## attempt conversion
  x_num <- suppressWarnings(as.numeric(x_clean))

  ## proportion successfully converted
  prop_numeric <- mean(
    is.na(x_clean) | !is.na(x_num)
  )

  ## convert only if mostly numeric
  if (prop_numeric >= threshold) {
    x_num_full <- as.numeric(gsub(",", ".", x, fixed = TRUE))
    return(x_num_full)
  }

  return(x)
}
auto_detect_plot_type <- function(
  x,
  max_levels_numeric = 5,
  numeric_threshold = 0.9
) {
  x <- auto_numeric_convert(x)

  ## character / factor → categorical
  if (is.character(x) || is.factor(x)) {
    return(list(
      type = "categorical",
      data = x
    ))
  }

  ## numeric handling
  if (is.numeric(x)) {
    n_unique <- data.table::uniqueN(
      stats::na.omit(x)
    )

    ## binary / small discrete numeric
    if (n_unique <= max_levels_numeric) {
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
  module_keep = c("fTFC1", "fTFC2"),
  outfile = "TCGA_clinical_associations.pdf",
  print_plot = TRUE,
  min_group_n_for_test = 3
) {
  checkmate::assert_data_frame(data)
  checkmate::assert_character(vars, min.len = 1)
  checkmate::assert_string(score_col)
  checkmate::assert_flag(print_plot)
  checkmate::assert_int(min_group_n_for_test, lower = 1)

  dt <- data.table::as.data.table(data)

  checkmate::assert_subset(
    c(score_col, "moduleType"),
    colnames(dt)
  )

  vars <- vars[vars %in% colnames(dt)]

  dt <- dt[moduleType %in% module_keep]
  if (!is.null(outfile)) {
    pdf(
      outfile,
      width = 12,
      height = 8,
      onefile = TRUE
    )
    on.exit(dev.off())
  }

  for (vv in vars) {
    message("Processing: ", vv)

    dd <- copy(
      dt[
        !is.na(get(vv)) &
          !is.na(get(score_col))
      ]
    )

    if (nrow(dd) < 10) {
      next
    }
    det <- auto_detect_plot_type(x = dd[[vv]])
    plot_type <- det$type
    x <- det$data

    try(
      {
        ## -------------------------
        ## categorical variables
        ## -------------------------

        if (!is.numeric(x)) {
          dd$plot_var <- dd[[vv]]

          if (is.factor(dd$plot_var)) {
            dd$plot_var_display <- droplevels(dd$plot_var)
          } else {
            dd$plot_var_display <- reorder(
              dd$plot_var,
              dd[[score_col]],
              median,
              na.rm = TRUE
            )
          }

          if (nlevels(factor(dd$plot_var_display)) < 2) {
            next
          }

          level_order <- levels(factor(dd$plot_var_display))
          level_counts <- table(factor(
            dd$plot_var_display,
            levels = level_order
          ))
          level_labels <- paste0(
            level_order,
            "\n(n=",
            as.integer(level_counts),
            ")"
          )
          dd$plot_var_display_n <- factor(
            dd$plot_var_display,
            levels = level_order,
            labels = level_labels
          )

          panel_split <- split(
            dd,
            interaction(dd$moduleType, dd$cancerNormal, drop = TRUE)
          )
          panel_stats <- lapply(panel_split, function(subd) {
            grp <- droplevels(factor(subd$plot_var_display))
            y <- subd[[score_col]]

            grp_counts <- table(grp)
            keep_levels <- names(grp_counts)[grp_counts >= min_group_n_for_test]
            keep_idx <- grp %in% keep_levels

            grp_test <- droplevels(grp[keep_idx])
            y_test <- y[keep_idx]
            n_groups <- nlevels(grp_test)

            if (n_groups < 2) {
              return(data.frame(
                moduleType = as.character(subd$moduleType[1]),
                cancerNormal = as.character(subd$cancerNormal[1]),
                test_name = "Insufficient n",
                p_raw = NA_real_,
                posthoc_str = "",
                use_for_bh = FALSE,
                stringsAsFactors = FALSE
              ))
            }

            if (n_groups == 2) {
              grp_levels <- levels(grp_test)
              x1 <- y_test[grp_test == grp_levels[1]]
              x2 <- y_test[grp_test == grp_levels[2]]
              test_res <- tryCatch(
                suppressWarnings(stats::wilcox.test(x = x1, y = x2)),
                error = function(e) NULL
              )
              test_name <- "Wilcoxon"
              posthoc_str <- ""
            } else {
              test_res <- tryCatch(
                suppressWarnings(stats::kruskal.test(x = y_test, g = grp_test)),
                error = function(e) NULL
              )
              test_name <- "Kruskal"

              ## pairwise post-hoc (run now; display conditioned on q later)
              posthoc_str <- tryCatch(
                {
                  pw <- suppressWarnings(
                    stats::pairwise.wilcox.test(
                      x = y_test,
                      g = grp_test,
                      p.adjust.method = "BH"
                    )
                  )
                  pm <- pw$p.value
                  idx <- which(!is.na(pm), arr.ind = TRUE)
                  pairs <- apply(idx, 1, function(i) {
                    rn <- rownames(pm)[i[1]]
                    cn <- colnames(pm)[i[2]]
                    pv <- signif(pm[i[1], i[2]], 2)
                    paste0(cn, " vs ", rn, ": q=", pv)
                  })
                  paste(pairs, collapse = "\n")
                },
                error = function(e) ""
              )
            }

            if (is.null(test_res)) {
              return(NULL)
            }

            data.frame(
              moduleType = as.character(subd$moduleType[1]),
              cancerNormal = as.character(subd$cancerNormal[1]),
              test_name = test_name,
              p_raw = test_res$p.value,
              posthoc_str = if (exists("posthoc_str")) posthoc_str else "",
              use_for_bh = TRUE,
              stringsAsFactors = FALSE
            )
          })
          panel_stats <- Filter(Negate(is.null), panel_stats)

          if (length(panel_stats) > 0) {
            panel_stats <- do.call(rbind, panel_stats)
            panel_stats$q_bh <- NA_real_
            test_idx <- panel_stats$use_for_bh
            panel_stats$q_bh[test_idx] <- p.adjust(
              panel_stats$p_raw[test_idx],
              method = "BH"
            )
            panel_stats$label <- ifelse(
              panel_stats$use_for_bh,
              paste0(
                panel_stats$test_name,
                "\nq=",
                signif(panel_stats$q_bh, 3),
                ifelse(
                  panel_stats$test_name == "Kruskal" &
                    !is.na(panel_stats$q_bh) &
                    panel_stats$q_bh < 0 &
                    nchar(panel_stats$posthoc_str) > 0,
                  paste0("\n", panel_stats$posthoc_str),
                  ""
                )
              ),
              paste0("Insufficient n\n(<", min_group_n_for_test, "/group)")
            )
          } else {
            panel_stats <- NULL
          }

          p <- ggplot(
            dd,
            aes(
              .data$plot_var_display_n,
              .data[[score_col]],
              fill = .data$moduleType
            )
          ) +
            geom_hline(yintercept = 0, linewidth = 0.5) +
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
            scale_fill_manual(
              values = c(
                "fTFC1" = grey(0.6),
                "fTFC2" = "orange",
                "TDS" = "#4D9DE0"
              )
            ) +
            facet_grid(
              ~ moduleType + cancerNormal,
              scales = "free_x",
              space = "free"
            ) +
            {
              if (!is.null(panel_stats)) {
                geom_text(
                  data = panel_stats,
                  aes(x = Inf, y = Inf, label = .data$label),
                  inherit.aes = FALSE,
                  hjust = 1.05,
                  vjust = 1.2,
                  size = 3.1
                )
              }
            } +
            labs(
              title = vv,
              x = NULL,
              y = "Enrichment score"
            ) +
            theme_classic(base_size = 12) +
            theme(
              axis.line = element_blank(),
              panel.border = element_rect(
                fill = NA,
                colour = "black",
                linewidth = 0.7
              ),
              strip.background = element_blank(),
              axis.text.x = element_text(
                angle = 90,
                hjust = 1,
                vjust = 0.5,
                size = 8
              ),
              axis.ticks = element_line(linewidth = 0.5),
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
            geom_hline(yintercept = 0, linewidth = 0.5) +
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
              axis.line = element_blank(),
              panel.border = element_rect(
                fill = NA,
                colour = "black",
                linewidth = 0.7
              ),
              axis.ticks = element_line(linewidth = 0.5),
              strip.background = element_blank()
            )
        }

        if (isTRUE(print_plot)) {
          print(p)
        }
      },
      silent = TRUE
    )
  }
  return(p)
}

# Import datasets --------------------------------------------------------------
## fTFC1/fTFC2 signature in bulk RNAseq data -----------------------------------
allScore <- data.table::fread(file.path(
  outDir,
  "fTFC1.2_top100_geneSignatures_SingScore_bulkRNAseq.csv"
))
## Add cancerNormal and normalised_score to allScore globally
allScore <- allScore |>
  mutate(
    cancerNormal = case_when(
      cancerType %in% c("Normal") ~ "Normal",
      cancerType == "Normal.adj" ~ "Normal.adj",
      TRUE ~ "Tumour"
    )
  )

## Normalise within each dataset + module
## (tumour score - median normal score for that dataset/module)
allScore <- allScore |>
  group_by(source, moduleType) |>
  mutate(
    normalised_score = TotalScore -
      median(TotalScore[cancerNormal == "Normal"], na.rm = TRUE)
  ) |>
  ungroup()

## TCGA metadata ---------------------------------------------------------------
tcga_se_path <- "Data/published_bulkRNAseq/TCGA_Thyroid/TCGA_Thyroid_bulkRNA_se.RDS"
tcga_se <- readRDS(tcga_se_path)
tcga_mdat <- as.data.frame(colData(tcga_se))
tcga_mdat$sampleID <- rownames(tcga_mdat)
tcga_mdat$source <- "TCGA_Thyroid"
tcga_mdat$sampleName <- rownames(tcga_mdat)
tcga_mdat$cancerType <- ifelse(
  tcga_mdat$tissue_type == "Normal",
  "Normal",
  paste0("PTC_", tcga_mdat$classification_of_tumor)
)
tcga_mdat$age <- tcga_mdat$age_at_diagnosis
tcga_mdat$sex <- tcga_mdat$gender


table(tcga_mdat$sampleID %in% allScore$sampleID)


# Correlation analysis ---------------------------------------------------------
samples_to_keep = intersect(allScore$sampleID, tcga_mdat$sampleID)
tcga_data = merge(
  tcga_mdat[match(samples_to_keep, tcga_mdat$sampleID), ],
  allScore_merged[
    allScore_merged$sampleID %in% samples_to_keep,
    c(
      "sampleID",
      colnames(allScore_merged)[
        !colnames(allScore_merged) %in% colnames(tcga_mdat)
      ]
    )
  ],
  by = 'sampleID',
  all = TRUE
)

checkmate::assert_true(
  nrow(tcga_data) == 526 * n_distinct(allScore_merged$moduleType)
)
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
# # Tumour sub-types
# "definition", "tumor_descriptor","sample_type",
# "ajcc_pathologic_stage","synchronous_malignancy","ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m"
# # survival metrics
# "vital_status", "days_to_diagnosis", "days_to_last_follow_up"
# "age_at_diagnosis",
# "year_of_diagnosis","days_to_death"
# # clinical treatment
# "treatments","primary_diagnosis", "prior_malignancy","prior_treatment","diagnosis_is_primary_disease",
# "residual_disease"
# "classification_of_tumor", "tumor_focality"
# # patient metadata
# "sex", "age"

ggplot(
  tcga_data[tcga_data$moduleType %in% c("fTFC1", "fTFC2"), ],
  aes(reorder(paper_BRAF, TotalScore, median), TotalScore)
) +
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


vars_of_interest <- c(
  "definition",
  "tumor_descriptor",
  "sample_type",
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
  score_col = "normalised_score",
  outfile = file.path(
    outDir,
    "TCGA_reviewer4_clinical_associations_normScore.pdf"
  )
)

list_cols <- names(tcga_data)[vapply(tcga_data, is.list, logical(1))]
data.table::fwrite(
  tcga_data %>% dplyr::select(-all_of(list_cols)),
  file = file.path(outDir, 'TCGA_singscore_signatures.csv')
)

# Refined plots ----------------------------------------------------------------
tcga_data <- read.csv("../tmp_thyroid/TCGA_singscore_signatures.csv")
clinical_outcome_var <- c(
  "ajcc_pathologic_n",
  "paper_pathologic_N",
  "paper_N0vsN1b",
  "paper_Follow_up_New_Tumor_Event",
  "paper_Extrathyroidal_extension",
  "vital_status",
  "days_to_death",
  "days_to_last_follow_up",
  "ajcc_pathologic_t",
  "paper_pathologic_T",
  "ajcc_pathologic_stage",
  "paper_Neoplasm_Disease_Stage"
)
sapply(
  tcga_mdat[tcga_mdat$cancerType != "Normal", ][, clinical_outcome_var],
  function(x) sum(is.na(x))
)
lapply(tcga_mdat[, clinical_outcome_var], table)

plot_tcga_associations(
  data = tcga_data,
  vars = clinical_outcome_var,
  score_col = "TDS_normalised",
  outfile = file.path(
    outDir,
    "TCGA_reviewer4_clinical_associations_selectedVar_TDS.pdf"
  )
)

# Reviewer comment figure (single page layout) --------------------------------
reviewer_plots <- list()

# A: Nodal metastasis status (AJCC N)
d <- tcga_data[
  tcga_data$moduleType %in%
    c("fTFC1", "fTFC2") &
    tcga_data$cancerNormal != "Normal" &
    tcga_data$ajcc_pathologic_n != "",
]
d$ajcc_pathologic_n <- factor(
  d$ajcc_pathologic_n,
  c("N0", "N1", "N1a", "N1b")
)
reviewer_plots[["A"]] <- plot_tcga_associations(
  data = d,
  vars = "ajcc_pathologic_n",
  score_col = "normalised_score",
  outfile = NULL,
  print_plot = FALSE
) +
  labs(
    title = "Nodal Metastasis Status (AJCC N)",
    y = "Centralised enrichment score"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    title = element_text(size = 10.5)
  )

# B: Extrathyroidal extension category
d <- tcga_data[
  tcga_data$moduleType %in%
    c("fTFC1", "fTFC2") &
    tcga_data$cancerNormal != "Normal" &
    !tcga_data$paper_Extrathyroidal_extension %in%
      c("", "[Not Available]", "[Unknown]"),
]
d$paper_Extrathyroidal_extension <- factor(
  d$paper_Extrathyroidal_extension,
  c("None", "Minimal (T3)", "Moderate/Advanced (T4a)", "Very Advanced (T4b)")
)
reviewer_plots[["B"]] <- plot_tcga_associations(
  data = d,
  vars = "paper_Extrathyroidal_extension",
  score_col = "normalised_score",
  outfile = NULL,
  print_plot = FALSE
) +
  labs(
    title = "Extrathyroidal Extension Severity",
    y = "Centralised enrichment score"
  ) +
  theme(title = element_text(size = 10.5))

# C: Disease stage
d <- tcga_data[
  tcga_data$moduleType %in%
    c("fTFC1", "fTFC2") &
    tcga_data$cancerNormal != "Normal" &
    tcga_data$paper_Neoplasm_Disease_Stage != "",
]
d$paper_Neoplasm_Disease_Stage <- factor(
  d$paper_Neoplasm_Disease_Stage,
  c("Stage I", "Stage II", "Stage III", "Stage IV", "Stage IVA", "Stage IVC")
)
reviewer_plots[["C"]] <- plot_tcga_associations(
  data = d,
  vars = "paper_Neoplasm_Disease_Stage",
  score_col = "normalised_score",
  outfile = NULL,
  print_plot = FALSE
) +
  labs(title = "Neoplasm Disease Stage", y = "Centralised enrichment score") +
  theme(title = element_text(size = 10.5))

# D: Patient vital status
d <- tcga_data[
  tcga_data$moduleType %in%
    c("fTFC1", "fTFC2") &
    tcga_data$cancerNormal != "Normal",
]
d$vital_status <- factor(
  d$vital_status,
  c("Alive", "Dead")
)
reviewer_plots[["D"]] <- plot_tcga_associations(
  data = d,
  vars = "vital_status",
  score_col = "normalised_score",
  outfile = NULL,
  print_plot = FALSE
) +
  labs(title = "Vital Status", y = "Centralised enrichment score") +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    title = element_text(size = 10.5)
  )

# Add driver annotation
driver_annot <- data.table::fread(
  "../tmp_thyroid/SupFig15_fTFC1.2_moduleScore_bulkSamples_driver_rawData.tsv"
)
driver_annot <- driver_annot[moduleType == "fTFC2"][, c(
  "sampleID",
  "driver",
  "driver_category",
  "group"
)]
table(tcga_data$sampleID %in% driver_annot$sampleID)
table(driver_annot$sampleID %in% tcga_data$sampleID)
dim(tcga_data)
tcga_data <- merge(tcga_data, driver_annot, by = "sampleID", all.x = TRUE)
dim(tcga_data)
tcga_data$group[
  tcga_data$group == "unknown" & tcga_data$paper_RAS == 1
] <- "RAS"

p <- reviewer_plots[["C"]] +
  facet_grid(
    ~ moduleType + cancerNormal + group,
    scales = "free_x",
    space = "free"
  )
p

# E: Molecular driver group
d <- tcga_data[
  tcga_data$moduleType %in%
    c("fTFC1", "fTFC2") &
    tcga_data$cancerNormal != "Normal",
]
d$group <- factor(
  d$group,
  c(
    'BRAF',
    'NCOA4_RET',
    'CCDC6_RET',
    'RET-OTHER',
    'NTRK',
    'RAS',
    'others',
    'unknown'
  )
)
reviewer_plots[["E"]] <- plot_tcga_associations(
  data = d,
  vars = "group",
  score_col = "normalised_score",
  outfile = NULL,
  print_plot = FALSE
) +
  labs(title = "Molecular Driver Group", y = "Centralised enrichment score") +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    title = element_text(size = 10.5)
  )

d$group[!is.na(d$paper_RAS) & d$paper_RAS == 1] <- "RAS"
table(d$group, d$ajcc_pathologic_n)


# extra plot: Disease stage by driver group
d <- tcga_data[
  tcga_data$moduleType %in%
    c("fTFC1", "fTFC2") &
    tcga_data$cancerNormal != "Normal" &
    !tcga_data$paper_Neoplasm_Disease_Stage %in% c("", '[Not Available]'),
]
# d$paper_Neoplasm_Disease_Stage <- factor(
#   d$paper_Neoplasm_Disease_Stage,
#   c("Stage I", "Stage II", "Stage III", "Stage IV", "Stage IVA", "Stage IVC")
# )
d$group2 = paste0(as.character(d$paper_Neoplasm_Disease_Stage), '_', d$group)
p <- plot_tcga_associations(
  data = d,
  vars = "group2",
  score_col = "normalised_score",
  outfile = NULL,
  print_plot = FALSE
) +
  labs(title = "Neoplasm Disease Stage", y = "Centralised enrichment score") +
  theme(title = element_text(size = 10.5)) +
  facet_grid(
    moduleType ~ cancerNormal + group,
    scales = "free_x",
    space = "free"
  )
p


min_group_n_for_test = 3
vars = "group2"
module_keep <- c('fTFC1', 'fTFC2')
score_col <- "normalised_score"
dt <- data.table::as.data.table(d)
dt <- dt[moduleType %in% module_keep]
dd <- copy(
  dt[
    !is.na(get(vars)) &
      !is.na(get(score_col))
  ]
)

if (nrow(dd) < 10) {
  next
}
det <- auto_detect_plot_type(x = dd[[vars]])
plot_type <- det$type
x <- det$data
dd$plot_var <- dd[[vars]]

if (is.factor(dd$plot_var)) {
  dd$plot_var_display <- droplevels(dd$plot_var)
} else {
  dd$plot_var_display <- reorder(
    dd$plot_var,
    dd[[score_col]],
    median,
    na.rm = TRUE
  )
}

if (nlevels(factor(dd$plot_var_display)) < 2) {
  next
}

level_order <- levels(factor(dd$plot_var_display))
level_counts <- table(factor(
  dd$plot_var_display,
  levels = level_order
))
level_labels <- paste0(
  level_order,
  "\n(n=",
  as.integer(level_counts),
  ")"
)
dd$plot_var_display_n <- factor(
  dd$plot_var_display,
  levels = level_order,
  labels = level_labels
)
panel_split <- split(
  dd,
  interaction(dd$moduleType, dd$group, drop = TRUE)
)
panel_stats <- lapply(panel_split, function(subd) {
  grp <- droplevels(factor(subd$plot_var_display))
  y <- subd[[score_col]]

  grp_counts <- table(grp)
  keep_levels <- names(grp_counts)[grp_counts >= min_group_n_for_test]
  keep_idx <- grp %in% keep_levels

  grp_test <- droplevels(grp[keep_idx])
  y_test <- y[keep_idx]
  n_groups <- nlevels(grp_test)

  if (n_groups < 2) {
    return(data.frame(
      moduleType = as.character(subd$moduleType[1]),
      cancerNormal = as.character(subd$cancerNormal[1]),
      group = as.character(subd$group[1]),
      test_name = "Insufficient n",
      p_raw = NA_real_,
      posthoc_str = "",
      use_for_bh = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  if (n_groups == 2) {
    grp_levels <- levels(grp_test)
    x1 <- y_test[grp_test == grp_levels[1]]
    x2 <- y_test[grp_test == grp_levels[2]]
    test_res <- tryCatch(
      suppressWarnings(stats::wilcox.test(x = x1, y = x2)),
      error = function(e) NULL
    )
    test_name <- "Wilcoxon"
    posthoc_str <- ""
  } else {
    test_res <- tryCatch(
      suppressWarnings(stats::kruskal.test(x = y_test, g = grp_test)),
      error = function(e) NULL
    )
    test_name <- "Kruskal"

    ## pairwise post-hoc (run now; display conditioned on q later)
    posthoc_str <- tryCatch(
      {
        pw <- suppressWarnings(
          stats::pairwise.wilcox.test(
            x = y_test,
            g = grp_test,
            p.adjust.method = "BH"
          )
        )
        pm <- pw$p.value
        idx <- which(!is.na(pm), arr.ind = TRUE)
        pairs <- apply(idx, 1, function(i) {
          rn <- rownames(pm)[i[1]]
          cn <- colnames(pm)[i[2]]
          pv <- signif(pm[i[1], i[2]], 2)
          paste0(cn, " vs ", rn, ": q=", pv)
        })
        paste(pairs, collapse = "\n")
      },
      error = function(e) ""
    )
  }

  if (is.null(test_res)) {
    return(NULL)
  }

  data.frame(
    moduleType = as.character(subd$moduleType[1]),
    cancerNormal = as.character(subd$cancerNormal[1]),
    group = as.character(subd$group[1]),
    test_name = test_name,
    p_raw = test_res$p.value,
    posthoc_str = if (exists("posthoc_str")) posthoc_str else "",
    use_for_bh = TRUE,
    stringsAsFactors = FALSE
  )
})
panel_stats <- Filter(Negate(is.null), panel_stats)

if (length(panel_stats) > 0) {
  panel_stats <- do.call(rbind, panel_stats)
  panel_stats$q_bh <- NA_real_
  test_idx <- panel_stats$use_for_bh
  panel_stats$q_bh[test_idx] <- p.adjust(
    panel_stats$p_raw[test_idx],
    method = "BH"
  )
  panel_stats$label <- ifelse(
    panel_stats$use_for_bh,
    paste0(
      panel_stats$test_name,
      "\nq=",
      signif(panel_stats$q_bh, 3),
      ifelse(
        panel_stats$test_name == "Kruskal" &
          !is.na(panel_stats$q_bh) &
          panel_stats$q_bh < 0 &
          nchar(panel_stats$posthoc_str) > 0,
        paste0("\n", panel_stats$posthoc_str),
        ""
      )
    ),
    paste0("Insufficient n\n(<", min_group_n_for_test, "/group)")
  )
} else {
  panel_stats <- NULL
}
dd$group = factor(
  dd$group,
  c(
    'BRAF',
    'NCOA4_RET',
    'CCDC6_RET',
    'RET-OTHER',
    'NTRK',
    'RAS',
    'others',
    'unknown'
  )
)
panel_stats$group = factor(
  panel_stats$group,
  c(
    'BRAF',
    'NCOA4_RET',
    'CCDC6_RET',
    'RET-OTHER',
    'NTRK',
    'RAS',
    'others',
    'unknown'
  )
)
dd$plot_var_display_n <- sub(
  paste(
    paste0('_', unique(as.character(dd$group))),
    collapse = "|"
  ),
  '',
  dd$plot_var_display_n
)
p <- ggplot(
  dd,
  aes(
    .data$plot_var_display_n,
    .data[[score_col]],
    fill = .data$moduleType
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
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
  scale_fill_manual(
    values = c(
      "fTFC1" = grey(0.6),
      "fTFC2" = "orange",
      "TDS" = "#4D9DE0"
    )
  ) +
  facet_grid(
    moduleType ~ group,
    scales = "free_x",
    space = "free"
  ) +
  {
    if (!is.null(panel_stats)) {
      geom_text(
        data = panel_stats,
        aes(x = Inf, y = Inf, label = .data$label),
        inherit.aes = FALSE,
        hjust = 1.05,
        vjust = 1.2,
        size = 3.1
      )
    }
  } +
  labs(
    title = vars,
    x = NULL,
    y = "Enrichment score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.7
    ),
    strip.background = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 8
    ),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "none"
  ) +
  labs(title = "Neoplasm Disease Stage", y = "Centralised enrichment score") +
  theme(title = element_text(size = 10.5))
p
reviewer_plots[["F"]] <- p

ggsave(
  filename = file.path(
    plotDir,
    "TCGA_reviewer4_figure_diseaseStage_Driver.pdf"
  ),
  plot = p,
  width = 10,
  height = 4.7
)


reviewer_comment_figure <- patchwork::wrap_plots(
  reviewer_plots,
  design = "AB\nCD\nEE\nFF",
  heights = c(1, 1, 1.2, 2),
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "fTFC Signatures and TCGA Clinical Outcome",
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 11, face = "bold"))
  )

print(reviewer_comment_figure)

ggsave(
  filename = file.path(
    plotDir,
    "TCGA_reviewer4_figure_selected_clinical_outcomes.pdf"
  ),
  plot = reviewer_comment_figure,
  width = 11,
  height = 18
)


# fTFC1 / fTFC2 / TDS scores by molecular driver group ------------------------
## Reshape to long format: fTFC1 and fTFC2 come from normalised_score;
## TDS comes from TDS_normalised (one value per sample, so deduplicate via
## the fTFC1 rows before binding).
d_ftfc <- tcga_data |>
  dplyr::filter(
    moduleType %in% c("fTFC1", "fTFC2"),
    cancerNormal != "Normal",
    !is.na(group),
    group != ""
  ) |>
  dplyr::select(sampleID, module = moduleType, score = normalised_score, group)

d_tds <- tcga_data |>
  dplyr::filter(
    moduleType == "fTFC1", # one row per sample
    cancerNormal != "Normal",
    !is.na(group),
    group != ""
  ) |>
  dplyr::transmute(
    sampleID,
    module = "TDS",
    score = TDS_normalised,
    group
  )

df_modules <- dplyr::bind_rows(d_ftfc, d_tds) |>
  dplyr::mutate(module = factor(module, c("fTFC1", "fTFC2", "TDS")))

## per-module n labels for the x-axis
module_n <- df_modules |>
  dplyr::count(module) |>
  dplyr::mutate(module_label = paste0(module, "\n(n=", n, ")"))
df_modules <- dplyr::left_join(
  df_modules,
  module_n[, c("module", "module_label")],
  by = "module"
)
df_modules$module_label <- factor(
  df_modules$module_label,
  levels = module_n$module_label
)

module_colors <- c(
  "fTFC1" = grey(0.6),
  "fTFC2" = "orange",
  "TDS" = "#4D9DE0"
)

p_drivers <- ggplot(
  df_modules,
  aes(x = module_label, y = score, fill = module)
) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  ggbeeswarm::geom_quasirandom(width = 0.15, alpha = 0.45, size = 0.5) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.75,
    colour = "black",
    linewidth = 0.35
  ) +
  scale_fill_manual(values = module_colors) +
  facet_grid(~group, scales = "free_x", space = "free") +
  labs(
    title = "fTFC and TDS Enrichment Scores by Molecular Driver",
    x = NULL,
    y = "Normalised enrichment score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.7),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 9),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 9),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "none"
  )

print(p_drivers)

ggsave(
  filename = file.path(plotDir, "fTFC_TDS_scores_by_driver_group.pdf"),
  plot = p_drivers,
  width = 8,
  height = 5
)

# Reviewer figure v2: fTFC1 + fTFC2 + TDS as three module types --------------
## Helper: duplicate each filtered dataset with a synthetic TDS moduleType row
add_tds_rows <- function(d) {
  tds_rows <- d[d$moduleType == "fTFC1", ]
  tds_rows$moduleType <- "TDS"
  tds_rows$normalised_score <- tds_rows$TDS_normalised
  rbind(d, tds_rows)
}

reviewer_plots_tds <- list()

# A
d <- tcga_data[
  tcga_data$moduleType %in%
    c("fTFC1", "fTFC2") &
    tcga_data$cancerNormal != "Normal" &
    tcga_data$ajcc_pathologic_n != "",
]
d$ajcc_pathologic_n <- factor(d$ajcc_pathologic_n, c("N0", "N1", "N1a", "N1b"))
reviewer_plots_tds[["A"]] <- plot_tcga_associations(
  data = add_tds_rows(d),
  vars = "ajcc_pathologic_n",
  score_col = "normalised_score",
  module_keep = c("fTFC1", "fTFC2", "TDS"),
  outfile = NULL,
  print_plot = FALSE
) +
  labs(
    title = "Nodal Metastasis Status (AJCC N)",
    y = "Centralised enrichment score"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    title = element_text(size = 10.5)
  )

# B
d <- tcga_data[
  tcga_data$moduleType %in%
    c("fTFC1", "fTFC2") &
    tcga_data$cancerNormal != "Normal" &
    !tcga_data$paper_Extrathyroidal_extension %in%
      c("", "[Not Available]", "[Unknown]"),
]
d$paper_Extrathyroidal_extension <- factor(
  d$paper_Extrathyroidal_extension,
  c("None", "Minimal (T3)", "Moderate/Advanced (T4a)", "Very Advanced (T4b)")
)
reviewer_plots_tds[["B"]] <- plot_tcga_associations(
  data = add_tds_rows(d),
  vars = "paper_Extrathyroidal_extension",
  score_col = "normalised_score",
  module_keep = c("fTFC1", "fTFC2", "TDS"),
  outfile = NULL,
  print_plot = FALSE
) +
  labs(
    title = "Extrathyroidal Extension Severity",
    y = "Centralised enrichment score"
  ) +
  theme(title = element_text(size = 10.5))

# C
d <- tcga_data[
  tcga_data$moduleType %in%
    c("fTFC1", "fTFC2") &
    tcga_data$cancerNormal != "Normal",
]
d$vital_status <- factor(d$vital_status, c("Alive", "Dead"))
reviewer_plots_tds[["C"]] <- plot_tcga_associations(
  data = add_tds_rows(d),
  vars = "vital_status",
  score_col = "normalised_score",
  module_keep = c("fTFC1", "fTFC2", "TDS"),
  outfile = NULL,
  print_plot = FALSE
) +
  labs(title = "Vital Status", y = "Centralised enrichment score") +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    title = element_text(size = 10.5)
  )

# D
d <- tcga_data[
  tcga_data$moduleType %in%
    c("fTFC1", "fTFC2") &
    tcga_data$cancerNormal != "Normal" &
    tcga_data$paper_Neoplasm_Disease_Stage != "",
]
d$paper_Neoplasm_Disease_Stage <- factor(
  d$paper_Neoplasm_Disease_Stage,
  c("Stage I", "Stage II", "Stage III", "Stage IV", "Stage IVA", "Stage IVC")
)
reviewer_plots_tds[["C"]] <- plot_tcga_associations(
  data = add_tds_rows(d),
  vars = "paper_Neoplasm_Disease_Stage",
  score_col = "normalised_score",
  module_keep = c("fTFC1", "fTFC2", "TDS"),
  outfile = NULL,
  print_plot = FALSE
) +
  labs(title = "Neoplasm Disease Stage", y = "Centralised enrichment score") +
  theme(title = element_text(size = 10.5))

# E
d <- tcga_data[
  tcga_data$moduleType %in%
    c("fTFC1", "fTFC2") &
    tcga_data$cancerNormal != "Normal",
]
reviewer_plots_tds[["E"]] <- plot_tcga_associations(
  data = add_tds_rows(d),
  vars = "group",
  score_col = "normalised_score",
  module_keep = c("fTFC1", "fTFC2", "TDS"),
  outfile = NULL,
  print_plot = FALSE
) +
  labs(title = "Molecular Driver Group", y = "Centralised enrichment score") +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    title = element_text(size = 10.5)
  )

reviewer_comment_figure_tds <- patchwork::wrap_plots(
  reviewer_plots_tds,
  design = "AB\nCD\nEE",
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "fTFC Signatures and TCGA Clinical Outcome (incl. TDS)",
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 11, face = "bold"))
  )

print(reviewer_comment_figure_tds)

ggsave(
  filename = file.path(
    plotDir,
    "TCGA_reviewer4_figure_selected_clinical_outcomes_withTDS.pdf"
  ),
  plot = reviewer_comment_figure_tds,
  width = 16,
  height = 14
)


dd = add_tds_rows(d)
dd$group = ifelse(dd$group == 'RAS', 'RAS', 'not-RAS')
dd$ajcc_pathologic_n <- factor(
  dd$ajcc_pathologic_n,
  c("N0", "N1", "N1a", "N1b")
)
plot_tcga_associations(
  data = dd,
  vars = "ajcc_pathologic_n",
  score_col = "normalised_score",
  module_keep = c("fTFC1", "fTFC2", "TDS"),
  outfile = NULL,
  print_plot = FALSE
) +
  facet_grid(
    ~ moduleType + group,
    scales = "free_x",
    space = "free"
  ) +
  labs(title = "Molecular Driver Group", y = "Centralised enrichment score") +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    title = element_text(size = 10.5)
  )

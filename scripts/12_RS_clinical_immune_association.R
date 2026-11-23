###############################################################################
## Stage 12 Association of RS with Clinical and Immune Characteristics
##
## Inputs:
## data/expr_ESCA_tumor.csv (RSEM log2(x+1))
## data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv (clinical)
## tables/Tab10_risk_scores_per_sample.csv (RiskScore, RiskGroup)
## tables/Tab8A_ESTIMATE_scores.csv
## tables/Tab8B_CIBERSORT_results.csv
## tables/Tab8C_ssGSEA_Bindea.csv (long: cell_type x samples)
##
## Outputs:
## tables/Tab12_RS_clinical_distribution.csv
## tables/Tab12_RS_clinical_test.csv
## tables/Tab12_RS_immune_spearman.csv
## tables/Tab12_RS_immune_wilcoxon.csv
## figures/Fig8A_RS_clinical_barplot.{pdf,png}
## figures/Fig8B_RS_immune_heatmap.{pdf,png}
## figures/Fig8C_RS_immune_spearman_forest.{pdf,png}
## figures/Fig8D_RS_immune_wilcoxon_boxplot.{pdf,png}
## logs/12_RS_association_YYYYMMDD.log
## logs/sessionInfo_stage12.txt
###############################################################################

# run from the repository root (see README.md)

INPUT_EXPR <- "data/expr_ESCA_tumor.csv"
INPUT_CLI <- "data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv"
INPUT_RS <- "tables/Tab10_risk_scores_per_sample.csv"
INPUT_EST <- "tables/Tab8A_ESTIMATE_scores.csv"
INPUT_CIB <- "tables/Tab8B_CIBERSORT_results.csv"
INPUT_SSG <- "tables/Tab8C_ssGSEA_Bindea.csv"

TBLDIR <- "tables"
FIGDIR <- "figures"
LOGDIR <- "logs"
for (d in c(TBLDIR, FIGDIR, LOGDIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

SEED <- 1262118388
set.seed(SEED)

# Logging
log_file <- file.path(LOGDIR, "12_diagnostics.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
cat(sprintf("=== Stage 12 RS-Association Pipeline start (%s) ===\n\n", Sys.time()))

# Packages
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
  library(ragg)
})

# Color palette
COL_RS <- c(Low = "#377EB8", High = "#E41A1C") # ColorBrewer Set1

# 1. Load inputs
cat(">>> STEP 1: Load inputs\n")
expr <- fread(INPUT_EXPR)
cli <- fread(INPUT_CLI)
rs <- fread(INPUT_RS)
est <- fread(INPUT_EST)
cib <- fread(INPUT_CIB)
ssg <- fread(INPUT_SSG)

cat(sprintf(" expr : %d genes %d cols (incl gene)\n", nrow(expr), ncol(expr)))
cat(sprintf(" clinical : %d rows %d cols\n", nrow(cli), ncol(cli)))
cat(sprintf(" RS samples : %d (RiskGroup levels: %s)\n",
            nrow(rs), paste(unique(rs$RiskGroup), collapse=",")))
cat(sprintf(" ESTIMATE : %d %d\n", nrow(est), ncol(est)))
cat(sprintf(" CIBERSORT : %d %d\n", nrow(cib), ncol(cib)))
cat(sprintf(" ssGSEA Bindea : %d %d (long format)\n\n", nrow(ssg), ncol(ssg)))

# 2. Sample alignment (RS ESTIMATE CIBERSORT ssGSEA clinical)
cat(">>> STEP 2: Sample alignment\n")
ssg_samples <- setdiff(colnames(ssg), "cell_type")
samples_keep <- Reduce(intersect, list(rs$sampleID,
                                       est$sampleID,
                                       cib$sampleID,
                                       ssg_samples,
                                       cli$sampleID))
N <- length(samples_keep)
cat(sprintf(" RS samples : %d\n", nrow(rs)))
cat(sprintf(" ESTIMATE samples : %d\n", nrow(est)))
cat(sprintf(" CIBERSORT samples : %d\n", nrow(cib)))
cat(sprintf(" ssGSEA Bindea samples: %d\n", length(ssg_samples)))
cat(sprintf(" Clinical samples : %d\n", nrow(cli)))
cat(sprintf(" 5-way intersection (N): %d\n\n", N))
stopifnot(N > 50)

# Order
rs <- rs[match(samples_keep, sampleID)]
est <- est[match(samples_keep, sampleID)]
cib <- cib[match(samples_keep, sampleID)]

# Risk group factor (ensure Low first)
rs[, RiskGroup := factor(RiskGroup, levels = c("Low", "High"))]
n_low <- sum(rs$RiskGroup == "Low")
n_high <- sum(rs$RiskGroup == "High")
cat(sprintf(" Risk group: Low=%d, High=%d\n\n", n_low, n_high))

# 3. STEP 3: RS group clinical distribution
cat(">>> STEP 3: RS group clinical variable distribution\n")

CLIN_CONT <- c("age_at_initial_pathologic_diagnosis")
CLIN_CAT <- c("gender", "pathologic_stage", "neoplasm_histologic_grade",
               "residual_tumor", "histological_type",
               "alcohol_history_documented", "tobacco_smoking_history",
               "reflux_history", "barretts_esophagus")

cli_keep <- cli[match(samples_keep, sampleID), c("sampleID", CLIN_CONT, CLIN_CAT), with = FALSE]
cli_keep <- merge(cli_keep, rs[, .(sampleID, RiskScore, RiskGroup)], by = "sampleID")

# Distribution table (Table 1)
build_dist_row_cont <- function(varname, dt) {
  v_low <- as.numeric(dt[RiskGroup == "Low", get(varname)])
  v_high <- as.numeric(dt[RiskGroup == "High", get(varname)])
  data.table(
    variable = varname,
    level = "median (IQR)",
    Low = sprintf("%.1f (%.1f-%.1f)",
                       median(v_low, na.rm=TRUE),
                       quantile(v_low, 0.25, na.rm=TRUE),
                       quantile(v_low, 0.75, na.rm=TRUE)),
    High = sprintf("%.1f (%.1f-%.1f)",
                       median(v_high, na.rm=TRUE),
                       quantile(v_high, 0.25, na.rm=TRUE),
                       quantile(v_high, 0.75, na.rm=TRUE)),
    Total = sprintf("%.1f (%.1f-%.1f)",
                       median(dt[[varname]], na.rm=TRUE),
                       quantile(dt[[varname]], 0.25, na.rm=TRUE),
                       quantile(dt[[varname]], 0.75, na.rm=TRUE))
  )
}

build_dist_rows_cat <- function(varname, dt) {
  v <- dt[[varname]]
  v[is.na(v) | v == "" | v == "[Not Available]" | v == "[Not Applicable]"] <- NA
  lvls <- sort(unique(v[!is.na(v)]))
  rows <- lapply(lvls, function(lv) {
    n_low <- sum(dt$RiskGroup == "Low" & v == lv, na.rm = TRUE)
    n_high <- sum(dt$RiskGroup == "High" & v == lv, na.rm = TRUE)
    n_tot <- n_low + n_high
    data.table(
      variable = varname,
      level = as.character(lv),
      Low = sprintf("%d (%.1f%%)", n_low, 100*n_low / max(1, n_low_total <- sum(dt$RiskGroup == "Low" & !is.na(v)))),
      High = sprintf("%d (%.1f%%)", n_high, 100*n_high / max(1, sum(dt$RiskGroup == "High" & !is.na(v)))),
      Total = sprintf("%d (%.1f%%)", n_tot, 100*n_tot / max(1, sum(!is.na(v))))
    )
  })
  if (sum(is.na(v)) > 0) {
    rows[[length(rows)+1]] <- data.table(
      variable = varname, level = "NA",
      Low = sprintf("%d", sum(dt$RiskGroup == "Low" & is.na(v))),
      High = sprintf("%d", sum(dt$RiskGroup == "High" & is.na(v))),
      Total = sprintf("%d", sum(is.na(v)))
    )
  }
  rbindlist(rows)
}

dist_rows <- list()
for (v in CLIN_CONT) dist_rows[[v]] <- build_dist_row_cont(v, cli_keep)
for (v in CLIN_CAT) dist_rows[[v]] <- build_dist_rows_cat(v, cli_keep)
dist_dt <- rbindlist(dist_rows)
fwrite(dist_dt, file.path(TBLDIR, "Tab12_RS_clinical_distribution.csv"))
cat(sprintf(" [OK] Tab12_RS_clinical_distribution.csv (%d rows)\n", nrow(dist_dt)))

# Statistical tests
test_clin <- function(varname, dt, is_cont) {
  v <- dt[[varname]]
  if (is_cont) {
    v_low <- as.numeric(v[dt$RiskGroup == "Low"])
    v_high <- as.numeric(v[dt$RiskGroup == "High"])
    p <- tryCatch(wilcox.test(v_low, v_high, exact = FALSE)$p.value,
                  error = function(e) NA_real_)
    data.table(variable = varname, test_type = "Wilcoxon", n_levels = NA_integer_,
               p_value = p)
  } else {
    v[is.na(v) | v == "" | v == "[Not Available]" | v == "[Not Applicable]"] <- NA
    sub <- dt[!is.na(v)]
    sub_v <- v[!is.na(v)]
    tab <- table(sub$RiskGroup, sub_v)
    if (any(dim(tab) < 2)) {
      return(data.table(variable = varname, test_type = "Skipped (single level)",
                        n_levels = ncol(tab), p_value = NA_real_))
    }
    expected_min <- min(suppressWarnings(chisq.test(tab)$expected), na.rm=TRUE)
    if (expected_min < 5) {
      p <- tryCatch(fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value,
                    error = function(e) NA_real_)
      tt <- "Fisher (sim)"
    } else {
      p <- suppressWarnings(chisq.test(tab)$p.value)
      tt <- "Chi-square"
    }
    data.table(variable = varname, test_type = tt, n_levels = ncol(tab), p_value = p)
  }
}

test_rows <- list()
for (v in CLIN_CONT) test_rows[[v]] <- test_clin(v, cli_keep, is_cont = TRUE)
for (v in CLIN_CAT) test_rows[[v]] <- test_clin(v, cli_keep, is_cont = FALSE)
test_dt <- rbindlist(test_rows)
test_dt[, p_adj_BH := p.adjust(p_value, method = "BH")]
test_dt[, sig := fifelse(is.na(p_value), "ns",
                  fifelse(p_value < 0.001, "***",
                  fifelse(p_value < 0.01, "**",
                  fifelse(p_value < 0.05, "*", "ns"))))]
fwrite(test_dt, file.path(TBLDIR, "Tab12_RS_clinical_test.csv"))
cat(" [OK] Tab12_RS_clinical_test.csv\n")
print(test_dt)
cat("\n")

# Fig 8A: stacked barplot of categorical vars
cat(" >>> Fig8A: clinical categorical stacked barplot\n")
cat_long_list <- lapply(CLIN_CAT, function(v) {
  vec <- cli_keep[[v]]
  vec[is.na(vec) | vec == "" | vec == "[Not Available]" | vec == "[Not Applicable]"] <- NA
  data.table(variable = v,
             level = ifelse(is.na(vec), "NA", as.character(vec)),
             RiskGroup= cli_keep$RiskGroup)
})
cat_long <- rbindlist(cat_long_list)
cat_long <- cat_long[level != "NA"] # drop NA bars

# Sig label per variable
sig_lbl <- test_dt[, .(variable, sig, p_adj_BH)]
cat_long_sum <- cat_long[, .N, by = .(variable, level, RiskGroup)]
cat_long_sum <- cat_long_sum[, .(level = level, N = N,
                                  prop = N / sum(N)),
                             by = .(variable, RiskGroup)]

p_clin <- ggplot(cat_long_sum, aes(x = RiskGroup, y = prop, fill = level)) +
  geom_col(width = 0.8, color = "black", linewidth = 0.2) +
  facet_wrap(~ variable, scales = "free_x", ncol = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_fill_brewer(palette = "Set3") +
  geom_text(data = sig_lbl, aes(x = 1.5, y = 1.05, label = sig),
            inherit.aes = FALSE, size = 3.5, fontface = "bold") +
  labs(x = NULL, y = "Proportion") +
  theme_classic(base_size = 10) +
  theme(legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 7),
        strip.text = element_text(face = "bold", size = 9))

cairo_pdf(file.path(FIGDIR, "Fig8A_RS_clinical_barplot.pdf"), width = 12, height = 8)
print(p_clin); dev.off()
ragg::agg_png(file.path(FIGDIR, "Fig8A_RS_clinical_barplot.png"), width = 12, height = 8, units = "in", res = 600)
print(p_clin); dev.off()
cat(" [OK] Fig8A_RS_clinical_barplot.{pdf,png}\n\n")

# 4. STEP 4 Build immune feature matrix (102 features)
cat(">>> STEP 4: Build 102-feature immune matrix\n")

# 4.1 ESTIMATE 3 (already aligned)
est_mat <- as.matrix(est[, .(StromalScore, ImmuneScore, ESTIMATEScore)])
rownames(est_mat) <- est$sampleID
est_features <- colnames(est_mat)
cat(sprintf(" ESTIMATE features : %d\n", length(est_features)))

# 4.2 CIBERSORT 22 cell types (columns 4..25 in cib; metadata = sampleID, cluster, histology)
cib_meta_cols <- c("sampleID", "cluster", "histology", "P-value", "Correlation", "RMSE")
cib_cell_cols <- setdiff(colnames(cib), cib_meta_cols)
stopifnot(length(cib_cell_cols) == 22)
cib_mat <- as.matrix(cib[, ..cib_cell_cols])
rownames(cib_mat) <- cib$sampleID
# only confident deconvolutions (permutation p < 0.05; Methods 2.6) enter the risk-score analyses, as in the
# cluster comparison: the fractions of the other samples (one borderline, p = 0.06) are set missing
cib_mat[cib$`P-value` >= 0.05, ] <- NA
cat(sprintf(" CIBERSORT features : %d (confident samples: %d of %d)\n", length(cib_cell_cols),
            sum(cib$`P-value` < 0.05), nrow(cib)))

# 4.3 ssGSEA Bindea 28 cell types transpose
ssg_keep <- ssg[, c("cell_type", samples_keep), with = FALSE]
ssg_mat <- t(as.matrix(ssg_keep[, -1]))
colnames(ssg_mat) <- ssg_keep$cell_type
# Disambiguate Neutrophils / Eosinophils (also in CIBERSORT 22) append " (Bindea)"
dup_names <- intersect(colnames(ssg_mat), cib_cell_cols)
if (length(dup_names) > 0) {
  cat(sprintf(" ssGSEA Bindea CIBERSORT cell type %ddisambiguation: %s\n",
              length(dup_names), paste(dup_names, collapse = ", ")))
  for (nm in dup_names) {
    colnames(ssg_mat)[colnames(ssg_mat) == nm] <- paste0(nm, " (Bindea)")
  }
}
ssg_features <- colnames(ssg_mat)
cat(sprintf(" ssGSEA Bindea cell types: %d\n", length(ssg_features)))

# 4.4 Immune checkpoint 30 genes
liu_to_ncbi <- c(
  "PD-L1" = "CD274",
  "SLAMF4" = "CD244",
  "TIM3" = "HAVCR2",
  "PD-1" = "PDCD1",
  "PD-L2" = "PDCD1LG2",
  "B7-H3" = "CD276",
  "B7.2" = "CD86",
  "GITR" = "TNFRSF18",
  "DR3" = "TNFRSF25",
  "OX40" = "TNFRSF4",
  "CD30" = "TNFRSF8",
  "4-1BB" = "TNFRSF9",
  "APRIL" = "TNFSF13",
  "BAFF" = "TNFSF13B",
  "LIGHT" = "TNFSF14",
  "TL1A" = "TNFSF15",
  "4-1BB-L" = "TNFSF9"
)
ckpt_liu <- c("PD-L1","SLAMF4","CD96","CSF1R","TIM3","IL10RB","KIR2DL1","KIR2DL3",
              "LGALS9","PD-1","PD-L2","TGFB1","B7-H3","CD48","CD70","B7.2",
              "CXCR4","IL6R","MICB","TMIGD2","GITR","DR3","OX40","CD30","4-1BB",
              "APRIL","BAFF","LIGHT","TL1A","4-1BB-L")
stopifnot(length(ckpt_liu) == 30)
ckpt_ncbi <- ifelse(ckpt_liu %in% names(liu_to_ncbi),
                    unname(liu_to_ncbi[ckpt_liu]),
                    ckpt_liu)
# expression matrix lookup
expr_genes <- expr$sample
ckpt_match <- match(ckpt_ncbi, expr_genes)
n_match <- sum(!is.na(ckpt_match))
cat(sprintf(" Immune checkpoint matched: %d / 30\n", n_match))
if (n_match < 30) {
  miss <- ckpt_liu[is.na(ckpt_match)]
  cat(sprintf(" Missing: %s\n", paste(miss, collapse = ", ")))
}
ckpt_keep_idx <- !is.na(ckpt_match)
ckpt_mat <- t(as.matrix(expr[ckpt_match[ckpt_keep_idx], ..samples_keep]))
colnames(ckpt_mat) <- ckpt_liu[ckpt_keep_idx]
ckpt_features <- colnames(ckpt_mat)
cat(sprintf(" Immune checkpoint features: %d\n", length(ckpt_features)))

# 4.5 HLA family 19 genes
hla_genes <- c("HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB",
               "HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2",
               "HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-E","HLA-F","HLA-G")
stopifnot(length(hla_genes) == 19)
hla_match <- match(hla_genes, expr_genes)
n_hla <- sum(!is.na(hla_match))
cat(sprintf(" HLA family matched: %d / 19\n", n_hla))
if (n_hla < 19) {
  miss <- hla_genes[is.na(hla_match)]
  cat(sprintf(" Missing: %s\n", paste(miss, collapse = ", ")))
}
hla_keep_idx <- !is.na(hla_match)
hla_mat <- t(as.matrix(expr[hla_match[hla_keep_idx], ..samples_keep]))
colnames(hla_mat) <- hla_genes[hla_keep_idx]
hla_features <- colnames(hla_mat)
cat(sprintf(" HLA features: %d\n", length(hla_features)))

# 4.6 Combined feature matrix
feat_mat <- cbind(est_mat[samples_keep, , drop=FALSE],
                  cib_mat[samples_keep, , drop=FALSE],
                  ssg_mat[samples_keep, , drop=FALSE],
                  ckpt_mat[samples_keep, , drop=FALSE],
                  hla_mat[samples_keep, , drop=FALSE])
storage.mode(feat_mat) <- "numeric"

feat_category <- c(rep("ESTIMATE", length(est_features)),
                   rep("CIBERSORT", length(cib_cell_cols)),
                   rep("ssGSEA_Bindea", length(ssg_features)),
                   rep("Immune_checkpoint", length(ckpt_features)),
                   rep("HLA", length(hla_features)))
names(feat_category) <- colnames(feat_mat)

cat(sprintf("\n Combined feature matrix: %d samples %d features\n", nrow(feat_mat), ncol(feat_mat)))
cat(" Per-category counts:\n"); print(table(feat_category))
cat("\n")

# RS vector aligned
rs_vec <- rs$RiskScore
rg_vec <- rs$RiskGroup
names(rs_vec) <- names(rg_vec) <- rs$sampleID

# 5. STEP 5: RS immune Spearman
cat(">>> STEP 5: RS immune Spearman correlation\n")

spear_rows <- lapply(colnames(feat_mat), function(f) {
  x <- feat_mat[, f]
  ct <- suppressWarnings(cor.test(rs_vec, x, method = "spearman", exact = FALSE))
  data.table(
    feature_category = unname(feat_category[f]),
    feature = f,
    rho = unname(ct$estimate),
    p_value = ct$p.value,
    n = sum(!is.na(x))
  )
})
spear_dt <- rbindlist(spear_rows)
spear_dt[, p_adj_BH := p.adjust(p_value, method = "BH")]
spear_dt[, sig := fifelse(is.na(p_value), "ns",
                   fifelse(p_value < 0.001, "***",
                   fifelse(p_value < 0.01, "**",
                   fifelse(p_value < 0.05, "*", "ns"))))]
setorder(spear_dt, p_value)
fwrite(spear_dt, file.path(TBLDIR, "Tab12_RS_immune_spearman.csv"))
cat(sprintf(" [OK] Tab12_RS_immune_spearman.csv\n"))
cat(sprintf(" Significant (raw p<0.05): %d / %d\n", sum(spear_dt$p_value < 0.05, na.rm=TRUE), nrow(spear_dt)))
cat(" Top 10 by |rho|:\n")
print(head(spear_dt[order(-abs(rho))], 10))
cat("\n")

# 6. STEP 6 Fig 8B: heatmap (RS-sorted samples top features)
cat(">>> STEP 6: Fig8B RS-immune heatmap\n")
TOP_HM <- 50
top_feats <- head(spear_dt[order(-abs(rho))]$feature, TOP_HM)
mat_hm <- t(scale(feat_mat[, top_feats], center = TRUE, scale = TRUE)) # row z-score
ord_rs <- order(rs_vec) # samples sorted by RS

ha_top <- HeatmapAnnotation(
  RiskScore = rs_vec[ord_rs],
  RiskGroup = factor(rg_vec[ord_rs], levels = c("Low", "High")),
  col = list(
    RiskScore = colorRamp2(c(min(rs_vec), 0, max(rs_vec)), c("#377EB8","white","#E41A1C")),
    RiskGroup = COL_RS
  ),
  annotation_legend_param = list(RiskGroup = list(title = "RiskGroup"))
)
ha_left <- rowAnnotation(
  Category = feat_category[top_feats],
  col = list(Category = c(
    ESTIMATE = "#1B9E77", CIBERSORT = "#D95F02", ssGSEA_Bindea = "#7570B3",
    Immune_checkpoint = "#E7298A", HLA = "#66A61E"
  )),
  show_annotation_name = FALSE
)

ht <- Heatmap(
  mat_hm[, ord_rs],
  name = "z-score",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC","white","#B2182B")),
  top_annotation = ha_top,
  left_annotation = ha_left,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 7),
  column_title = NULL,
  heatmap_legend_param = list(title = "z-score")
)

cairo_pdf(file.path(FIGDIR, "Fig8B_RS_immune_heatmap.pdf"), width = 11, height = 9)
draw(ht, merge_legend = TRUE); dev.off()
ragg::agg_png(file.path(FIGDIR, "Fig8B_RS_immune_heatmap.png"), width = 11, height = 9, units = "in", res = 600)
draw(ht, merge_legend = TRUE); dev.off()
cat(" [OK] Fig8B_RS_immune_heatmap.{pdf,png}\n\n")

# 7. STEP 7 Fig 8C: Spearman forest plot
cat(">>> STEP 7: Fig8C Spearman forest plot (top 30 |rho|)\n")

# 95% CI for Spearman via Fisher z transform
fisher_ci <- function(rho, n) {
  if (is.na(rho) || n < 4) return(c(NA, NA))
  z <- atanh(rho); se <- 1 / sqrt(n - 3)
  ci <- tanh(z + c(-1, 1) * 1.96 * se)
  ci
}
spear_top30 <- head(spear_dt[order(-abs(rho))], 30)
cis <- t(mapply(fisher_ci, spear_top30$rho, spear_top30$n))
spear_top30[, c("rho_low","rho_high") := list(cis[, 1], cis[, 2])]
spear_top30[, feature := factor(feature, levels = rev(feature))] # for forest order

p_for <- ggplot(spear_top30, aes(x = rho, y = feature, color = feature_category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = rho_low, xmax = rho_high), height = 0.25, linewidth = 0.4) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c(
    ESTIMATE = "#1B9E77", CIBERSORT = "#D95F02", ssGSEA_Bindea = "#7570B3",
    Immune_checkpoint = "#E7298A", HLA = "#66A61E"
  ), name = "Category") +
  labs(x = "Spearman (RS vs feature)", y = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(size = 7))

cairo_pdf(file.path(FIGDIR, "Fig8C_RS_immune_spearman_forest.pdf"), width = 8.5, height = 9)
print(p_for); dev.off()
ragg::agg_png(file.path(FIGDIR, "Fig8C_RS_immune_spearman_forest.png"), width = 8.5, height = 9, units = "in", res = 600)
print(p_for); dev.off()
cat(" [OK] Fig8C_RS_immune_spearman_forest.{pdf,png}\n\n")

# 8. STEP 8: RS group immune Wilcoxon
cat(">>> STEP 8: RS group immune Wilcoxon\n")

wilcox_rows <- lapply(colnames(feat_mat), function(f) {
  x <- feat_mat[, f]
  v_low <- x[rg_vec == "Low"]
  v_high <- x[rg_vec == "High"]
  med_low <- median(v_low, na.rm = TRUE)
  med_high <- median(v_high, na.rm = TRUE)
  # log2FC: only meaningful for non-negative scores. CIBERSORT/ssGSEA/Bindea/checkpoint/HLA all log2(x+1) or proportions/scores; ESTIMATE can be negative.
  fc <- if (med_low > 0 && med_high > 0) log2(med_high / med_low) else (med_high - med_low)
  p <- tryCatch(wilcox.test(v_low, v_high, exact = FALSE)$p.value,
                 error = function(e) NA_real_)
  data.table(
    feature_category = unname(feat_category[f]),
    feature = f,
    n_low = sum(!is.na(v_low)),
    n_high = sum(!is.na(v_high)),
    median_low = med_low,
    median_high = med_high,
    log2FC_or_diff = fc,
    p_value = p
  )
})
wilcox_dt <- rbindlist(wilcox_rows)
wilcox_dt[, p_adj_BH := p.adjust(p_value, method = "BH")]
wilcox_dt[, sig := fifelse(is.na(p_value), "ns",
                    fifelse(p_value < 0.001, "***",
                    fifelse(p_value < 0.01, "**",
                    fifelse(p_value < 0.05, "*", "ns"))))]
setorder(wilcox_dt, p_value)
fwrite(wilcox_dt, file.path(TBLDIR, "Tab12_RS_immune_wilcoxon.csv"))
n_sig_w <- sum(wilcox_dt$p_value < 0.05, na.rm=TRUE)
cat(sprintf(" [OK] Tab12_RS_immune_wilcoxon.csv\n"))
cat(sprintf(" Significant (raw p<0.05): %d / %d\n", n_sig_w, nrow(wilcox_dt)))
cat(" Top 10 by raw p:\n")
print(head(wilcox_dt, 10))
cat("\n")

# 9. STEP 9 Fig 8D: Wilcoxon boxplot for raw p-significant features
cat(">>> STEP 9: Fig8D raw-p significant feature boxplot\n")
sig_feats <- wilcox_dt[!is.na(p_value) & p_value < 0.05, feature]
TOP_BOX <- min(30, length(sig_feats))
if (TOP_BOX == 0) {
  cat(" No raw-p significant features skipping Fig8D\n")
} else {
  show_feats <- head(wilcox_dt[!is.na(p_value) & p_value < 0.05][order(p_value), feature], TOP_BOX)
  long_dt <- rbindlist(lapply(show_feats, function(f) {
    data.table(feature = f, RiskGroup = rg_vec, value = feat_mat[, f])
  }))
  long_dt[, feature := factor(feature, levels = show_feats)]

  # significance label per facet
  lbl_dt <- wilcox_dt[feature %in% show_feats, .(feature, sig, p_adj_BH)]
  lbl_dt[, feature := factor(feature, levels = show_feats)]
  ymax_dt <- long_dt[, .(y_max = max(value, na.rm = TRUE)), by = feature]
  lbl_dt <- merge(lbl_dt, ymax_dt, by = "feature")
  lbl_dt[, y_pos := y_max + 0.08 * abs(y_max) + 0.02]

  p_box <- ggplot(long_dt, aes(x = RiskGroup, y = value, fill = RiskGroup)) +
    geom_boxplot(outlier.size = 0.5, color = "black", linewidth = 0.3) +
    geom_text(data = lbl_dt, aes(x = 1.5, y = y_pos, label = sig),
              inherit.aes = FALSE, size = 3, fontface = "bold") +
    facet_wrap(~ feature, scales = "free_y", ncol = 6) +
    scale_fill_manual(values = COL_RS) +
    labs(x = NULL, y = "Feature value") +
    theme_bw(base_size = 9) +
    theme(legend.position = "top",
          strip.text = element_text(face = "bold", size = 7),
          axis.text.x = element_text(face = "bold", size = 7))

  W <- 14
  H <- ceiling(TOP_BOX / 6) * 1.8 + 1
  cairo_pdf(file.path(FIGDIR, "Fig8D_RS_immune_wilcoxon_boxplot.pdf"), width = W, height = H)
  print(p_box); dev.off()
  ragg::agg_png(file.path(FIGDIR, "Fig8D_RS_immune_wilcoxon_boxplot.png"), width = W, height = H, units = "in", res = 600)
  print(p_box); dev.off()
  cat(sprintf(" [OK] Fig8D_RS_immune_wilcoxon_boxplot.{pdf,png} (%d sig features shown)\n", TOP_BOX))
}
cat("\n")

# 10. Summary report
cat("========================================================================\n")
cat(" Stage 12 RS-Association Summary Report \n")
cat("========================================================================\n")
cat(sprintf(" N samples (5-way intersect) : %d (Low=%d, High=%d)\n", N, n_low, n_high))
cat(sprintf(" Clinical variables tested : %d (cont=%d, cat=%d)\n",
            length(CLIN_CONT)+length(CLIN_CAT), length(CLIN_CONT), length(CLIN_CAT)))
cat(sprintf(" Significant (raw p<0.05) : %d\n", sum(test_dt$p_value < 0.05, na.rm=TRUE)))
cat(sprintf(" Immune features : %d (3+22+28+%d+%d)\n",
            ncol(feat_mat), length(ckpt_features), length(hla_features)))
cat(sprintf(" Spearman significant (raw p<0.05) : %d / %d\n",
            sum(spear_dt$p_value < 0.05, na.rm=TRUE), nrow(spear_dt)))
cat(sprintf(" Wilcoxon significant (raw p<0.05) : %d / %d\n",
            n_sig_w, nrow(wilcox_dt)))
cat("========================================================================\n\n")

# 11. sessionInfo
sess_file <- file.path(LOGDIR, "sessionInfo/sessionInfo_stage12.txt")
sink(sess_file)
cat(sprintf("Session info captured: %s\n\n", Sys.time()))
print(sessionInfo())
sink()
cat(sprintf(">>> sessionInfo saved: %s\n", sess_file))

cat(sprintf("\n=== Stage 12 RS-Association Pipeline done (%s) ===\n", Sys.time()))

sink(type = "message")
sink()
close(log_con)

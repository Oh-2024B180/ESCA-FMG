###############################################################################
## Stage 7: Molecular Subtyping (ConsensusClusterPlus, 18 FMG)
##
## Reference:
## 7A: Input feature -- 18 FMG TCGA-ESCA expression, z-score (gene-wise)
## 7B: CCP parameters -- maxK=9, reps=1000, pItem=0.8, pFeature=1.0,
## clusterAlg="pam", distance="spearman",
## innerLinkage="ward.D2", finalLinkage="ward.D2",
## seed=1262118388
## 7C: Optimal k -- multi-criteria decision (computed/logged in STEP 6):
##   S1 cluster balance  -- min cluster >= 10% of N        [Senbabaoglu 2014, ref 43]
##   S2 PAC local minimum -- smallest k that is a local min [Senbabaoglu 2014, ref 43]
##   corroborated by monotone-decreasing relative CDF delta-area (Monti 2003)
## Expected/observed: k=2 (global-lowest-PAC k=9 rejected: min cluster < 10%)
## 7D: Subtype validation -- PCA + KM (log-rank p<0.05) + 18 FMG heatmap
##
## Reference PMIDs:
## 20427518 : Wilkerson & Hayes 2010, Bioinformatics (ConsensusClusterPlus)
## 25158761 : Senbabaoglu 2014, Sci Rep (PAC critique)
## 27207943 : Gu 2016, Bioinformatics (ComplexHeatmap)
## 28398314 : Vivian 2017, Nat Biotechnol (UCSC Toil expr)
##
## Inputs:
## data/expr_ESCA_tumor.csv (UCSC Xena log2(RSEM+1), 58,581 x 181)
## tables/Tab1_2_FMG_list_primary.csv (18 Primary FMG)
## data/survival_ESCA_survival.csv (OS, OS.time)
## data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv (sample_type, histological_type)
##
## Outputs:
## tables/Tab7_cluster_assignment.csv (sampleID, cluster, OS/OS.time, histology)
## tables/Tab7_cluster_summary.csv (per-cluster: n, ESCC/EAC/Other, median OS, log-rank p)
## figures/Fig3A_CDF_curves.{pdf,png}
## figures/Fig3B_PAC_score.{pdf,png}
## figures/Fig3C_consensus_matrix_k2.{pdf,png}
## figures/Fig3D_PCA_clusters.{pdf,png}
## figures/Fig3E_KM_clusters.{pdf,png}
## figures/Fig3F_FMG_heatmap_clusters.{pdf,png}
## figures/CCP_native/ (CCP native output dir)
##
## logs/07_consensus_clustering_YYYYMMDD.log
## logs/sessionInfo_stage7.txt
##
## Notes:
## - z-score: gene-wise (each FMG row scaled to mean 0, sd 1); redundant with the Spearman distance,
## but applied as specified in 7A.
## - ConsensusClusterPlus saves its own PDFs to plot_dir; Fig3A/3B/3C are drawn separately here.
## - Color palette: ColorBrewer Dark2 (color-blind-safe, CLAUDE.md §7-4)
###############################################################################

# 0. setup
# run from the repository root (see README.md)

INPUT_EXPR <- "data/expr_ESCA_tumor.csv"
INPUT_FMG <- "tables/Tab1_2_FMG_list_primary.csv"
INPUT_SURV <- "data/survival_ESCA_survival.csv"
INPUT_CLI <- "data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv"

TBLDIR <- "tables"
FIGDIR <- "figures"
CCPDIR <- file.path(FIGDIR, "CCP_native")
LOGDIR <- "logs"

for (d in c(TBLDIR, FIGDIR, CCPDIR, LOGDIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

SEED <- 1262118388
MAX_K <- 9
REPS <- 1000
P_ITEM <- 0.8
P_FEATURE <- 1.0
CLUSTER_ALG <- "pam"
DISTANCE <- "spearman"
INNER_LINK <- "ward.D2"
FINAL_LINK <- "ward.D2"
K_EXPECTED <- 2 # expected k; K_FINAL is computed data-driven in STEP 6

# Color palette: ColorBrewer Dark2 (color-blind-safe)
CLUSTER_COLORS <- c("C1" = "#1B9E77", "C2" = "#D95F02")
HIST_COLORS <- c("ESCC" = "#7570B3", "EAC" = "#E7298A", "Other" = "#66A61E")

set.seed(SEED)

# Logging start
log_file <- file.path(LOGDIR, "07_diagnostics.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat(sprintf("=== Stage 7 Pipeline start (%s) ===\n\n", Sys.time()))

# 1. load
suppressPackageStartupMessages({
  library(data.table)
  library(ConsensusClusterPlus)
  library(survival)
  library(survminer)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggplot2)
  library(RColorBrewer)
  library(ragg)
})

cat(sprintf("ConsensusClusterPlus version : %s\n", packageVersion("ConsensusClusterPlus")))
cat(sprintf("ComplexHeatmap version : %s\n", packageVersion("ComplexHeatmap")))
cat(sprintf("survival version : %s\n\n", packageVersion("survival")))

# 2. Load input data
cat(">>> STEP 2: Load inputs\n")
expr <- fread(INPUT_EXPR) # 1st col = "sample" (gene name), rest = TCGA sample IDs
fmg <- fread(INPUT_FMG)
surv <- fread(INPUT_SURV)
cli <- fread(INPUT_CLI)

cat(sprintf(" expr rows (genes) : %d\n", nrow(expr)))
cat(sprintf(" expr cols (samples+1) : %d\n", ncol(expr)))
cat(sprintf(" FMG count : %d\n", nrow(fmg)))
cat(sprintf(" survival rows : %d\n", nrow(surv)))
cat(sprintf(" clinical rows : %d\n\n", nrow(cli)))

# 3. Sample alignment
cat(">>> STEP 3: Sample alignment\n")
expr_samples <- setdiff(colnames(expr), "sample")
primary_tumor_samples <- cli[sample_type == "Primary Tumor", sampleID]
surv_samples <- surv$sample

samples_keep <- Reduce(intersect, list(expr_samples, primary_tumor_samples, surv_samples))
N <- length(samples_keep)
cat(sprintf(" expr samples : %d\n", length(expr_samples)))
cat(sprintf(" Primary Tumor (clinical): %d\n", length(primary_tumor_samples)))
cat(sprintf(" Survival samples : %d\n", length(surv_samples)))
cat(sprintf(" 3-way intersection (N) : %d\n\n", N))

if (N < 50) stop("Sample intersection too small; check inputs.")

# 4. The 18 FMG matrix + z-score
cat(">>> STEP 4: Build 18-FMG expression matrix and z-score\n")
fmg_genes <- fmg$gene
expr_fmg <- expr[match(fmg_genes, expr$sample), c("sample", samples_keep), with = FALSE]

# missingFMG
missing_fmg <- fmg_genes[!fmg_genes %in% expr$sample]
if (length(missing_fmg) > 0) {
  cat(sprintf(" Missing FMG (not in expr): %s\n", paste(missing_fmg, collapse = ", ")))
} else {
  cat(" All 18 FMG present in expr matrix.\n")
}

mat_raw <- as.matrix(expr_fmg[, -1])
rownames(mat_raw) <- expr_fmg$sample
storage.mode(mat_raw) <- "numeric"

# Gene-wise z-score
mat_z <- t(scale(t(mat_raw), center = TRUE, scale = TRUE))
cat(sprintf(" Matrix dim (FMG x sample): %d x %d\n", nrow(mat_z), ncol(mat_z)))
cat(sprintf(" Per-gene z-score sanity (mean ~ 0, sd ~ 1):\n"))
print(summary(rowMeans(mat_z)))
print(summary(apply(mat_z, 1, sd)))
cat("\n")

# 5. Run ConsensusClusterPlus
cat(sprintf(">>> STEP 5: ConsensusClusterPlus (maxK=%d, reps=%d) -- this may take several minutes\n",
            MAX_K, REPS))
ccp_res <- ConsensusClusterPlus(
  d = mat_z,
  maxK = MAX_K,
  reps = REPS,
  pItem = P_ITEM,
  pFeature = P_FEATURE,
  clusterAlg = CLUSTER_ALG,
  distance = DISTANCE,
  innerLinkage = INNER_LINK,
  finalLinkage = FINAL_LINK,
  seed = SEED,
  plot = "pdf",
  writeTable = FALSE,
  title = CCPDIR,
  verbose = FALSE
)
cat(sprintf(" CCP done. Native output -> %s\n\n", CCPDIR))

# 6. CDF + PAC
cat(">>> STEP 6: CDF / Delta-area / PAC\n")

# CDF: consensus matrix lower-triangle (excluding diagonal) values distribution
compute_cdf <- function(M) {
  v <- M[lower.tri(M, diag = FALSE)]
  ecdf_fn <- ecdf(v)
  data.frame(x = sort(v), y = ecdf_fn(sort(v)))
}

# PAC: P(consensus index in [0.1, 0.9])
compute_pac <- function(M, low = 0.1, high = 0.9) {
  v <- M[lower.tri(M, diag = FALSE)]
  e <- ecdf(v)
  e(high) - e(low)
}

ks <- 2:MAX_K
cdf_list <- lapply(ks, function(k) {
  d <- compute_cdf(ccp_res[[k]]$consensusMatrix)
  d$k <- k
  d
})
cdf_df <- do.call(rbind, cdf_list)

pac_vec <- vapply(ks, function(k) compute_pac(ccp_res[[k]]$consensusMatrix), numeric(1))
names(pac_vec) <- ks

# Delta area (Monti et al. 2003) -- area under CDF
auc_vec <- vapply(ks, function(k) {
  d <- compute_cdf(ccp_res[[k]]$consensusMatrix)
  # trapezoidal AUC
  sum(diff(c(0, d$x, 1)) * c(0, d$y))
}, numeric(1))
delta_auc <- c(auc_vec[1], diff(auc_vec) / auc_vec[-length(auc_vec)])
names(delta_auc) <- ks

cat(" PAC by k:\n"); print(round(pac_vec, 4))
cat(" Delta AUC by k:\n"); print(round(delta_auc, 4))

# --- choice of k: multi-criteria decision (balance + PAC local-min + parsimony) ---
# S1 cluster balance: min cluster >= 10% of N           [Senbabaoglu 2014, Sci Rep]
# S2 PAC local minimum (boundary-aware)                 [Senbabaoglu 2014, Sci Rep]
# Select smallest k satisfying S1 & S2; corroborated by monotone-decreasing
# relative CDF delta-area (diminishing returns of extra clusters; Monti 2003).
frac_vec <- vapply(ks, function(k) min(table(ccp_res[[k]]$consensusClass)) / N, numeric(1))
names(frac_vec) <- ks
is_pac_locmin <- vapply(seq_along(ks), function(i) {
  left  <- if (i == 1) TRUE else pac_vec[i] <= pac_vec[i - 1]
  right <- if (i == length(ks)) TRUE else pac_vec[i] <= pac_vec[i + 1]
  left && right
}, logical(1))
balanced_k   <- ks[frac_vec >= 0.10]
pac_locmin_k <- ks[is_pac_locmin]
admissible_k <- intersect(pac_locmin_k, balanced_k)
if (length(admissible_k) == 0L) stop("Choice of k: no admissible k (S1 & S2 intersection empty).")
K_FINAL <- min(admissible_k)
cat(" Min cluster fraction by k:\n"); print(round(frac_vec, 3))
cat(sprintf(" PAC local minima            : {%s}\n", paste(pac_locmin_k, collapse = ", ")))
cat(sprintf(" Balanced (min_frac >= 0.10) : {%s}\n", paste(balanced_k, collapse = ", ")))
cat(sprintf(" Admissible (S1 & S2)        : {%s}\n", paste(admissible_k, collapse = ", ")))
cat(sprintf(" Rel. delta-area monotone-decreasing: %s\n", all(diff(delta_auc) < 0)))
cat(sprintf(" Global-lowest-PAC k = %s (min_frac = %.3f)\n",
            names(which.min(pac_vec)), frac_vec[which.min(pac_vec)]))
if (K_FINAL != K_EXPECTED)
  warning(sprintf("Choice of k: data-driven k (%d) != expected k (%d); review inputs.", K_FINAL, K_EXPECTED))
cat(sprintf(" Choice of k -- final k = %d (expected %d)\n\n", K_FINAL, K_EXPECTED))

# 7. Fig 3A: CDF curves
cat(">>> STEP 7: Fig3A CDF curves\n")
p_cdf <- ggplot(cdf_df, aes(x = x, y = y, color = factor(k))) +
  geom_step(linewidth = 0.6) +
  scale_color_brewer(palette = "Dark2", name = "k") +
  labs(x = "Consensus index", y = "CDF") +
  theme_classic(base_size = 12)

cairo_pdf(file.path(FIGDIR, "Fig3A_CDF_curves.pdf"), width = 7, height = 5)
print(p_cdf); dev.off()
ragg::agg_png(file.path(FIGDIR, "Fig3A_CDF_curves.png"), width = 7, height = 5, units = "in", res = 600)
print(p_cdf); dev.off()

# 8. Fig 3B: PAC + Delta AUC
cat(">>> STEP 8: Fig3B PAC + Delta AUC\n")
metric_df <- data.frame(
  k = rep(ks, 2),
  value = c(pac_vec, delta_auc),
  metric = rep(c("PAC (lower = better)", "Delta AUC (elbow)"), each = length(ks))
)
p_pac <- ggplot(metric_df, aes(x = k, y = value)) +
  geom_line(color = "#1B9E77", linewidth = 0.8) +
  geom_point(color = "#1B9E77", size = 2.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = ks) +
  geom_vline(xintercept = K_FINAL, linetype = "dashed", color = "#D95F02") +
  labs(x = "k", y = NULL) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

cairo_pdf(file.path(FIGDIR, "Fig3B_PAC_score.pdf"), width = 9, height = 4)
print(p_pac); dev.off()
ragg::agg_png(file.path(FIGDIR, "Fig3B_PAC_score.png"), width = 9, height = 4, units = "in", res = 600)
print(p_pac); dev.off()

# 9. k=K_FINAL cluster assignment
cat(sprintf(">>> STEP 9: cluster assignment at k = %d\n", K_FINAL))
cluster_vec <- ccp_res[[K_FINAL]]$consensusClass
# standardized labels: C1, C2, ...
cluster_lbl <- paste0("C", as.integer(cluster_vec))
names(cluster_lbl) <- names(cluster_vec)

# Cluster size balance check
n_per <- table(cluster_lbl)
cat(" Cluster sizes:\n"); print(n_per)
min_frac <- min(n_per) / N
cat(sprintf(" Min cluster fraction: %.3f (criterion: >= 0.10)\n",
            min_frac))
if (min_frac < 0.10) {
  cat(" Cluster size imbalance below 10% -- review k choice.\n")
}

# Save Tab7 assignment
hist_map <- cli[, .(sampleID,
                    histology = fifelse(histological_type == "Esophagus Squamous Cell Carcinoma", "ESCC",
                                fifelse(histological_type == "Esophagus Adenocarcinoma, NOS", "EAC", "Other")))]
surv_sub <- surv[, .(sample, OS, OS.time)]

assign_dt <- data.table(
  sampleID = samples_keep,
  cluster = cluster_lbl[samples_keep]
)
assign_dt <- merge(assign_dt, hist_map, by = "sampleID", all.x = TRUE)
assign_dt <- merge(assign_dt, surv_sub, by.x = "sampleID", by.y = "sample", all.x = TRUE)
fwrite(assign_dt, file.path(TBLDIR, "Tab7_cluster_assignment.csv"))
cat(sprintf(" [OK] Tab7_cluster_assignment.csv\n\n"))

# 10. Fig 3C: consensus matrix heatmap @ k=K_FINAL
cat(">>> STEP 10: Fig3C consensus matrix heatmap\n")
cm <- ccp_res[[K_FINAL]]$consensusMatrix
# sample order given by CCP (within groups)
ord <- order(cluster_lbl[colnames(mat_z)])
cm_ord <- cm[ord, ord]
cluster_ord <- cluster_lbl[colnames(mat_z)][ord]

ha_top <- HeatmapAnnotation(
  Cluster = cluster_ord,
  col = list(Cluster = CLUSTER_COLORS),
  show_annotation_name = FALSE
)
ht_cm <- Heatmap(
  cm_ord,
  name = "Consensus",
  col = colorRamp2(c(0, 1), c("white", "#08519C")),
  top_annotation = ha_top,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = NULL,
  heatmap_legend_param = list(title = "Consensus")
)

cairo_pdf(file.path(FIGDIR, sprintf("Fig3C_consensus_matrix_k%d.pdf", K_FINAL)),
          width = 6.5, height = 6)
draw(ht_cm); dev.off()
ragg::agg_png(file.path(FIGDIR, sprintf("Fig3C_consensus_matrix_k%d.png", K_FINAL)),
              width = 6.5, height = 6, units = "in", res = 600)
draw(ht_cm); dev.off()

# 11. Fig 3D: PCA
cat(">>> STEP 11: Fig3D PCA on 18 FMG z-score\n")
pca <- prcomp(t(mat_z), center = FALSE, scale. = FALSE)
pca_df <- data.frame(
  sampleID = colnames(mat_z),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  cluster = cluster_lbl[colnames(mat_z)]
)
ve <- (pca$sdev^2) / sum(pca$sdev^2)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(linewidth = 0.5) +
  scale_color_manual(values = CLUSTER_COLORS) +
  labs(x = sprintf("PC1 (%.1f%%)", 100 * ve[1]),
       y = sprintf("PC2 (%.1f%%)", 100 * ve[2])) +
  theme_classic(base_size = 12)

cairo_pdf(file.path(FIGDIR, "Fig3D_PCA_clusters.pdf"), width = 6, height = 5)
print(p_pca); dev.off()
ragg::agg_png(file.path(FIGDIR, "Fig3D_PCA_clusters.png"), width = 6, height = 5, units = "in", res = 600)
print(p_pca); dev.off()

# 12. Fig 3E: KM survival
cat(">>> STEP 12: Fig3E KM survival (log-rank)\n")
km_dt <- assign_dt[!is.na(OS.time) & !is.na(OS)]
fit <- survfit(Surv(OS.time, OS) ~ cluster, data = km_dt)
sd_lr <- survdiff(Surv(OS.time, OS) ~ cluster, data = km_dt)
lr_p <- 1 - pchisq(sd_lr$chisq, length(sd_lr$n) - 1)
cat(sprintf(" log-rank p = %.4g (n = %d)\n", lr_p, nrow(km_dt)))

p_km <- ggsurvplot(
  fit,
  data = km_dt,
  pval = TRUE,
  pval.method = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = unname(CLUSTER_COLORS),
  legend.title= "Cluster",
  legend.labs = sort(unique(km_dt$cluster)),
  xlab = "Time (days)",
  ylab = "Overall survival",
  title = NULL,
  ggtheme = theme_classic(base_size = 12),
  risk.table.height = 0.28
)

cairo_pdf(file.path(FIGDIR, "Fig3E_KM_clusters.pdf"), width = 7, height = 6)
print(p_km, newpage = FALSE); dev.off()
ragg::agg_png(file.path(FIGDIR, "Fig3E_KM_clusters.png"), width = 7, height = 6, units = "in", res = 600)
print(p_km, newpage = FALSE); dev.off()

# 13. Fig 3F: 18 FMG heatmap by cluster
cat(">>> STEP 13: Fig3F 18 FMG heatmap (z-score, cluster + histology annotated)\n")
ord2 <- order(cluster_lbl[colnames(mat_z)])
mat_ord <- mat_z[, ord2]
cluster_ord2 <- cluster_lbl[colnames(mat_z)][ord2]
hist_ord2 <- assign_dt[match(colnames(mat_ord), sampleID), histology]

ha_top2 <- HeatmapAnnotation(
  Cluster = cluster_ord2,
  Histology = hist_ord2,
  col = list(Cluster = CLUSTER_COLORS,
             Histology = HIST_COLORS),
  annotation_legend_param = list(
    Cluster = list(title = "Cluster"),
    Histology = list(title = "Histology")
  )
)

ht_fmg <- Heatmap(
  mat_ord,
  name = "z-score",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  top_annotation = ha_top2,
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  column_split = factor(cluster_ord2, levels = sort(unique(cluster_ord2))),
  row_names_gp = gpar(fontsize = 10),
  column_title = NULL,
  heatmap_legend_param = list(title = "z-score")
)

cairo_pdf(file.path(FIGDIR, "Fig3F_FMG_heatmap_clusters.pdf"), width = 10, height = 6)
draw(ht_fmg, merge_legend = TRUE); dev.off()
ragg::agg_png(file.path(FIGDIR, "Fig3F_FMG_heatmap_clusters.png"), width = 10, height = 6, units = "in", res = 600)
draw(ht_fmg, merge_legend = TRUE); dev.off()

# 14. Cluster summary table
cat(">>> STEP 14: Tab7_cluster_summary.csv\n")
summary_dt <- assign_dt[, .(
  n = .N,
  n_ESCC = sum(histology == "ESCC", na.rm = TRUE),
  n_EAC = sum(histology == "EAC", na.rm = TRUE),
  n_Other = sum(histology == "Other",na.rm = TRUE),
  pct_ESCC = round(100 * mean(histology == "ESCC", na.rm = TRUE), 2),
  pct_EAC = round(100 * mean(histology == "EAC", na.rm = TRUE), 2),
  median_OS_days = as.double(median(OS.time, na.rm = TRUE)),
  death_rate = round(mean(OS, na.rm = TRUE), 3)
), by = cluster][order(cluster)]
summary_dt[, log_rank_p_overall := lr_p]
fwrite(summary_dt, file.path(TBLDIR, "Tab7_cluster_summary.csv"))
print(summary_dt)
cat("\n")

# 15. Histology x Cluster cross-tab (chi-square)
cat(">>> STEP 15: Histology x Cluster association (chi-square)\n")
ct <- table(assign_dt$cluster, assign_dt$histology)
print(ct)
chi <- suppressWarnings(chisq.test(ct))
cat(sprintf(" Chi-square p = %.4g\n\n", chi$p.value))

# 16. sessionInfo
sess_file <- file.path(LOGDIR, "sessionInfo/sessionInfo_stage07.txt")
sink(sess_file)
cat(sprintf("Session info captured: %s\n\n", Sys.time()))
print(sessionInfo())
sink()
cat(sprintf(">>> sessionInfo saved: %s\n", sess_file))

cat(sprintf("\n=== Stage 7 Pipeline done (%s) ===\n", Sys.time()))

# end
sink(type = "message")
sink()
close(log_con)

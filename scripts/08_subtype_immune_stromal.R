###############################################################################
## Stage 8: Subtype Characterization (Immune & Stromal)
##
## - 8A: ESTIMATE [PMID 24113773 -- Yoshihara 2013 Nat Commun]
## - 8B: CIBERSORT (LM22, perm=1000, QN=TRUE) [PMID 25822800 -- Newman 2015 Nat Methods]
## Inclusion criterion: CIBERSORT p < 0.05
## - 8C: ssGSEA + Bindea 28-cell-type signature [PMID 19847166 -- Barbie 2009 Nature
## PMID 24138885 -- Bindea 2013 Immunity]
## GSVA >=v1.50 ssgseaParam() API
## - 8D: Immune checkpoint genes (Liu 2024, n=30 after SLAMF4=CD244 alias)
## Wilcoxon rank-sum, BH FDR < 0.05
## - 8E: HLA family genes (Liu 2024, n=19)
## Wilcoxon rank-sum, BH FDR < 0.05
##
## Inputs:
## data/expr_ESCA_tumor.csv (TCGA log2(RSEM+1), 58581 x 181)
## tables/Tab7_cluster_assignment.csv (sampleID, cluster {C1,C2}, histology, OS)
## data/LM22.txt (CIBERSORT 22 cell type signature matrix)
## data/CIBERSORT.R (Newman 2015 v1.02 local script)
## data/Bindea_28_signatures.tsv (28 cell types x 578 unique genes)
##
## Outputs:
## Tab8A_ESTIMATE_scores.csv (per-sample 3 scores + cluster)
## Tab8A_ESTIMATE_stats.csv (per-score C1 vs C2 Wilcoxon)
## Tab8B_CIBERSORT_results.csv (per-sample 22 cell types + p-value)
## Tab8B_CIBERSORT_celltype_stats.csv (per-cell-type C1 vs C2 Wilcoxon)
## Tab8C_ssGSEA_Bindea.csv (per-sample 28 cell type score)
## Tab8C_ssGSEA_Bindea_celltype_stats.csv (per-cell-type C1 vs C2 Wilcoxon)
## Tab8D_checkpoint_wilcoxon.csv (per-gene C1 vs C2 + log2FC_C1_minus_C2)
## Tab8E_HLA_wilcoxon.csv (per-gene C1 vs C2 + log2FC_C1_minus_C2)
##
## Outputs (figures/, PDF cairo + PNG ragg @ 600 dpi):
## Fig4A_ESTIMATE_clusters.{pdf,png} (3-panel boxplot)
## Fig4B_CIBERSORT_clusters.{pdf,png} (22-panel per-cell-type boxplot)
## Fig4B_supp_CIBERSORT_stackedbar.{pdf,png} (per-sample stacked barplot)
## Fig4C_ssGSEA_Bindea_clusters.{pdf,png} (heatmap 28 x 181)
## Fig4D_immune_checkpoint_clusters.{pdf,png} (heatmap 30 x 181)
## Fig4E_HLA_clusters.{pdf,png} (heatmap 19 x 181)
##
## Random seed: 1262118388 (Liu et al. 2024)
##
## Notes:
## - Cluster colors (as in stage 7): C1 = "#1B9E77" (Dark2 teal), C2 = "#D95F02" (Dark2 orange)
## - Histology colors: ESCC = "#1B9E77", EAC = "#D95F02", Other = "#7570B3"
## - Heatmap z-score: row-wise (per cell type/gene)
## - Heatmap palette: RdBu reversed (low blue, high red), midpoint 0
##
## SLAMF4 -> CD244 alias mapped (NCBI Gene ID 51744). n_checkpoint = 30
##
## Bindea 28 source: hybrid (IOBR v2.2.0 27 + NK CD56bright direct from
## Bindea 2013 supp Table S1 [PMID 24138885]).
###############################################################################

# 0. + directory
# run from the repository root (see README.md)

INPUT_EXPR <- "data/expr_ESCA_tumor.csv"
INPUT_CLUSTER <- "tables/Tab7_cluster_assignment.csv"
INPUT_LM22 <- "data/LM22.txt"
INPUT_CIBSCR <- "data/CIBERSORT.R"
INPUT_BINDEA <- "data/Bindea_28_signatures.tsv"

TBLDIR <- "tables"
FIGDIR <- "figures"
LOGDIR <- "logs"
CACHE <- "cache"

for (d in c(TBLDIR, FIGDIR, LOGDIR, CACHE)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

ADJ_P <- 0.05 # Stage 8 BH FDR cutoff (8B/8C/8D/8E)
SEED <- 1262118388

# Cluster + histology palettes (as in stage 7)
CLUSTER_COL <- c(C1 = "#1B9E77", C2 = "#D95F02") # Dark2; same colours as stages 7 and 9
HIST_COL <- c(ESCC = "#1B9E77", EAC = "#D95F02", Other = "#7570B3")

set.seed(SEED)

# 0.6. Logging
log_file <- file.path(LOGDIR, "08_diagnostics.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat(sprintf("=== Stage 8 Pipeline start (%s) ===\n\n", Sys.time()))

# 1. load
suppressPackageStartupMessages({
  library(data.table)
  library(estimate)
  library(GSVA)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(ragg)
  library(grid)
  # CIBERSORT dependencies
  library(e1071)
  library(preprocessCore)
  library(parallel)
})

cat(sprintf("estimate : %s\n", packageVersion("estimate")))
cat(sprintf("GSVA : %s\n", packageVersion("GSVA")))
cat(sprintf("ComplexHeatmap : %s\n", packageVersion("ComplexHeatmap")))
cat(sprintf("ggplot2 / ragg : %s / %s\n",
            packageVersion("ggplot2"), packageVersion("ragg")))
cat("\n")

# 2. Load input
cat(">>> STEP 2: Load inputs\n")
expr_dt <- fread(INPUT_EXPR)
cluster_dt <- fread(INPUT_CLUSTER)
bindea_dt <- fread(INPUT_BINDEA)
cat(sprintf(" expr : %d genes x %d cols\n", nrow(expr_dt), ncol(expr_dt)))
cat(sprintf(" cluster: %d samples (%s)\n",
            nrow(cluster_dt),
            paste(sprintf("%s=%d", names(table(cluster_dt$cluster)),
                          as.integer(table(cluster_dt$cluster))),
                  collapse = ", ")))
cat(sprintf(" Bindea: %d (cell_type, gene) pairs, %d cell types, %d unique genes\n",
            nrow(bindea_dt),
            length(unique(bindea_dt$cell_type)),
            length(unique(bindea_dt$gene_symbol))))

# 3. Expression matrix: the 181 samples only
cat("\n>>> STEP 3: Subset expression to Tab7 181 samples\n")
key_g <- colnames(expr_dt)[1]
setnames(expr_dt, key_g, "gene_id")

samples_keep <- cluster_dt$sampleID
sample_in_expr <- intersect(samples_keep, colnames(expr_dt))
cat(sprintf(" Tab7 samples in expr matrix: %d / %d\n",
            length(sample_in_expr), length(samples_keep)))
stopifnot(length(sample_in_expr) == nrow(cluster_dt))

# Build matrix: gene rows, sample cols
expr_mat <- as.matrix(expr_dt[, ..sample_in_expr])
rownames(expr_mat) <- expr_dt$gene_id

# Cluster vector aligned to expr columns
cluster_dt <- cluster_dt[match(colnames(expr_mat), sampleID)]
stopifnot(all(cluster_dt$sampleID == colnames(expr_mat)))
cluster_vec <- factor(cluster_dt$cluster, levels = c("C1","C2"))
hist_vec <- factor(cluster_dt$histology, levels = c("ESCC","EAC","Other"))
N_TOTAL <- length(cluster_vec)
cat(sprintf(" N_TOTAL (analysis): %d samples (C1=%d, C2=%d)\n",
            N_TOTAL, sum(cluster_vec=="C1"), sum(cluster_vec=="C2")))

# helper: Wilcoxon C1 vs C2 + BH
wilcox_c1c2 <- function(values, cluster) {
  # returns p-value of Wilcoxon rank-sum (two-sided)
  v1 <- values[cluster == "C1"]
  v2 <- values[cluster == "C2"]
  if (length(v1) < 2 || length(v2) < 2) return(NA_real_)
  if (all(is.na(v1)) || all(is.na(v2))) return(NA_real_)
  tryCatch(
    wilcox.test(v1, v2, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}
sig_label <- function(p) {
  ifelse(is.na(p), "ns",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01, "**",
  ifelse(p < 0.05, "*", "ns"))))
}

# helper: save an output figure
save_fig <- function(plot_obj, stem, w = 12, h = 5, is_heatmap = FALSE) {
  pdf_f <- file.path(FIGDIR, paste0(stem, ".pdf"))
  png_f <- file.path(FIGDIR, paste0(stem, ".png"))
  if (is_heatmap) {
    cairo_pdf(pdf_f, width = w, height = h)
    draw(plot_obj, merge_legend = TRUE)
    dev.off()
    ragg::agg_png(png_f, width = w, height = h, units = "in", res = 600)
    draw(plot_obj, merge_legend = TRUE)
    dev.off()
  } else {
    ggsave(pdf_f, plot = plot_obj, width = w, height = h, device = cairo_pdf)
    ggsave(png_f, plot = plot_obj, width = w, height = h, dpi = 600,
           device = ragg::agg_png)
  }
  cat(sprintf(" [OK] %s.{pdf,png}\n", stem))
}

# 4. STEP A: ESTIMATE (8A)
cat("\n>>> STEP A: ESTIMATE\n")

# ESTIMATE workflow: write expr to GCT, run filterCommonGenes -> estimateScore
# Use GCT format (rows = gene, cols = sample), platform = "illumina" for RNA-seq
expr_for_est <- file.path(CACHE, "expr_for_estimate.txt")
filt_for_est <- file.path(CACHE, "expr_filtered_estimate.gct")
out_for_est <- file.path(CACHE, "estimate_scores.gct")

# tab-separated, first column = gene symbol (col name = "GeneSymbol")
expr_for_est_dt <- data.table(GeneSymbol = rownames(expr_mat))
expr_for_est_dt <- cbind(expr_for_est_dt, as.data.table(expr_mat))
fwrite(expr_for_est_dt, expr_for_est, sep = "\t", quote = FALSE)

estimate::filterCommonGenes(input.f = expr_for_est,
                             output.f = filt_for_est,
                             id = "GeneSymbol")
estimate::estimateScore(input.ds = filt_for_est,
                         output.ds = out_for_est,
                         platform = "illumina") # RNA-seq -> illumina

# Read ESTIMATE GCT result
gct_lines <- readLines(out_for_est)
# GCT: line 1 = "#1.2", line 2 = N_rows N_cols, line 3 = header, line 4+ = data
hdr <- strsplit(gct_lines[3], "\t")[[1]]
data_lines <- gct_lines[-(1:3)]
data_split <- do.call(rbind, strsplit(data_lines, "\t"))
score_dt <- as.data.table(data_split)
setnames(score_dt, hdr)
# First two cols = NAME, Description; rest = sample cols
sample_cols_est <- setdiff(hdr, c("NAME","Description"))
# Transpose: long format
estimate_long <- score_dt[, .(score = NAME, value = lapply(.SD, as.numeric)),
                          .SDcols = sample_cols_est]
# Build wide -> long manually
score_mat <- as.matrix(score_dt[, ..sample_cols_est])
storage.mode(score_mat) <- "numeric"
rownames(score_mat) <- score_dt$NAME

# Match samples to cluster
match_idx <- match(colnames(expr_mat), sample_cols_est)
# Sometimes ESTIMATE alters sample names (replaces - with .) -- handle both
if (any(is.na(match_idx))) {
  alt_cols <- gsub("\\.", "-", sample_cols_est)
  match_idx <- match(colnames(expr_mat), alt_cols)
  colnames(score_mat) <- alt_cols
}
stopifnot(!any(is.na(match_idx)))
score_mat <- score_mat[, match_idx]

# Output table: sampleID, cluster, histology, StromalScore, ImmuneScore, ESTIMATEScore
estimate_tbl <- data.table(
  sampleID = colnames(expr_mat),
  cluster = as.character(cluster_vec),
  histology = as.character(hist_vec),
  StromalScore = score_mat["StromalScore", ],
  ImmuneScore = score_mat["ImmuneScore", ],
  ESTIMATEScore = score_mat["ESTIMATEScore", ]
)
fwrite(estimate_tbl, file.path(TBLDIR, "Tab8A_ESTIMATE_scores.csv"))
cat(sprintf(" [OK] Tab8A_ESTIMATE_scores.csv (%d samples)\n", nrow(estimate_tbl)))

# Wilcoxon C1 vs C2 per score
build_estimate_stat <- function(score_name, vals, cluster) {
  data.table(
    score = score_name,
    n_C1 = sum(cluster == "C1" & !is.na(vals)),
    n_C2 = sum(cluster == "C2" & !is.na(vals)),
    median_C1 = median(vals[cluster == "C1"], na.rm = TRUE),
    median_C2 = median(vals[cluster == "C2"], na.rm = TRUE),
    mean_C1 = mean( vals[cluster == "C1"], na.rm = TRUE),
    mean_C2 = mean( vals[cluster == "C2"], na.rm = TRUE),
    p_value = wilcox_c1c2(vals, cluster)
  )
}
estimate_stats <- rbindlist(list(
  build_estimate_stat("StromalScore", estimate_tbl$StromalScore, estimate_tbl$cluster),
  build_estimate_stat("ImmuneScore", estimate_tbl$ImmuneScore, estimate_tbl$cluster),
  build_estimate_stat("ESTIMATEScore", estimate_tbl$ESTIMATEScore, estimate_tbl$cluster)
))
estimate_stats[, p_adj_BH := p.adjust(p_value, method = "BH")]
estimate_stats[, sig_label := sig_label(p_adj_BH)]
estimate_stats[, significant_BH_05 := !is.na(p_adj_BH) & p_adj_BH < ADJ_P]

# Backward-compat: keep simple 'sig' column for downstream plotting code
estimate_stats[, sig := sig_label]

fwrite(estimate_stats, file.path(TBLDIR, "Tab8A_ESTIMATE_stats.csv"))
cat(sprintf(" [OK] Tab8A_ESTIMATE_stats.csv\n"))
cat(" ESTIMATE C1 vs C2 stats:\n")
print(estimate_stats)

# Fig4A: 3-panel boxplot
plot_dt_a <- melt(estimate_tbl,
                  id.vars = c("sampleID","cluster","histology"),
                  measure.vars = c("StromalScore","ImmuneScore","ESTIMATEScore"),
                  variable.name = "score", value.name = "value")
plot_dt_a[, score := factor(score, levels = c("StromalScore","ImmuneScore","ESTIMATEScore"))]
plot_dt_a[, cluster := factor(cluster, levels = c("C1","C2"))]

estimate_stats_plot <- copy(estimate_stats)
estimate_stats_plot[, score := factor(score, levels = c("StromalScore","ImmuneScore","ESTIMATEScore"))]
ymax_a <- plot_dt_a[, .(y_max = max(value, na.rm = TRUE)), by = score]
estimate_stats_plot <- merge(estimate_stats_plot, ymax_a, by = "score")
estimate_stats_plot[, y_pos := y_max + 0.05 * abs(y_max)]

p_a <- ggplot(plot_dt_a, aes(x = cluster, y = value, fill = cluster)) +
  geom_violin(width = 0.8, alpha = 0.6, color = "black", linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.4, color = "black", linewidth = 0.3) +
  geom_text(data = estimate_stats_plot,
            aes(x = 1.5, y = y_pos, label = sig),
            inherit.aes = FALSE, size = 4.5, fontface = "bold") +
  facet_wrap(~ score, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = CLUSTER_COL, name = "Cluster") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = "ESTIMATE score") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face="bold", size=10),
        axis.text.x = element_text(face="bold"),
        legend.position = "top")

save_fig(p_a, "Fig4A_ESTIMATE_clusters", w = 10, h = 5)

# 5. STEP B: CIBERSORT (8B)
cat("\n>>> STEP B: CIBERSORT\n")
cat(" Writing mixture file (TCGA log2(RSEM+1) directly, QN=TRUE)...\n")

mixture_file <- file.path(CACHE, "tcga_mixture_for_cibersort.txt")
mix_dt <- data.table(`Gene symbol` = rownames(expr_mat))
mix_dt <- cbind(mix_dt, as.data.table(expr_mat))
fwrite(mix_dt, mixture_file, sep = "\t", quote = FALSE)

# Source CIBERSORT and run
cat(sprintf(" Sourcing %s ...\n", INPUT_CIBSCR))
source(INPUT_CIBSCR) # defines CIBERSORT() function

cat(" Running CIBERSORT(perm=1000, QN=TRUE) ... (~5-15 min single-threaded)\n")
cs_t0 <- Sys.time()
ciber_res <- CIBERSORT(INPUT_LM22, mixture_file, perm = 1000, QN = TRUE)
cs_dt <- Sys.time() - cs_t0
cat(sprintf(" CIBERSORT done in %.1f min\n", as.numeric(cs_dt, units = "mins")))

# ciber_res = matrix, rows = samples, cols = 22 cell types + P-value + Correlation + RMSE
ciber_dt <- as.data.table(ciber_res, keep.rownames = "sampleID")
ciber_dt[, cluster := cluster_dt$cluster[match(sampleID, cluster_dt$sampleID)]]
ciber_dt[, histology := cluster_dt$histology[match(sampleID, cluster_dt$sampleID)]]
fwrite(ciber_dt, file.path(TBLDIR, "Tab8B_CIBERSORT_results.csv"))
cat(sprintf(" [OK] Tab8B_CIBERSORT_results.csv (%d samples)\n", nrow(ciber_dt)))

# Filter samples with CIBERSORT p < 0.05
p_col <- grep("^P-value$|^p-value$|^Pvalue$|^P\\.value$", colnames(ciber_dt),
              value = TRUE, ignore.case = TRUE)
if (length(p_col) == 0) p_col <- "P-value"
ciber_dt[, ciber_p := as.numeric(get(p_col))]
n_pass <- sum(ciber_dt$ciber_p < 0.05, na.rm = TRUE)
cat(sprintf(" CIBERSORT p<0.05 inclusion: %d / %d samples\n",
            n_pass, nrow(ciber_dt)))

ciber_pass <- ciber_dt[ciber_p < 0.05]

# 22 cell type columns
ct22 <- colnames(ciber_res)
ct22 <- setdiff(ct22, c("P-value","Correlation","RMSE","Pearson Correlation","Pearson_Correlation"))
cat(sprintf(" CIBERSORT cell types: %d\n", length(ct22)))

# Per-cell-type Wilcoxon C1 vs C2 (only on pass samples)
b_stats <- rbindlist(lapply(ct22, function(ct) {
  vals <- as.numeric(ciber_pass[[ct]])
  data.table(
    cell_type = ct,
    n_C1 = sum(ciber_pass$cluster == "C1"),
    n_C2 = sum(ciber_pass$cluster == "C2"),
    median_C1 = median(vals[ciber_pass$cluster == "C1"], na.rm = TRUE),
    median_C2 = median(vals[ciber_pass$cluster == "C2"], na.rm = TRUE),
    p_value = wilcox_c1c2(vals, ciber_pass$cluster)
  )
}))
b_stats[, p_adj_BH := p.adjust(p_value, method = "BH")]
b_stats[, sig := sig_label(p_adj_BH)]
b_stats[, significant_BH_05 := !is.na(p_adj_BH) & p_adj_BH < ADJ_P]
fwrite(b_stats, file.path(TBLDIR, "Tab8B_CIBERSORT_results.csv"),
       append = FALSE) # overwrite with extended file
# Actually we want both: the per-sample matrix (already saved) AND the per-cell-type stats
# Save stats separately to avoid clobbering
fwrite(b_stats, file.path(TBLDIR, "Tab8B_CIBERSORT_celltype_stats.csv"))
# Re-save per-sample
fwrite(ciber_dt, file.path(TBLDIR, "Tab8B_CIBERSORT_results.csv"))
cat(sprintf(" [OK] Tab8B_CIBERSORT_celltype_stats.csv\n"))
cat(sprintf(" Significant cell types (BH<%.2f): %d / %d\n",
            ADJ_P, sum(b_stats$significant_BH_05), nrow(b_stats)))

# Fig4B: per-cell-type boxplot, 22 panels (4x6 layout)
long_b <- melt(ciber_pass[, c("sampleID","cluster", ct22), with = FALSE],
               id.vars = c("sampleID","cluster"),
               variable.name = "cell_type", value.name = "proportion")
long_b[, cell_type := factor(cell_type, levels = ct22)]
long_b[, cluster := factor(cluster, levels = c("C1","C2"))]
long_b[, proportion := as.numeric(proportion)]

ymax_b <- long_b[, .(y_max = max(proportion, na.rm = TRUE)), by = cell_type]
b_stats_plot <- copy(b_stats)
b_stats_plot[, cell_type := factor(cell_type, levels = ct22)]
b_stats_plot <- merge(b_stats_plot, ymax_b, by = "cell_type")
b_stats_plot[, y_pos := y_max + 0.05 * abs(y_max) + 0.005]

p_b <- ggplot(long_b, aes(x = cluster, y = proportion, fill = cluster)) +
  geom_boxplot(outlier.size = 0.4, color = "black", linewidth = 0.3) +
  geom_text(data = b_stats_plot,
            aes(x = 1.5, y = y_pos, label = sig),
            inherit.aes = FALSE, size = 3, fontface = "bold") +
  facet_wrap(~ cell_type, ncol = 6, scales = "free_y") +
  scale_fill_manual(values = CLUSTER_COL, name = "Cluster") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.18)),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Cell-type proportion") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face="bold", size=8),
        axis.text.x = element_text(face="bold", size=8),
        axis.text.y = element_text(size = 7),
        legend.position = "top")

save_fig(p_b, "Fig4B_CIBERSORT_clusters", w = 14, h = 9)

# Fig4B_supp: per-sample stacked barplot (col by cluster order)
# Order samples by cluster then by ESTIMATE Stromal Score (or similar)
stack_df <- copy(long_b)
sample_order <- ciber_pass[order(cluster, sampleID), sampleID]
stack_df[, sampleID := factor(sampleID, levels = sample_order)]

# Color palette for 22 cell types: ColorBrewer Set3 + Paired (combined for 22)
ct_colors <- c(
  RColorBrewer::brewer.pal(12, "Set3"),
  RColorBrewer::brewer.pal(10, "Paired")
)[seq_len(length(ct22))]
names(ct_colors) <- ct22

cluster_strip <- ciber_pass[order(cluster, sampleID), .(sampleID, cluster)]
p_b_supp <- ggplot(stack_df, aes(x = sampleID, y = proportion, fill = cell_type)) +
  geom_col(width = 1) +
  scale_fill_manual(values = ct_colors, name = "Cell type") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0,0))) +
  labs(x = NULL, y = "Proportion") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")) +
  guides(fill = guide_legend(ncol = 1))

save_fig(p_b_supp, "Fig4B_supp_CIBERSORT_stackedbar", w = 14, h = 6)

# 6. STEP C: ssGSEA Bindea (8C)
cat("\n>>> STEP C: ssGSEA Bindea 28-cell-type\n")

# Build gene sets (named list of character vectors)
gene_sets <- split(bindea_dt$gene_symbol, bindea_dt$cell_type)
# Filter to genes present in expr matrix
gene_sets_filt <- lapply(gene_sets, function(g) intersect(g, rownames(expr_mat)))
n_per_set <- vapply(gene_sets_filt, length, integer(1))
cat(sprintf(" Bindea gene sets: %d cell types, gene counts after expr filter:\n",
            length(gene_sets_filt)))
print(sort(n_per_set))

# GSVA 2.0+ ssgseaParam API
ssgsea_params <- ssgseaParam(exprData = expr_mat, geneSets = gene_sets_filt)
cat(" Running ssGSEA via GSVA::gsva(ssgseaParam(...)) ...\n")
ssgsea_t0 <- Sys.time()
ssgsea_res <- gsva(ssgsea_params, verbose = FALSE)
cat(sprintf(" ssGSEA done in %.1f s\n",
            as.numeric(Sys.time() - ssgsea_t0, units = "secs")))

# ssgsea_res: matrix [cell_type x sample]
cat(sprintf(" ssGSEA result: %d cell types x %d samples\n",
            nrow(ssgsea_res), ncol(ssgsea_res)))

# Per-cell-type Wilcoxon
c_stats <- rbindlist(lapply(rownames(ssgsea_res), function(ct) {
  vals <- ssgsea_res[ct, ]
  data.table(
    cell_type = ct,
    median_C1 = median(vals[cluster_vec == "C1"]),
    median_C2 = median(vals[cluster_vec == "C2"]),
    p_value = wilcox_c1c2(vals, as.character(cluster_vec))
  )
}))
c_stats[, p_adj_BH := p.adjust(p_value, method = "BH")]
c_stats[, sig := sig_label(p_adj_BH)]
c_stats[, significant_BH_05 := !is.na(p_adj_BH) & p_adj_BH < ADJ_P]

# Save ssGSEA: cell_type rows, sample cols + cluster + histology header at top via separate file
ssgsea_save <- as.data.table(ssgsea_res, keep.rownames = "cell_type")
fwrite(ssgsea_save, file.path(TBLDIR, "Tab8C_ssGSEA_Bindea.csv"))
fwrite(c_stats, file.path(TBLDIR, "Tab8C_ssGSEA_Bindea_celltype_stats.csv"))
cat(sprintf(" [OK] Tab8C_ssGSEA_Bindea.csv + Tab8C_ssGSEA_Bindea_celltype_stats.csv\n"))
cat(sprintf(" Significant cell types (BH<%.2f): %d / %d\n",
            ADJ_P, sum(c_stats$significant_BH_05), nrow(c_stats)))

# Fig4C: Heatmap (row z-score, columns ordered by cluster then histology)
col_order <- order(cluster_vec, hist_vec)
ssg_z <- t(scale(t(ssgsea_res[, col_order]), center = TRUE, scale = TRUE))

col_anno_c <- HeatmapAnnotation(
  Cluster = cluster_vec[col_order],
  Histology = hist_vec[col_order],
  col = list(Cluster = CLUSTER_COL, Histology = HIST_COL),
  annotation_name_side = "left",
  annotation_legend_param = list(
    Cluster = list(title = "Cluster"),
    Histology = list(title = "Histology")
  )
)

# Significance row annotation
sig_row <- c_stats$sig[match(rownames(ssg_z), c_stats$cell_type)]
row_anno_c <- rowAnnotation(
  BH_FDR = anno_text(sig_row, gp = gpar(fontsize = 8, fontface = "bold"),
                     just = "left", location = unit(2, "mm")),
  show_annotation_name = FALSE,
  width = unit(8, "mm")
)

heat_col <- circlize::colorRamp2(c(-2, 0, 2),
                                  c("#0571B0", "white", "#CA0020"))

ht_c <- Heatmap(ssg_z,
                name = "z-score",
                col = heat_col,
                top_annotation = col_anno_c,
                right_annotation = row_anno_c,
                cluster_columns = FALSE,
                cluster_rows = TRUE,
                show_column_names = FALSE,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 9),
                column_title = NULL,
                heatmap_legend_param = list(title = "Row z-score"))

save_fig(ht_c, "Fig4C_ssGSEA_Bindea_clusters", w = 14, h = 8, is_heatmap = TRUE)

# 7. STEP D: Immune checkpoint (8D)
cat("\n>>> STEP D: Immune checkpoint genes (n=30 after alias)\n")

# Liu 2024 checkpoint list (alias mapped per NCBI Gene; SLAMF4 = CD244)
checkpoint_alias <- c(
  "PD-L1" = "CD274",
  "SLAMF4" = "CD244",
  "CD96" = "CD96",
  "CSF1R" = "CSF1R",
  "TIM3" = "HAVCR2",
  "IL10RB" = "IL10RB",
  "KIR2DL1" = "KIR2DL1",
  "KIR2DL3" = "KIR2DL3",
  "LGALS9" = "LGALS9",
  "PD-1" = "PDCD1",
  "PD-L2" = "PDCD1LG2",
  "TGFB1" = "TGFB1",
  "B7-H3" = "CD276",
  "CD48" = "CD48",
  "CD70" = "CD70",
  "B7.2" = "CD86",
  "CXCR4" = "CXCR4",
  "IL6R" = "IL6R",
  "MICB" = "MICB",
  "TMIGD2" = "TMIGD2",
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
# Note: 4-1BB and 4-1BB-L both -> TNFRSF9 / TNFSF9. Liu 2024 may treat as separate
# (receptor + ligand). But TNFSF9 = 4-1BB-L, TNFRSF9 = 4-1BB. Already correct.

stopifnot(length(checkpoint_alias) == 30)
cp_present <- checkpoint_alias %in% rownames(expr_mat)
cat(sprintf(" Checkpoint genes resolved in expr matrix: %d / %d\n",
            sum(cp_present), length(checkpoint_alias)))
miss <- names(checkpoint_alias)[!cp_present]
if (length(miss) > 0) {
  cat(sprintf(" WARNING: missing: %s\n", paste(miss, collapse = ", ")))
}

cp_genes <- checkpoint_alias[cp_present]
cp_mat <- expr_mat[as.character(cp_genes), ]
rownames(cp_mat) <- names(cp_genes) # use Liu 2024 display names

# Wilcoxon C1 vs C2 per gene
d_stats <- rbindlist(lapply(rownames(cp_mat), function(g) {
  vals <- cp_mat[g, ]
  data.table(
    gene_display = g,
    gene_symbol = as.character(checkpoint_alias[g]),
    median_C1 = median(vals[cluster_vec == "C1"]),
    median_C2 = median(vals[cluster_vec == "C2"]),
    log2FC_C1_minus_C2 = median(vals[cluster_vec == "C1"]) -
                             median(vals[cluster_vec == "C2"]),
    p_value = wilcox_c1c2(vals, as.character(cluster_vec))
  )
}))
d_stats[, p_adj_BH := p.adjust(p_value, method = "BH")]
d_stats[, sig := sig_label(p_adj_BH)]
d_stats[, significant_BH_05 := !is.na(p_adj_BH) & p_adj_BH < ADJ_P]
fwrite(d_stats, file.path(TBLDIR, "Tab8D_checkpoint_wilcoxon.csv"))
cat(sprintf(" [OK] Tab8D_checkpoint_wilcoxon.csv\n"))
cat(sprintf(" Significant checkpoint genes (BH<%.2f): %d / %d\n",
            ADJ_P, sum(d_stats$significant_BH_05), nrow(d_stats)))

# Fig4D heatmap
cp_z <- t(scale(t(cp_mat[, col_order]), center = TRUE, scale = TRUE))
sig_row_d <- d_stats$sig[match(rownames(cp_z), d_stats$gene_display)]
row_anno_d <- rowAnnotation(
  BH_FDR = anno_text(sig_row_d, gp = gpar(fontsize = 8, fontface = "bold"),
                     just = "left", location = unit(2, "mm")),
  show_annotation_name = FALSE,
  width = unit(8, "mm")
)
ht_d <- Heatmap(cp_z,
                name = "z-score",
                col = heat_col,
                top_annotation = col_anno_c,
                right_annotation = row_anno_d,
                cluster_columns = FALSE,
                cluster_rows = TRUE,
                show_column_names = FALSE,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 9),
                column_title = NULL,
                heatmap_legend_param = list(title = "Row z-score"))
save_fig(ht_d, "Fig4D_immune_checkpoint_clusters", w = 14, h = 7, is_heatmap = TRUE)

# 8. STEP E: HLA family (8E)
cat("\n>>> STEP E: HLA family genes (n=19)\n")

hla_list <- c(
  "HLA-A","HLA-B","HLA-C",
  "HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB",
  "HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2",
  "HLA-DRA","HLA-DRB1","HLA-DRB5",
  "HLA-E","HLA-F","HLA-G"
)
hla_present <- hla_list %in% rownames(expr_mat)
cat(sprintf(" HLA genes resolved: %d / %d\n", sum(hla_present), length(hla_list)))
hla_genes <- hla_list[hla_present]
hla_mat <- expr_mat[hla_genes, ]

e_stats <- rbindlist(lapply(rownames(hla_mat), function(g) {
  vals <- hla_mat[g, ]
  data.table(
    gene_symbol = g,
    median_C1 = median(vals[cluster_vec == "C1"]),
    median_C2 = median(vals[cluster_vec == "C2"]),
    log2FC_C1_minus_C2 = median(vals[cluster_vec == "C1"]) -
                          median(vals[cluster_vec == "C2"]),
    p_value = wilcox_c1c2(vals, as.character(cluster_vec))
  )
}))
e_stats[, p_adj_BH := p.adjust(p_value, method = "BH")]
e_stats[, sig := sig_label(p_adj_BH)]
e_stats[, significant_BH_05 := !is.na(p_adj_BH) & p_adj_BH < ADJ_P]
fwrite(e_stats, file.path(TBLDIR, "Tab8E_HLA_wilcoxon.csv"))
cat(sprintf(" [OK] Tab8E_HLA_wilcoxon.csv\n"))
cat(sprintf(" Significant HLA genes (BH<%.2f): %d / %d\n",
            ADJ_P, sum(e_stats$significant_BH_05), nrow(e_stats)))

# Fig4E heatmap
hla_z <- t(scale(t(hla_mat[, col_order]), center = TRUE, scale = TRUE))
sig_row_e <- e_stats$sig[match(rownames(hla_z), e_stats$gene_symbol)]
row_anno_e <- rowAnnotation(
  BH_FDR = anno_text(sig_row_e, gp = gpar(fontsize = 8, fontface = "bold"),
                     just = "left", location = unit(2, "mm")),
  show_annotation_name = FALSE,
  width = unit(8, "mm")
)
ht_e <- Heatmap(hla_z,
                name = "z-score",
                col = heat_col,
                top_annotation = col_anno_c,
                right_annotation = row_anno_e,
                cluster_columns = FALSE,
                cluster_rows = TRUE,
                show_column_names = FALSE,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 9),
                column_title = NULL,
                heatmap_legend_param = list(title = "Row z-score"))
save_fig(ht_e, "Fig4E_HLA_clusters", w = 14, h = 6, is_heatmap = TRUE)

# 9. Summary Report
cat("\n========================================================================\n")
cat(" Stage 8 Subtype Characterization Summary Report \n")
cat("========================================================================\n")
cat(sprintf(" N_TOTAL : %d (C1=%d, C2=%d)\n",
            N_TOTAL, sum(cluster_vec=="C1"), sum(cluster_vec=="C2")))
cat(sprintf(" BH FDR cutoff : %.3f\n\n", ADJ_P))

cat(" --- 8A ESTIMATE C1 vs C2 ---\n")
print(estimate_stats)

cat(sprintf("\n --- 8B CIBERSORT ---\n"))
cat(sprintf(" Samples passing p<0.05 : %d / %d\n",
            sum(ciber_dt$ciber_p < 0.05, na.rm=TRUE), nrow(ciber_dt)))
cat(sprintf(" Significant cell types : %d / %d\n",
            sum(b_stats$significant_BH_05), nrow(b_stats)))
cat(" Top 5 by smallest p_adj_BH:\n")
print(head(b_stats[order(p_adj_BH), .(cell_type, median_C1, median_C2, p_adj_BH, sig)], 5))

cat(sprintf("\n --- 8C ssGSEA Bindea 28 ---\n"))
cat(sprintf(" Significant cell types : %d / %d\n",
            sum(c_stats$significant_BH_05), nrow(c_stats)))
cat(" Top 5 by smallest p_adj_BH:\n")
print(head(c_stats[order(p_adj_BH), .(cell_type, median_C1, median_C2, p_adj_BH, sig)], 5))

cat(sprintf("\n --- 8D Immune checkpoint (n=%d) ---\n", nrow(d_stats)))
cat(sprintf(" Significant genes (BH<%.2f) : %d / %d\n",
            ADJ_P, sum(d_stats$significant_BH_05), nrow(d_stats)))
cat(" Top 5 by abs(log2FC_C1_minus_C2) among sig:\n")
top_d <- d_stats[significant_BH_05 == TRUE]
top_d[, abs_lfc := abs(log2FC_C1_minus_C2)]
setorder(top_d, -abs_lfc)
print(head(top_d[, .(gene_display, gene_symbol, log2FC_C1_minus_C2, p_adj_BH, sig)], 5))

cat(sprintf("\n --- 8E HLA family (n=%d) ---\n", nrow(e_stats)))
cat(sprintf(" Significant genes (BH<%.2f) : %d / %d\n",
            ADJ_P, sum(e_stats$significant_BH_05), nrow(e_stats)))
cat(" All HLA genes sorted by p_adj_BH:\n")
print(e_stats[order(p_adj_BH), .(gene_symbol, median_C1, median_C2, log2FC_C1_minus_C2, p_adj_BH, sig)])

cat("\n========================================================================\n")

# 10. sessionInfo
sess_file <- file.path(LOGDIR, "sessionInfo/sessionInfo_stage08.txt")
sink(sess_file)
cat(sprintf("Session info captured: %s\n\n", Sys.time()))
print(sessionInfo())
sink()
cat(sprintf("\n>>> sessionInfo saved: %s\n", sess_file))

cat(sprintf("\n=== Stage 8 Pipeline done (%s) ===\n", Sys.time()))

sink(type = "message")
sink()
close(log_con)

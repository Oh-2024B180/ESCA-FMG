###############################################################################
##  Stage 8 / Step B (CIBERSORT) — RERUN with linear-space input + QN=FALSE
##
##  Reason for rerun:
##    Previous run (log2(RSEM+1) input + QN=TRUE) showed deconvolution failure:
##      mean Correlation ≈ 0.024 (normal ≥ 0.5)
##      Only 3 / 181 samples passed ciber_p < 0.05 → cell-type proportions
##      essentially noise.
##    Cause: log2(RSEM+1) + QN=TRUE mismatches LM22 reference space
##           (LM22 was built on Affymetrix HGU133+2.0 microarray, linear scale).
##    Fix (per Newman et al. 2015 Nat Methods + CIBERSORT FAQ for RNA-seq):
##      use linear-space input + QN=FALSE.
##
##  CRITICAL CHANGES vs original Step B (exactly 2):
##    1. Mixture matrix: log2(RSEM+1) -> linear via 2^x - 1, floor negatives to 0
##    2. CIBERSORT option: QN = TRUE -> QN = FALSE
##    perm = 1000  (unchanged)
##
##  Inputs / Outputs are exposed as variables at the top of this script so that
##  the user can adjust paths without editing the body.
##
##  Outputs (overwrite previous Step B):
##    Tab8B_CIBERSORT_results.csv          (per-sample 22 cell types + p/corr/RMSE + cluster + histology)
##    Tab8B_CIBERSORT_celltype_stats.csv   (per-cell-type C1 vs C2 Wilcoxon)
##    Fig4B_CIBERSORT_clusters.{pdf,png}   (per-cell-type boxplot, 6-col facet)
##    Fig4B_supp_CIBERSORT_stackedbar.{pdf,png}  (per-sample stacked bar, cluster-faceted)
##    08B_rerun_diagnostics.txt            (Correlation/p/RMSE summary + verdict)
##
##  NOT touched (downstream of original Step A/C/D/E):
##    Tab8A, Tab8C, Tab8D, Tab8E and Fig4A / Fig4C / Fig4D / Fig4E
##

# ── 1. PATHS / VARIABLES (user-editable) ─────────────────────────────────
# Working directory
PROJECT_ROOT <- getwd()   # run from the repository root (see README.md)

# Inputs
EXPR_FILE      <- "data/expr_ESCA_tumor.csv"            # rows=genes, cols=samples; first col = gene symbol
CLUSTER_FILE   <- "tables/Tab7_cluster_assignment.csv"
LM22_FILE      <- "data/LM22.txt"                       # CIBERSORT signature matrix
CIBERSORT_R    <- "data/CIBERSORT.R"                    # Newman lab v1.02 local script

# Cluster file column names (Tab7 schema)
CLUSTER_COL_SAMPLE  <- "sampleID"
CLUSTER_COL_CLUSTER <- "cluster"
CLUSTER_COL_HIST    <- "histology"

# Outputs
OUT_TBL_DIR    <- "tables"
OUT_FIG_DIR    <- "figures"
OUT_LOG_DIR    <- "logs"
CACHE_DIR      <- "cache"

TAB_RESULTS    <- file.path(OUT_TBL_DIR, "Tab8B_CIBERSORT_results.csv")
TAB_CELLSTAT   <- file.path(OUT_TBL_DIR, "Tab8B_CIBERSORT_celltype_stats.csv")
FIG_BOX_PDF    <- file.path(OUT_FIG_DIR, "Fig4B_CIBERSORT_clusters.pdf")
FIG_BOX_PNG    <- file.path(OUT_FIG_DIR, "Fig4B_CIBERSORT_clusters.png")
FIG_STACK_PDF  <- file.path(OUT_FIG_DIR, "Fig4B_supp_CIBERSORT_stackedbar.pdf")
FIG_STACK_PNG  <- file.path(OUT_FIG_DIR, "Fig4B_supp_CIBERSORT_stackedbar.png")
DIAG_LOG       <- file.path(OUT_LOG_DIR, "08B_rerun_diagnostics.txt")

# Verdict thresholds (3 tiers)
PASS_CORR     <- 0.30
PASS_PASSRATE <- 0.30
MARGINAL_CORR <- 0.15

# Cluster + histology palettes (same as stages 7 and 8)
CLUSTER_COL <- c(C1 = "#E41A1C", C2 = "#377EB8")

# ── 2. PACKAGES + SEED ──────────────────────────────────────────────────────
setwd(PROJECT_ROOT)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(tibble)
  library(readr)
})

# Note: this seed differs from the main pipeline seed 1262118388
# (Liu 2024). 20260505 marks this as a re-run experiment distinct from the
# canonical run; documented in CLAUDE.md §4.3 deviation rationale (rerun-only).
set.seed(20260505)

# Ensure output dirs exist
for (d in c(OUT_TBL_DIR, OUT_FIG_DIR, OUT_LOG_DIR, CACHE_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Diagnostic log file connection (we'll write incrementally)
diag_lines <- c(
  sprintf("=== Stage 8 / Step B CIBERSORT RERUN diagnostics ===  %s",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "Inputs:",
  sprintf("  EXPR_FILE    : %s", EXPR_FILE),
  sprintf("  CLUSTER_FILE : %s", CLUSTER_FILE),
  sprintf("  LM22_FILE    : %s", LM22_FILE),
  sprintf("  CIBERSORT_R  : %s", CIBERSORT_R),
  "Critical changes vs original Step B:",
  "  1. Mixture: log2(RSEM+1) -> linear (2^x - 1, floor 0)",
  "  2. CIBERSORT: QN = TRUE -> QN = FALSE",
  "  perm = 1000 (unchanged)",
  ""
)

# ── 3. LOAD INPUTS ──────────────────────────────────────────────────────────
cat("[1] Loading expression matrix and cluster assignment...\n")

# Auto-detect TSV vs CSV by file extension
expr_raw <- if (grepl("\\.tsv$|\\.txt$", EXPR_FILE, ignore.case = TRUE)) {
  read_tsv(EXPR_FILE, show_col_types = FALSE)
} else {
  read_csv(EXPR_FILE, show_col_types = FALSE)
}

cluster_raw <- if (grepl("\\.tsv$|\\.txt$", CLUSTER_FILE, ignore.case = TRUE)) {
  read_tsv(CLUSTER_FILE, show_col_types = FALSE)
} else {
  read_csv(CLUSTER_FILE, show_col_types = FALSE)
}

cat(sprintf("    expression: %d rows x %d cols\n", nrow(expr_raw), ncol(expr_raw)))
cat(sprintf("    cluster   : %d rows  (cols: %s)\n",
            nrow(cluster_raw), paste(colnames(cluster_raw), collapse = ", ")))

# Validate required cluster columns
req_cols <- c(CLUSTER_COL_SAMPLE, CLUSTER_COL_CLUSTER, CLUSTER_COL_HIST)
if (!all(req_cols %in% colnames(cluster_raw))) {
  stop(sprintf("Cluster file missing required columns. Expected: %s. Found: %s",
               paste(req_cols, collapse=","), paste(colnames(cluster_raw), collapse=",")))
}
cluster_dt <- cluster_raw %>%
  select(sample    = all_of(CLUSTER_COL_SAMPLE),
         cluster   = all_of(CLUSTER_COL_CLUSTER),
         histology = all_of(CLUSTER_COL_HIST))

# ── 4. GENE SYMBOL CLEANUP ──────────────────────────────────────────────────
cat("[2] Gene symbol cleanup: uppercase, strip version suffix, deduplicate...\n")

gene_col <- colnames(expr_raw)[1]
expr <- expr_raw %>%
  rename(GeneSymbol = !!gene_col) %>%
  mutate(GeneSymbol = toupper(GeneSymbol),
         GeneSymbol = sub("\\..*$", "", GeneSymbol)) %>%
  distinct(GeneSymbol, .keep_all = TRUE) %>%
  filter(!is.na(GeneSymbol) & GeneSymbol != "")

cat(sprintf("    gene rows after cleanup: %d (was %d)\n",
            nrow(expr), nrow(expr_raw)))

# ── 5. === CRITICAL CHANGE 1: LINEAR-SPACE CONVERSION === ───────────────────
cat("[3] === CRITICAL CHANGE 1: Linear-space conversion (2^x - 1, floor 0) ===\n")

sample_cols <- setdiff(colnames(expr), "GeneSymbol")
mat_log <- as.matrix(expr[, sample_cols])
rownames(mat_log) <- expr$GeneSymbol

# log-space stats
ls_min <- min(mat_log, na.rm = TRUE)
ls_med <- median(mat_log, na.rm = TRUE)
ls_max <- max(mat_log, na.rm = TRUE)
cat(sprintf("    log2 space  : min=%.3f  median=%.3f  max=%.3f\n", ls_min, ls_med, ls_max))

# Convert
mat_lin <- 2^mat_log - 1
mat_lin[mat_lin < 0] <- 0  # floor negatives (rare, from numerical edge)

lin_min <- min(mat_lin, na.rm = TRUE)
lin_med <- median(mat_lin, na.rm = TRUE)
lin_max <- max(mat_lin, na.rm = TRUE)
cat(sprintf("    linear space: min=%.3f  median=%.3f  max=%.3f\n", lin_min, lin_med, lin_max))

diag_lines <- c(diag_lines,
  "Mixture-matrix space stats:",
  sprintf("  log2(RSEM+1) input  : min=%.3f  median=%.3f  max=%.3f", ls_min, ls_med, ls_max),
  sprintf("  linear (2^x - 1)    : min=%.3f  median=%.3f  max=%.3f", lin_min, lin_med, lin_max),
  ""
)

# ── 6. SAMPLE INTERSECT WITH CLUSTER ASSIGNMENT ─────────────────────────────
cat("[4] Intersecting expression samples with cluster assignment...\n")

samples_keep <- intersect(colnames(mat_lin), cluster_dt$sample)
cat(sprintf("    samples in expr matrix: %d\n", ncol(mat_lin)))
cat(sprintf("    samples in cluster    : %d\n", nrow(cluster_dt)))
cat(sprintf("    intersection (analysis): %d\n", length(samples_keep)))
stopifnot(length(samples_keep) > 0)

mat_lin_keep <- mat_lin[, samples_keep, drop = FALSE]

cluster_aligned <- cluster_dt %>%
  filter(sample %in% samples_keep) %>%
  arrange(match(sample, samples_keep))
stopifnot(all(cluster_aligned$sample == samples_keep))

# ── 7. WRITE LINEAR-SPACE MIXTURE TO TSV (CIBERSORT input format) ───────────
cat("[5] Writing linear-space mixture file (CIBERSORT input)...\n")
mixture_file <- file.path(CACHE_DIR, "tcga_mixture_linear_for_cibersort.tsv")

mix_tbl <- as_tibble(mat_lin_keep, rownames = "GeneSymbol")
write_tsv(mix_tbl, mixture_file)
cat(sprintf("    wrote: %s (%d genes x %d samples)\n",
            mixture_file, nrow(mix_tbl), ncol(mix_tbl) - 1))

# ── 8. SOURCE CIBERSORT.R + RUN ─────────────────────────────────────────────
cat("[6] Sourcing CIBERSORT.R and running deconvolution ...\n")
cat(sprintf("    source: %s\n", CIBERSORT_R))
source(CIBERSORT_R)  # defines CIBERSORT() function in global env

# === CRITICAL CHANGE 2: QN = FALSE ===
cat("[7] === CRITICAL CHANGE 2: CIBERSORT(perm=1000, QN=FALSE) ===\n")
cat("    (single-threaded on Windows; expect ~5-15 min)\n")
t0 <- Sys.time()
ciber_res <- CIBERSORT(LM22_FILE, mixture_file, perm = 1000, QN = FALSE)
elapsed <- as.numeric(Sys.time() - t0, units = "mins")
cat(sprintf("    CIBERSORT done in %.1f min\n", elapsed))

diag_lines <- c(diag_lines,
  sprintf("CIBERSORT runtime  : %.1f min", elapsed),
  sprintf("CIBERSORT options  : perm = 1000, QN = FALSE"),
  ""
)

# ── 9. DIAGNOSTICS (Correlation, ciber_p, RMSE) ─────────────────────────────
cat("[8] Computing diagnostics...\n")

# ciber_res: rows = sample, cols = 22 cell types + P-value + Correlation + RMSE
ciber_df <- as.data.frame(ciber_res, stringsAsFactors = FALSE)
ciber_df$sample <- rownames(ciber_res)
ciber_df <- ciber_df %>% select(sample, everything())

# Locate metric columns flexibly (some CIBERSORT.R variants name slightly differently)
nm_lower <- tolower(colnames(ciber_df))
col_p    <- colnames(ciber_df)[which(nm_lower %in% c("p-value","p.value","pvalue","p_value"))][1]
col_corr <- colnames(ciber_df)[which(nm_lower %in% c("correlation","pearson correlation"))][1]
col_rmse <- colnames(ciber_df)[which(nm_lower %in% c("rmse","root mean squared error"))][1]
stopifnot(!is.na(col_p), !is.na(col_corr), !is.na(col_rmse))

p_vals    <- as.numeric(ciber_df[[col_p]])
corr_vals <- as.numeric(ciber_df[[col_corr]])
rmse_vals <- as.numeric(ciber_df[[col_rmse]])

diag_metric <- function(label, v) {
  c(sprintf("  %-12s : mean = %8.4f  median = %8.4f  min = %8.4f  max = %8.4f",
            label, mean(v, na.rm=TRUE), median(v, na.rm=TRUE),
            min(v, na.rm=TRUE), max(v, na.rm=TRUE)))
}

cat("    Diagnostic summary (printed to stdout AND to log file):\n")
cat(diag_metric("Correlation", corr_vals), sep = "\n"); cat("\n")
cat(diag_metric("ciber_p",     p_vals),    sep = "\n"); cat("\n")
cat(diag_metric("RMSE",        rmse_vals), sep = "\n"); cat("\n")

mean_corr   <- mean(corr_vals, na.rm = TRUE)
pass_rate_p <- mean(p_vals < 0.05, na.rm = TRUE)
n_pass_p05  <- sum(p_vals < 0.05, na.rm = TRUE)
n_pass_p10  <- sum(p_vals < 0.10, na.rm = TRUE)

diag_lines <- c(diag_lines,
  "Diagnostic stats:",
  diag_metric("Correlation", corr_vals),
  diag_metric("ciber_p",     p_vals),
  diag_metric("RMSE",        rmse_vals),
  sprintf("  ciber_p < 0.05 pass : %d / %d  (%.1f%%)",
          n_pass_p05, length(p_vals), 100 * pass_rate_p),
  sprintf("  ciber_p < 0.10 pass : %d / %d  (%.1f%%)",
          n_pass_p10, length(p_vals), 100 * mean(p_vals < 0.10, na.rm = TRUE)),
  ""
)

# ── 10. 3-TIER VERDICT ──────────────────────────────────────────────────────
cat("[9] 3-tier verdict assessment...\n")

if (mean_corr > PASS_CORR && pass_rate_p > PASS_PASSRATE) {
  verdict   <- "PASS"
  filter_p  <- 0.05
  msg_v     <- sprintf("PASS: mean_corr=%.3f > %.2f AND ciber_p<0.05 pass rate=%.1f%% > %.0f%%",
                       mean_corr, PASS_CORR, 100*pass_rate_p, 100*PASS_PASSRATE)
} else if (mean_corr > MARGINAL_CORR) {
  verdict   <- "MARGINAL"
  filter_p  <- 0.10
  msg_v     <- sprintf("MARGINAL: mean_corr=%.3f in (%.2f, %.2f], using ciber_p<0.10 with caveat",
                       mean_corr, MARGINAL_CORR, PASS_CORR)
} else {
  verdict   <- "FAIL"
  filter_p  <- NA_real_
  msg_v     <- sprintf("FAIL: mean_corr=%.3f <= %.2f. Skipping Wilcoxon and figures, saving raw only.",
                       mean_corr, MARGINAL_CORR)
}
cat(sprintf("    VERDICT: %s\n", verdict))
cat(sprintf("    %s\n", msg_v))

diag_lines <- c(diag_lines,
  "============================================================",
  sprintf("VERDICT: %s", verdict),
  sprintf("  %s", msg_v),
  if (!is.na(filter_p)) sprintf("  Downstream Wilcoxon filter: ciber_p < %.2f", filter_p) else
    "  No downstream Wilcoxon (verdict FAIL).",
  "============================================================",
  ""
)

# ── 11. SAVE per-sample results (always, regardless of verdict) ─────────────
cat("[10] Writing Tab8B_CIBERSORT_results.csv (per-sample, all 181)...\n")

results_out <- ciber_df %>%
  left_join(cluster_aligned, by = "sample") %>%
  rename(sampleID = sample) %>%
  select(sampleID, cluster, histology, everything())

write_csv(results_out, TAB_RESULTS)
cat(sprintf("    [OK] %s\n", TAB_RESULTS))

diag_lines <- c(diag_lines,
  sprintf("Saved: %s (%d samples)", TAB_RESULTS, nrow(results_out)),
  ""
)

# ── 12. If FAIL: write diagnostics, leave existing celltype_stats/figures stale, quit ─
if (verdict == "FAIL") {
  diag_lines <- c(diag_lines,
    "FAIL handling:",
    "  - Tab8B_CIBERSORT_celltype_stats.csv from previous run NOT updated.",
    "  - Fig4B_CIBERSORT_clusters.{pdf,png} from previous run NOT updated.",
    "  - Fig4B_supp_CIBERSORT_stackedbar.{pdf,png} from previous run NOT updated.",
    "  Investigate input space / signature reference compatibility before reuse.",
    ""
  )
  writeLines(diag_lines, DIAG_LOG)
  cat(sprintf("    [OK] %s (diagnostics)\n", DIAG_LOG))
  cat("\n[FAIL] Verdict FAIL — quitting cleanly.\n")
  quit(save = "no", status = 0)
}

# ── 13. PASS / MARGINAL — Wilcoxon C1 vs C2 per cell type ──────────────────
cat(sprintf("[11] Filtering samples (ciber_p < %.2f) and Wilcoxon C1 vs C2 per cell type...\n",
            filter_p))

ciber_pass <- results_out %>% filter(.data[[col_p]] < filter_p)
cat(sprintf("    samples passing filter: %d / %d\n", nrow(ciber_pass), nrow(results_out)))
diag_lines <- c(diag_lines,
  sprintf("Samples after ciber_p<%.2f filter: %d / %d",
          filter_p, nrow(ciber_pass), nrow(results_out)),
  ""
)

# Identify 22 cell type columns (everything except metadata + diagnostic cols)
metadata_cols <- c("sampleID", "cluster", "histology")
diag_cols     <- c(col_p, col_corr, col_rmse)
cell_cols     <- setdiff(colnames(ciber_pass), c(metadata_cols, diag_cols))
cat(sprintf("    cell-type columns: %d\n", length(cell_cols)))

# Wilcoxon per cell type (anonymous lapply, plain script-friendly)
stat_rows <- lapply(cell_cols, function(ct) {
  v1 <- as.numeric(ciber_pass[[ct]][ciber_pass$cluster == "C1"])
  v2 <- as.numeric(ciber_pass[[ct]][ciber_pass$cluster == "C2"])
  p_val <- if (length(v1) >= 3 && length(v2) >= 3) {
    tryCatch(wilcox.test(v1, v2, exact = FALSE)$p.value,
             error = function(e) NA_real_)
  } else NA_real_
  tibble(
    cell_type   = ct,
    n_C1        = length(v1),
    n_C2        = length(v2),
    median_C1   = median(v1, na.rm = TRUE),
    median_C2   = median(v2, na.rm = TRUE),
    mean_C1     = mean(v1,   na.rm = TRUE),
    mean_C2     = mean(v2,   na.rm = TRUE),
    p_value     = p_val
  )
})
cellstat <- bind_rows(stat_rows)
cellstat$p_adj_BH         <- p.adjust(cellstat$p_value, method = "BH")
cellstat$sig_label        <- ifelse(is.na(cellstat$p_adj_BH), "ns",
                              ifelse(cellstat$p_adj_BH < 0.001, "***",
                              ifelse(cellstat$p_adj_BH < 0.01,  "**",
                              ifelse(cellstat$p_adj_BH < 0.05,  "*",  "ns"))))
cellstat$significant_BH_05 <- !is.na(cellstat$p_adj_BH) & cellstat$p_adj_BH < 0.05

write_csv(cellstat, TAB_CELLSTAT)
n_sig <- sum(cellstat$significant_BH_05)
cat(sprintf("    [OK] %s\n", TAB_CELLSTAT))
cat(sprintf("    Significant cell types (BH<0.05): %d / %d\n", n_sig, nrow(cellstat)))

diag_lines <- c(diag_lines,
  sprintf("Wilcoxon results: %d significant cell types (BH<0.05) of %d",
          n_sig, nrow(cellstat)),
  if (verdict == "MARGINAL") "  CAVEAT: verdict MARGINAL — interpret cell-type stats with caution." else "",
  ""
)

# ── 14. FIG 4B: per-cell-type boxplot, 6-col facet, free_y, percent labels ──
cat("[12] Generating Fig4B (per-cell-type boxplot)...\n")

long_b <- ciber_pass %>%
  select(sampleID, cluster, all_of(cell_cols)) %>%
  pivot_longer(cols = all_of(cell_cols),
               names_to = "cell_type",
               values_to = "proportion") %>%
  mutate(cell_type  = factor(cell_type, levels = cell_cols),
         cluster    = factor(cluster, levels = c("C1","C2")),
         proportion = as.numeric(proportion))

# y-position for sig labels: per-facet max + small offset
ymax_b <- long_b %>%
  group_by(cell_type) %>%
  summarise(y_max = max(proportion, na.rm = TRUE), .groups = "drop")

cellstat_plot <- cellstat %>%
  mutate(cell_type = factor(cell_type, levels = cell_cols)) %>%
  left_join(ymax_b, by = "cell_type") %>%
  mutate(y_pos = y_max + 0.05 * abs(y_max) + 0.005)

caveat_subtitle <- if (verdict == "MARGINAL") {
  sprintf("VERDICT MARGINAL  (mean Correlation=%.3f, ciber_p<%.2f filter)", mean_corr, filter_p)
} else {
  sprintf("VERDICT PASS  (mean Correlation=%.3f, ciber_p<%.2f filter)", mean_corr, filter_p)
}

p_box <- ggplot(long_b, aes(x = cluster, y = proportion, fill = cluster)) +
  geom_boxplot(outlier.size = 0.4, color = "black", linewidth = 0.3) +
  geom_text(data = cellstat_plot,
            aes(x = 1.5, y = y_pos, label = sig_label),
            inherit.aes = FALSE, size = 3, fontface = "bold") +
  facet_wrap(~ cell_type, ncol = 6, scales = "free_y") +
  scale_fill_manual(values = CLUSTER_COL, name = "Cluster") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.18)),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Cell-type proportion",
       title    = sprintf("CIBERSORT cell-type proportions by FMG cluster (n = %d, RERUN linear+QN=FALSE)",
                          nrow(ciber_pass)),
       subtitle = caveat_subtitle) +
  theme_bw(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5,
                                     color = ifelse(verdict == "MARGINAL", "#D95F02", "black")),
        strip.text    = element_text(face = "bold", size = 8),
        axis.text.x   = element_text(face = "bold", size = 8),
        axis.text.y   = element_text(size = 7),
        legend.position = "top")

ggsave(FIG_BOX_PDF, plot = p_box, width = 14, height = 9, device = cairo_pdf)
ggsave(FIG_BOX_PNG, plot = p_box, width = 14, height = 9, dpi = 600, device = "png")
cat(sprintf("    [OK] %s\n", FIG_BOX_PDF))
cat(sprintf("    [OK] %s\n", FIG_BOX_PNG))

# ── 15. FIG 4B_supp: per-sample stacked bar, cluster-faceted ────────────────
cat("[13] Generating Fig4B_supp (per-sample stacked bar, cluster-faceted)...\n")

stack_df <- long_b %>%
  arrange(cluster, sampleID) %>%
  mutate(sampleID = factor(sampleID, levels = unique(sampleID)))

# 22-color palette: Set3 (12) + Paired (10)
ct_colors <- c(
  RColorBrewer::brewer.pal(12, "Set3"),
  RColorBrewer::brewer.pal(10, "Paired")
)[seq_len(length(cell_cols))]
names(ct_colors) <- cell_cols

p_stack <- ggplot(stack_df, aes(x = sampleID, y = proportion, fill = cell_type)) +
  geom_col(width = 1) +
  facet_grid(. ~ cluster, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = ct_colors, name = "Cell type") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Proportion",
       title    = sprintf("CIBERSORT per-sample composition (cluster-faceted, n = %d, RERUN linear+QN=FALSE)",
                          nrow(ciber_pass)),
       subtitle = caveat_subtitle) +
  theme_bw(base_size = 9) +
  theme(plot.title       = element_text(face = "bold", size = 11, hjust = 0.5),
        plot.subtitle    = element_text(size = 9, hjust = 0.5,
                                        color = ifelse(verdict == "MARGINAL", "#D95F02", "black")),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        strip.background = element_rect(fill = "#EAEAEA", color = NA),
        strip.text       = element_text(face = "bold"),
        legend.position  = "right",
        legend.text      = element_text(size = 8),
        legend.key.size  = unit(0.4, "cm")) +
  guides(fill = guide_legend(ncol = 1))

ggsave(FIG_STACK_PDF, plot = p_stack, width = 14, height = 6, device = cairo_pdf)
ggsave(FIG_STACK_PNG, plot = p_stack, width = 14, height = 6, dpi = 600, device = "png")
cat(sprintf("    [OK] %s\n", FIG_STACK_PDF))
cat(sprintf("    [OK] %s\n", FIG_STACK_PNG))

# ── 16. Finalize diagnostics log ────────────────────────────────────────────
diag_lines <- c(diag_lines,
  "Outputs written:",
  sprintf("  %s", TAB_RESULTS),
  sprintf("  %s", TAB_CELLSTAT),
  sprintf("  %s", FIG_BOX_PDF),
  sprintf("  %s", FIG_BOX_PNG),
  sprintf("  %s", FIG_STACK_PDF),
  sprintf("  %s", FIG_STACK_PNG),
  sprintf("  %s", DIAG_LOG),
  "",
  sprintf("=== END  %s ===", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
writeLines(diag_lines, DIAG_LOG)
cat(sprintf("[14] [OK] %s\n", DIAG_LOG))

# ── 17. Cleanup: CIBERSORT.R default raw output → cache/ ────────────────────
# Newman lab CIBERSORT.R writes "CIBERSORT-Results.txt" to the working dir by
# default. We keep it as a reproducibility artifact but in cache/, not project
# root.
if (file.exists("CIBERSORT-Results.txt")) {
  cache_path <- file.path(CACHE_DIR, "CIBERSORT-Results.txt")
  file.rename("CIBERSORT-Results.txt", cache_path)
  cat(sprintf("    [moved] CIBERSORT-Results.txt -> %s\n", cache_path))
}

cat(sprintf("\n=== Stage 8 / Step B RERUN done.  VERDICT: %s ===\n", verdict))

###############################################################################
##  Stage 8 / Step B — IMPORT CIBERSORTx (B-mode) result + downstream
##
##  Reason:
##    Vanilla CIBERSORT (both log+QN=TRUE and linear+QN=FALSE) produced
##    deconvolution failure (mean Correlation < 0.03, ~1-3 / 181 samples
##    passing ciber_p<0.05). Cause: TCGA RNA-seq RSEM mixture × LM22
##    Affymetrix microarray reference cross-platform mismatch.
##
##  Solution (B-mode batch correction; Newman et al. 2019 *Nat Biotechnol*):
##    Run CIBERSORTx web tool with Batch correction = enabled, Mode = B
##    (RNA-seq mixture × microarray signature cross-platform mode).
##    Diagnostic: mean Correlation ~0.479, ciber_p<0.05 pass ~99.4% — PASS.
##
##  This script:
##    - Imports CIBERSORTx web result (no local CIBERSORT call)
##    - Runs same 3-tier verdict + Wilcoxon downstream + figures
##    - Overwrites Stage 8 / Step B outputs with corrected results
##
##  CIBERSORTx run options (recorded in diagnostic log header):
##    Job ID            : 3 (Stanford CIBERSORTx)
##    Date              : 2026-05-05 05:05:37
##    Signature matrix  : LM22.update-gene-symbols.txt (CIBERSORTx updated)
##    Mixture           : tcga_mixture_linear_for_cibersort.tsv (linear-space)
##    Batch correction  : enabled, B-mode
##    Disable QN        : true
##    Run mode          : relative
##    Permutations      : 100
##
##  Outputs (overwrite):
##    Tab8B_CIBERSORT_results.csv           per-sample 22 cell types + p/corr/RMSE + cluster + histology
##    Tab8B_CIBERSORT_celltype_stats.csv    per-cell-type C1 vs C2 Wilcoxon
##    Fig4B_CIBERSORT_clusters.{pdf,png}    boxplot, 6-col facet  (title: "CIBERSORTx (B-mode)")
##    Fig4B_supp_CIBERSORT_stackedbar.{pdf,png}  stacked bar      (title: "CIBERSORTx (B-mode)")
##    08B_cibersortx_import_diagnostics.txt
##
##  Untouched: Tab8A, Tab8C, Tab8D, Tab8E, Fig4A, Fig4C, Fig4D, Fig4E
##

# ── 1. PATHS / VARIABLES (user-editable) ─────────────────────────────────
PROJECT_ROOT <- getwd()   # run from the repository root (see README.md)

# Inputs
CIBERX_FILE    <- "data/cibersortx_LM22_TCGA_ESCA_Bmode.csv"
CLUSTER_FILE   <- "tables/Tab7_cluster_assignment.csv"

CLUSTER_COL_SAMPLE  <- "sampleID"
CLUSTER_COL_CLUSTER <- "cluster"
CLUSTER_COL_HIST    <- "histology"

# Outputs
OUT_TBL_DIR    <- "tables"
OUT_FIG_DIR    <- "figures"
OUT_LOG_DIR    <- "logs"

TAB_RESULTS    <- file.path(OUT_TBL_DIR, "Tab8B_CIBERSORT_results.csv")
TAB_CELLSTAT   <- file.path(OUT_TBL_DIR, "Tab8B_CIBERSORT_celltype_stats.csv")
FIG_BOX_PDF    <- file.path(OUT_FIG_DIR, "Fig4B_CIBERSORT_clusters.pdf")
FIG_BOX_PNG    <- file.path(OUT_FIG_DIR, "Fig4B_CIBERSORT_clusters.png")
FIG_STACK_PDF  <- file.path(OUT_FIG_DIR, "Fig4B_supp_CIBERSORT_stackedbar.pdf")
FIG_STACK_PNG  <- file.path(OUT_FIG_DIR, "Fig4B_supp_CIBERSORT_stackedbar.png")
DIAG_LOG       <- file.path(OUT_LOG_DIR, "08B_cibersortx_import_diagnostics.txt")

# Verdict thresholds (same as the rerun script)
PASS_CORR     <- 0.30
PASS_PASSRATE <- 0.30
MARGINAL_CORR <- 0.15

# Cluster palette
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
  library(RColorBrewer)
})

set.seed(20260505)

for (d in c(OUT_TBL_DIR, OUT_FIG_DIR, OUT_LOG_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Diagnostics log header
diag_lines <- c(
  sprintf("=== Stage 8 / Step B — CIBERSORTx (B-mode) IMPORT diagnostics ===  %s",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "CIBERSORTx run options (Stanford web tool):",
  "  Job ID            : 3",
  "  Date              : 2026-05-05 05:05:37",
  "  Signature matrix  : LM22.update-gene-symbols.txt (CIBERSORTx updated)",
  "  Mixture           : tcga_mixture_linear_for_cibersort.tsv (linear-space)",
  "  Batch correction  : enabled, B-mode",
  "  Disable QN        : true",
  "  Run mode          : relative",
  "  Permutations      : 100",
  "",
  "Reference: Newman AM et al. 2019 Nat Biotechnol -- CIBERSORTx batch correction",
  "  (B-mode is for cross-platform: RNA-seq mixture x microarray signature)",
  "",
  "Inputs to this import step:",
  sprintf("  CIBERX_FILE   : %s", CIBERX_FILE),
  sprintf("  CLUSTER_FILE  : %s", CLUSTER_FILE),
  ""
)

# ── 3. LOAD CIBERSORTx RESULT + CLUSTER ─────────────────────────────────────
cat("[1] Loading CIBERSORTx result CSV and cluster assignment...\n")

cx_raw <- read_csv(CIBERX_FILE, show_col_types = FALSE)
cluster_raw <- read_csv(CLUSTER_FILE, show_col_types = FALSE)

cat(sprintf("    CIBERSORTx: %d rows x %d cols\n", nrow(cx_raw), ncol(cx_raw)))
cat(sprintf("    cluster   : %d rows  (cols: %s)\n",
            nrow(cluster_raw), paste(colnames(cluster_raw), collapse = ", ")))

# Validate cluster columns
req_cols <- c(CLUSTER_COL_SAMPLE, CLUSTER_COL_CLUSTER, CLUSTER_COL_HIST)
if (!all(req_cols %in% colnames(cluster_raw))) {
  stop(sprintf("Cluster file missing required columns. Expected: %s",
               paste(req_cols, collapse = ", ")))
}

cluster_dt <- cluster_raw %>%
  select(sampleID  = all_of(CLUSTER_COL_SAMPLE),
         cluster   = all_of(CLUSTER_COL_CLUSTER),
         histology = all_of(CLUSTER_COL_HIST))

# ── 4. RENAME Mixture -> sampleID + DETECT METRIC COLUMNS (case-insensitive) ─
cat("[2] Locating Mixture / P-value / Correlation / RMSE columns...\n")

# Mixture column: usually "Mixture" but be flexible
mix_col <- intersect(colnames(cx_raw), c("Mixture", "mixture", "Sample", "sample", "sampleID"))
if (length(mix_col) == 0) {
  # Default to first column
  mix_col <- colnames(cx_raw)[1]
  cat(sprintf("    [warn] No standard 'Mixture' column found, defaulting to first: %s\n",
              mix_col))
}
mix_col <- mix_col[1]
cat(sprintf("    sample column   : %s\n", mix_col))

cx <- cx_raw %>% rename(sampleID = !!mix_col)

nm_lower <- tolower(colnames(cx))
col_p    <- colnames(cx)[which(nm_lower %in% c("p-value","p.value","pvalue","p_value"))][1]
col_corr <- colnames(cx)[which(nm_lower %in% c("correlation","pearson correlation"))][1]
col_rmse <- colnames(cx)[which(nm_lower %in% c("rmse","root mean squared error"))][1]
if (is.na(col_p) || is.na(col_corr) || is.na(col_rmse)) {
  stop(sprintf("Required metric columns not found. Found: P=%s Corr=%s RMSE=%s",
               col_p, col_corr, col_rmse))
}
cat(sprintf("    P-value col     : %s\n", col_p))
cat(sprintf("    Correlation col : %s\n", col_corr))
cat(sprintf("    RMSE col        : %s\n", col_rmse))

# 22 cell type columns = everything except sampleID + 3 metric cols
cell_cols <- setdiff(colnames(cx), c("sampleID", col_p, col_corr, col_rmse))
cat(sprintf("    cell-type cols  : %d\n", length(cell_cols)))

# ── 5. SAMPLE-ID NORMALIZATION (handle . -> - reversion if needed) ──────────
cat("[3] Sample ID matching with cluster assignment...\n")

n_match_orig <- sum(cx$sampleID %in% cluster_dt$sampleID)
cat(sprintf("    direct match (cibersortx vs cluster) : %d / %d\n",
            n_match_orig, nrow(cx)))

if (n_match_orig < nrow(cx) * 0.5) {
  cat("    [auto-fix] Many samples unmatched; trying '.' -> '-' substitution...\n")
  cx_alt <- cx %>% mutate(sampleID = gsub("\\.", "-", sampleID))
  n_match_alt <- sum(cx_alt$sampleID %in% cluster_dt$sampleID)
  cat(sprintf("    after '.' -> '-': match = %d / %d\n", n_match_alt, nrow(cx_alt)))
  if (n_match_alt > n_match_orig) {
    cat("    [fixed] using corrected sample IDs.\n")
    cx <- cx_alt
  }
}

samples_keep <- intersect(cx$sampleID, cluster_dt$sampleID)
cat(sprintf("    final intersection: %d samples\n", length(samples_keep)))
stopifnot(length(samples_keep) > 0)

# ── 6. DIAGNOSTICS ──────────────────────────────────────────────────────────
cat("[4] Computing diagnostics on imported CIBERSORTx values...\n")

p_vals    <- as.numeric(cx[[col_p]])
corr_vals <- as.numeric(cx[[col_corr]])
rmse_vals <- as.numeric(cx[[col_rmse]])

diag_metric <- function(label, v) {
  c(sprintf("  %-12s : mean = %8.4f  median = %8.4f  min = %8.4f  max = %8.4f",
            label, mean(v, na.rm = TRUE), median(v, na.rm = TRUE),
            min(v, na.rm = TRUE), max(v, na.rm = TRUE)))
}

cat(diag_metric("Correlation", corr_vals), sep = "\n"); cat("\n")
cat(diag_metric("ciber_p",     p_vals),    sep = "\n"); cat("\n")
cat(diag_metric("RMSE",        rmse_vals), sep = "\n"); cat("\n")

mean_corr   <- mean(corr_vals, na.rm = TRUE)
pass_rate_p <- mean(p_vals < 0.05, na.rm = TRUE)
n_pass_p05  <- sum(p_vals < 0.05, na.rm = TRUE)
n_pass_p10  <- sum(p_vals < 0.10, na.rm = TRUE)

cat(sprintf("    ciber_p < 0.05 pass : %d / %d  (%.1f%%)\n",
            n_pass_p05, length(p_vals), 100 * pass_rate_p))
cat(sprintf("    ciber_p < 0.10 pass : %d / %d  (%.1f%%)\n",
            n_pass_p10, length(p_vals), 100 * mean(p_vals < 0.10, na.rm = TRUE)))

diag_lines <- c(diag_lines,
  "Diagnostic stats (imported CIBERSORTx values, all 181 samples):",
  diag_metric("Correlation", corr_vals),
  diag_metric("ciber_p",     p_vals),
  diag_metric("RMSE",        rmse_vals),
  sprintf("  ciber_p < 0.05 pass : %d / %d  (%.1f%%)",
          n_pass_p05, length(p_vals), 100 * pass_rate_p),
  sprintf("  ciber_p < 0.10 pass : %d / %d  (%.1f%%)",
          n_pass_p10, length(p_vals), 100 * mean(p_vals < 0.10, na.rm = TRUE)),
  ""
)

# ── 7. 3-TIER VERDICT ──────────────────────────────────────────────────────
cat("[5] 3-tier verdict assessment...\n")

if (mean_corr > PASS_CORR && pass_rate_p > PASS_PASSRATE) {
  verdict   <- "PASS"
  filter_p  <- 0.05
  msg_v     <- sprintf("PASS: mean_corr=%.3f > %.2f AND ciber_p<0.05 pass rate=%.1f%% > %.0f%%",
                       mean_corr, PASS_CORR, 100 * pass_rate_p, 100 * PASS_PASSRATE)
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
cat(sprintf("    VERDICT: %s\n    %s\n", verdict, msg_v))

diag_lines <- c(diag_lines,
  "============================================================",
  sprintf("VERDICT: %s", verdict),
  sprintf("  %s", msg_v),
  if (!is.na(filter_p)) sprintf("  Downstream Wilcoxon filter: ciber_p < %.2f", filter_p) else
    "  No downstream Wilcoxon (verdict FAIL).",
  "============================================================",
  ""
)

# ── 8. SAVE per-sample results (always) ─────────────────────────────────────
cat("[6] Writing Tab8B_CIBERSORT_results.csv (per-sample, all 181)...\n")

results_out <- cx %>%
  left_join(cluster_dt, by = "sampleID") %>%
  select(sampleID, cluster, histology, everything())

write_csv(results_out, TAB_RESULTS)
cat(sprintf("    [OK] %s\n", TAB_RESULTS))
diag_lines <- c(diag_lines,
  sprintf("Saved: %s (%d samples)", TAB_RESULTS, nrow(results_out)),
  ""
)

# ── 9. FAIL handling ────────────────────────────────────────────────────────
if (verdict == "FAIL") {
  diag_lines <- c(diag_lines,
    "FAIL handling: celltype_stats and figures from previous run NOT updated.",
    "")
  writeLines(diag_lines, DIAG_LOG)
  cat(sprintf("    [OK] %s (diagnostics)\n", DIAG_LOG))
  cat("[FAIL] Verdict FAIL — quitting cleanly.\n")
  quit(save = "no", status = 0)
}

# ── 10. PASS / MARGINAL: filter + Wilcoxon C1 vs C2 per cell type ──────────
cat(sprintf("[7] Filtering ciber_p < %.2f and Wilcoxon C1 vs C2 per cell type...\n",
            filter_p))

ciber_pass <- results_out %>% filter(.data[[col_p]] < filter_p)
cat(sprintf("    samples passing filter: %d / %d\n", nrow(ciber_pass), nrow(results_out)))
diag_lines <- c(diag_lines,
  sprintf("Samples after ciber_p<%.2f filter: %d / %d",
          filter_p, nrow(ciber_pass), nrow(results_out)),
  ""
)

# Per-cell-type Wilcoxon
stat_rows <- lapply(cell_cols, function(ct) {
  v1 <- as.numeric(ciber_pass[[ct]][ciber_pass$cluster == "C1"])
  v2 <- as.numeric(ciber_pass[[ct]][ciber_pass$cluster == "C2"])
  v1 <- v1[!is.na(v1)]; v2 <- v2[!is.na(v2)]
  p_val <- if (length(v1) >= 3 && length(v2) >= 3 &&
               (length(unique(c(v1, v2))) > 1)) {
    tryCatch(wilcox.test(v1, v2, exact = FALSE)$p.value,
             error = function(e) NA_real_)
  } else NA_real_
  tibble(
    cell_type = ct,
    n_C1      = length(v1),
    n_C2      = length(v2),
    median_C1 = median(v1, na.rm = TRUE),
    median_C2 = median(v2, na.rm = TRUE),
    mean_C1   = mean(  v1, na.rm = TRUE),
    mean_C2   = mean(  v2, na.rm = TRUE),
    p_value   = p_val
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

# Print top 5 by smallest p_adj_BH
top5 <- cellstat %>%
  arrange(p_adj_BH) %>%
  head(5) %>%
  select(cell_type, n_C1, n_C2, median_C1, median_C2, p_adj_BH, sig_label)
cat("    Top 5 by smallest p_adj_BH:\n")
print(as.data.frame(top5))

diag_lines <- c(diag_lines,
  sprintf("Wilcoxon results: %d significant cell types (BH<0.05) of %d",
          n_sig, nrow(cellstat)),
  if (verdict == "MARGINAL") "  CAVEAT: verdict MARGINAL — interpret with caution." else "",
  "Top 5 by smallest p_adj_BH:",
  capture.output(print(as.data.frame(top5))),
  ""
)

# ── 11. FIG 4B: per-cell-type boxplot, 6-col facet ──────────────────────────
cat("[8] Generating Fig4B (per-cell-type boxplot)...\n")

long_b <- ciber_pass %>%
  select(sampleID, cluster, all_of(cell_cols)) %>%
  pivot_longer(cols = all_of(cell_cols),
               names_to  = "cell_type",
               values_to = "proportion") %>%
  mutate(cell_type  = factor(cell_type, levels = cell_cols),
         cluster    = factor(cluster, levels = c("C1","C2")),
         proportion = as.numeric(proportion))

ymax_b <- long_b %>%
  group_by(cell_type) %>%
  summarise(y_max = max(proportion, na.rm = TRUE), .groups = "drop")

cellstat_plot <- cellstat %>%
  mutate(cell_type = factor(cell_type, levels = cell_cols)) %>%
  left_join(ymax_b, by = "cell_type") %>%
  mutate(y_pos = y_max + 0.05 * abs(y_max) + 0.005)

caveat_subtitle <- if (verdict == "MARGINAL") {
  sprintf("VERDICT MARGINAL  (mean Correlation=%.3f, ciber_p<%.2f filter)",
          mean_corr, filter_p)
} else {
  sprintf("VERDICT PASS  (mean Correlation=%.3f, ciber_p<%.2f filter, B-mode batch correction)",
          mean_corr, filter_p)
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
       title    = sprintf("CIBERSORTx (B-mode) cell-type proportions by FMG cluster (n = %d)",
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
cat(sprintf("    [OK] %s\n    [OK] %s\n", FIG_BOX_PDF, FIG_BOX_PNG))

# ── 12. FIG 4B_supp: per-sample stacked bar, cluster-faceted ───────────────
cat("[9] Generating Fig4B_supp (per-sample stacked bar, cluster-faceted)...\n")

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
       title    = sprintf("CIBERSORTx (B-mode) per-sample composition (cluster-faceted, n = %d)",
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
cat(sprintf("    [OK] %s\n    [OK] %s\n", FIG_STACK_PDF, FIG_STACK_PNG))

# ── 13. FINALIZE DIAGNOSTIC LOG ─────────────────────────────────────────────
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
cat(sprintf("[10] [OK] %s\n", DIAG_LOG))

cat(sprintf("\n=== Stage 8 / Step B IMPORT done.  VERDICT: %s ===\n", verdict))

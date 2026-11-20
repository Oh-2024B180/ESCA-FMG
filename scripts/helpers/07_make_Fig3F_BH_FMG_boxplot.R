###############################################################################
##  Fig 3F (alt) — 18 FMG mRNA expression: C1 vs C2 boxplot with BH-adj sig
##
##  Inputs:
##    data/expr_ESCA_tumor.csv                          (log2(RSEM+1), col1 = gene SYMBOL)
##    tables/Tab7_cluster_assignment.csv        (sampleID, cluster, histology, OS, OS.time)
##    tables/Tab1_2_FMG_list_primary.csv         (18 FMG, col `gene`)
##
##  Outputs:
##    tables/Tab7C_FMG_expr_diff_clusters_C1vsC2.csv
##    figures/Fig3F_supp_FMG_expr_diff_BH.{pdf,png}
##    logs/3F2_FMG_boxplot_diagnostics.txt

# ── 1. PATHS / VARIABLES ────────────────────────────────────────────────────
PROJECT_ROOT <- getwd()   # run from the repository root (see README.md)

INPUT_EXPR    <- "data/expr_ESCA_tumor.csv"
INPUT_CLUSTER <- "tables/Tab7_cluster_assignment.csv"
INPUT_FMG     <- "tables/Tab1_2_FMG_list_primary.csv"
# Use Stage 4 Tab4 mutation prevalence order for x-axis
INPUT_TAB4    <- "tables/Tab4_TCGA_mutation_frequency.csv"

OUT_TBL  <- "tables/Tab7C_FMG_expr_diff_clusters_C1vsC2.csv"
OUT_PDF  <- "figures/Fig3F_supp_FMG_expr_diff_BH.pdf"
OUT_PNG  <- "figures/Fig3F_supp_FMG_expr_diff_BH.png"
OUT_LOG  <- "logs/3F2_FMG_boxplot_diagnostics.txt"

# Stage 7 CLUSTER_COLORS (Dark2 palette)
CLUSTER_COLORS <- c("C1" = "#1B9E77", "C2" = "#D95F02")

ADJ_P_CUTOFF <- 0.05

# ── 2. PACKAGES + SEED ──────────────────────────────────────────────────────
setwd(PROJECT_ROOT)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ragg)
})

set.seed(1262118388)

# ── 3. LOAD INPUTS ──────────────────────────────────────────────────────────
cat("[1] Loading expression, cluster, FMG list, and Stage 4 gene order...\n")

expr_dt    <- fread(INPUT_EXPR)
cluster_dt <- fread(INPUT_CLUSTER)
fmg_dt     <- fread(INPUT_FMG)
tab4_dt    <- fread(INPUT_TAB4)

fmg_genes <- fmg_dt$gene
gene_order_x <- tab4_dt$gene  # Stage 4 mutation prevalence order

cat(sprintf("    expr      : %d genes x %d cols\n", nrow(expr_dt), ncol(expr_dt)))
cat(sprintf("    cluster   : %d samples (C1=%d, C2=%d)\n",
            nrow(cluster_dt), sum(cluster_dt$cluster=="C1"), sum(cluster_dt$cluster=="C2")))
cat(sprintf("    FMG       : %d\n", length(fmg_genes)))
cat(sprintf("    gene order (x-axis): %s\n", paste(gene_order_x, collapse=", ")))

# ── 4. SUBSET expression to 18 FMG x 181 cluster samples ────────────────────
cat("[2] Subsetting expression to 18 FMG and 181 Tab7 samples...\n")

key_g <- colnames(expr_dt)[1]
setnames(expr_dt, key_g, "gene_id")

samples_keep <- intersect(cluster_dt$sampleID, colnames(expr_dt))
stopifnot(length(samples_keep) == nrow(cluster_dt))

fmg_expr <- expr_dt[gene_id %in% fmg_genes, c("gene_id", samples_keep), with = FALSE]
stopifnot(all(fmg_genes %in% fmg_expr$gene_id))
cat(sprintf("    FMG expression rows: %d / %d  (samples = %d)\n",
            nrow(fmg_expr), length(fmg_genes), length(samples_keep)))

# ── 5. LONG FORMAT (gene, sample, cluster, expression) ──────────────────────
cat("[3] Reshaping to long format and joining cluster info...\n")

long_dt <- melt(fmg_expr, id.vars = "gene_id", variable.name = "sampleID",
                value.name = "expression", variable.factor = FALSE)
setnames(long_dt, "gene_id", "gene")

long_dt <- merge(long_dt,
                 cluster_dt[, .(sampleID, cluster)],
                 by = "sampleID", all.x = FALSE)

long_dt[, gene    := factor(gene, levels = gene_order_x)]
long_dt[, cluster := factor(cluster, levels = c("C1","C2"))]
long_dt[, expression := as.numeric(expression)]

cat(sprintf("    long_dt rows: %d  (= 18 genes x 181 samples = %d)\n",
            nrow(long_dt), 18 * 181))

# ── 6. WILCOXON C1 vs C2 PER GENE + BH ──────────────────────────────────────
cat("[4] Wilcoxon rank-sum per gene + BH adjustment...\n")

stat_rows <- lapply(gene_order_x, function(g) {
  v1 <- long_dt[gene == g & cluster == "C1", expression]
  v2 <- long_dt[gene == g & cluster == "C2", expression]
  v1 <- v1[!is.na(v1)]; v2 <- v2[!is.na(v2)]
  p_val <- if (length(v1) >= 3 && length(v2) >= 3) {
    tryCatch(wilcox.test(v1, v2, exact = FALSE)$p.value,
             error = function(e) NA_real_)
  } else NA_real_
  data.table(
    gene               = g,
    n_C1               = length(v1),
    n_C2               = length(v2),
    median_C1          = median(v1, na.rm = TRUE),
    median_C2          = median(v2, na.rm = TRUE),
    mean_C1            = mean(  v1, na.rm = TRUE),
    mean_C2            = mean(  v2, na.rm = TRUE),
    log2FC_C1minusC2   = median(v1, na.rm = TRUE) - median(v2, na.rm = TRUE),
    p_value            = p_val
  )
})
stat_dt <- rbindlist(stat_rows)
stat_dt[, p_adj_BH  := p.adjust(p_value, method = "BH")]

# Significance label (4-tier: ****/***/**/*/ns)
stat_dt[, sig_label := fifelse(is.na(p_adj_BH),         "ns",
                       fifelse(p_adj_BH < 0.0001,        "****",
                       fifelse(p_adj_BH < 0.001,         "***",
                       fifelse(p_adj_BH < 0.01,          "**",
                       fifelse(p_adj_BH < 0.05,          "*",   "ns")))))]

# Direction (C1>C2 / C2>C1 / ns)
stat_dt[, direction := fifelse(sig_label == "ns",                "ns",
                       fifelse(log2FC_C1minusC2 > 0,             "C1>C2",
                       fifelse(log2FC_C1minusC2 < 0,             "C2>C1", "ns")))]

# Reorder rows by Stage 4 gene order for output
stat_dt <- stat_dt[match(gene_order_x, gene)]

fwrite(stat_dt, OUT_TBL)
cat(sprintf("    [OK] %s\n", OUT_TBL))

n_sig    <- sum(stat_dt$p_adj_BH < ADJ_P_CUTOFF, na.rm = TRUE)
n_C1up   <- sum(stat_dt$direction == "C1>C2", na.rm = TRUE)
n_C2up   <- sum(stat_dt$direction == "C2>C1", na.rm = TRUE)
n_ns     <- sum(stat_dt$direction == "ns",    na.rm = TRUE)

cat(sprintf("    Significant (BH<%.2f): %d / %d\n", ADJ_P_CUTOFF, n_sig, nrow(stat_dt)))
cat(sprintf("    Direction: C1>C2 = %d,  C2>C1 = %d,  ns = %d\n", n_C1up, n_C2up, n_ns))

# ── 7. BUILD SINGLE-PANEL x-DODGE BOXPLOT ───────────────────────────────────
cat("[5] Building boxplot (single panel, 18 genes x C1/C2 dodge)...\n")

# Per-gene y_max for label positioning
ymax_dt <- long_dt[, .(y_max = max(expression, na.rm = TRUE)), by = gene]
stat_plot <- merge(stat_dt[, .(gene, sig_label)],
                   ymax_dt, by = "gene")
stat_plot[, gene  := factor(gene, levels = gene_order_x)]
stat_plot[, y_pos := y_max + 0.4]   # small offset above per-gene max (log2 scale)

p <- ggplot(long_dt, aes(x = gene, y = expression, fill = cluster)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7, color = "black", linewidth = 0.3,
               outlier.size = 0.5, outlier.alpha = 0.5) +
  geom_text(data = stat_plot,
            aes(x = gene, y = y_pos, label = sig_label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  scale_fill_manual(values = CLUSTER_COLORS, name = "Cluster") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(
    x     = NULL,
    y     = "mRNA expression level"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y      = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position  = "top"
  )

# ── 8. SAVE FIGURE (PDF cairo + PNG ragg @ 600 dpi) ─────────────────────────
cat("[6] Saving figure (PDF cairo + PNG ragg @ 600 dpi)...\n")

W <- 14
H <- 6

ggsave(OUT_PDF, plot = p, width = W, height = H, device = cairo_pdf)
cat(sprintf("    [OK] %s\n", OUT_PDF))

ggsave(OUT_PNG, plot = p, width = W, height = H, dpi = 600, device = ragg::agg_png)
cat(sprintf("    [OK] %s\n", OUT_PNG))

# ── 9. DIAGNOSTICS LOG ──────────────────────────────────────────────────────
cat("[7] Writing diagnostics log...\n")

diag_lines <- c(
  sprintf("=== Fig 3F (BH) — 18 FMG cluster expression diagnostics === %s",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Inputs:",
  sprintf("  EXPR    : %s", INPUT_EXPR),
  sprintf("  CLUSTER : %s", INPUT_CLUSTER),
  sprintf("  FMG     : %s", INPUT_FMG),
  sprintf("  ORDER   : Stage 4 mutation prevalence (%s)", INPUT_TAB4),
  "",
  "Sample distribution:",
  sprintf("  Total: %d  |  C1: %d  |  C2: %d",
          nrow(cluster_dt),
          sum(cluster_dt$cluster == "C1"),
          sum(cluster_dt$cluster == "C2")),
  "",
  "Three changes vs raw-p version:",
  "  1. raw p stars -> BH-adjusted p stars",
  "  2. 'Cluster 1'/'Cluster 2' -> 'C1'/'C2'",
  "  3. cluster colors -> Stage 7 Dark2 (C1 #1B9E77 green, C2 #D95F02 orange)",
  "",
  sprintf("BH<%.2f result summary:", ADJ_P_CUTOFF),
  sprintf("  Significant genes : %d / %d", n_sig, nrow(stat_dt)),
  sprintf("  Direction         : C1>C2 = %d  |  C2>C1 = %d  |  ns = %d",
          n_C1up, n_C2up, n_ns),
  "",
  "Per-gene stats (Stage 4 mutation prevalence order):",
  capture.output(print(stat_dt[, .(gene, n_C1, n_C2, median_C1, median_C2,
                                   log2FC_C1minusC2, p_value, p_adj_BH,
                                   sig_label, direction)])),
  "",
  "Outputs written:",
  sprintf("  %s", OUT_TBL),
  sprintf("  %s", OUT_PDF),
  sprintf("  %s", OUT_PNG),
  sprintf("  %s", OUT_LOG),
  "",
  sprintf("=== END  %s ===", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
writeLines(diag_lines, OUT_LOG)
cat(sprintf("    [OK] %s\n", OUT_LOG))

cat("\n=== Fig 3F (BH) generation done ===\n")

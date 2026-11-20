###############################################################################
## Stage 6: Tumor vs Normal Expression Difference (18 Primary FMG)
##
## Inputs:
## data/expr_ESCA_tumor.csv (TCGA tumor, 181 samples, log2(RSEM+1))
## data/expr_GTEx_esophagus.csv (GTEx esophagus all subtypes, 655 samples)
## data/GTEx_esophagus_subtype_map.csv (sample_expr -> subtype)
## tables/Tab1_2_FMG_list_primary.csv (18 Primary FMG)
##
## Outputs (tables/):
## Tab6_expression_tumor_vs_normal.csv
## columns: gene, n_tumor, n_normal,
## median_tumor, median_normal, log2FC,
## p_value (Wilcoxon), p_adj_BH, sig_label
##
## Outputs (figures/, PDF + PNG @ 600 dpi):
## Fig2C_TCGA_vs_GTEx_mucosa.{pdf,png}
###############################################################################

# 0. Working directory + output
# run from the repository root (see README.md)

INPUT_TUMOR <- "data/expr_ESCA_tumor.csv"
INPUT_GTEX <- "data/expr_GTEx_esophagus.csv"
INPUT_MAP <- "data/GTEx_esophagus_subtype_map.csv"
INPUT_FMG <- "tables/Tab1_2_FMG_list_primary.csv"
INPUT_TAB4 <- "tables/Tab4_TCGA_mutation_frequency.csv" # for panel order

TBLDIR <- "tables"
FIGDIR <- "figures"
LOGDIR <- "logs"

for (d in c(TBLDIR, FIGDIR, LOGDIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

FIG_STEM <- "Fig2C_TCGA_vs_GTEx_mucosa"
TBL_OUT <- file.path(TBLDIR, "Tab6_expression_tumor_vs_normal.csv")

# Sensitivity flag
SENSITIVITY_TCGA_NORMAL <- FALSE # TRUE TCGA matched normal addition (current datanone)

# Stats cutoffs
ADJ_P_CUTOFF <- 0.05 # Stage 5B

# 0.5. Reproducibility
set.seed(1262118388)

# 0.6. Logging start
log_file <- file.path(LOGDIR, "06_diagnostics.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat(sprintf("=== Stage 6 Pipeline start (%s) ===\n", Sys.time()))
cat(sprintf("SENSITIVITY_TCGA_NORMAL : %s\n", SENSITIVITY_TCGA_NORMAL))
cat(sprintf("BH adj.p cutoff : %.3f\n\n", ADJ_P_CUTOFF))

# 1. load
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ragg)
})

# 2. Load input
cat(">>> STEP 2: Load inputs\n")
expr_t <- fread(INPUT_TUMOR)
expr_g <- fread(INPUT_GTEX)
map_g <- fread(INPUT_MAP)
fmg <- fread(INPUT_FMG)
tab4 <- fread(INPUT_TAB4)

fmg_genes <- fmg$gene
gene_order_fig <- tab4$gene # Stage 4 mutation prevalence order

cat(sprintf(" Tumor expr : %d genes x %d cols\n", nrow(expr_t), ncol(expr_t)))
cat(sprintf(" GTEx expr : %d genes x %d cols\n", nrow(expr_g), ncol(expr_g)))
cat(sprintf(" GTEx subtype map : %d rows\n", nrow(map_g)))
cat(sprintf(" 18 FMG : %s\n", paste(fmg_genes, collapse=",")))
cat("\n")

# 3. Sample filter
cat(">>> STEP 3: Sample filtering\n")
# Tumor = all TCGA-* columns
tumor_samples <- grep("^TCGA-", colnames(expr_t), value = TRUE)
N_TUMOR <- length(tumor_samples)
cat(sprintf(" Tumor (TCGA-ESCA -01) : %d\n", N_TUMOR))

# Normal = GTEx Mucosa subtype, intersect with GTEx expr columns
mucosa_samples_map <- map_g[subtype == "Esophagus - Mucosa", sample_expr]
gtex_expr_cols <- colnames(expr_g)[-1]
normal_samples <- intersect(mucosa_samples_map, gtex_expr_cols)
N_NORMAL <- length(normal_samples)
cat(sprintf(" Normal (GTEx Esophagus-Mucosa): %d\n", N_NORMAL))

if (SENSITIVITY_TCGA_NORMAL) {
  stop("SENSITIVITY_TCGA_NORMAL=TRUE: TCGA matched normal expression file not yet available. Update script first.")
}
cat("\n")

# 4. The 18 FMG, long format
cat(">>> STEP 4: Extract 18 FMG expression and reshape to long format\n")

# column('sample')
key_t <- colnames(expr_t)[1]; setnames(expr_t, key_t, "gene_id")
key_g <- colnames(expr_g)[1]; setnames(expr_g, key_g, "gene_id")

t_fmg <- expr_t[gene_id %in% fmg_genes, c("gene_id", tumor_samples), with = FALSE]
g_fmg <- expr_g[gene_id %in% fmg_genes, c("gene_id", normal_samples), with = FALSE]

stopifnot(all(fmg_genes %in% t_fmg$gene_id))
stopifnot(all(fmg_genes %in% g_fmg$gene_id))
cat(sprintf(" FMG rows in tumor expr : %d / %d\n", nrow(t_fmg), length(fmg_genes)))
cat(sprintf(" FMG rows in GTEx expr : %d / %d\n", nrow(g_fmg), length(fmg_genes)))

# Long format
t_long <- melt(t_fmg, id.vars = "gene_id", variable.name = "sample",
               value.name = "value", variable.factor = FALSE)
t_long[, group := "Tumor"]

g_long <- melt(g_fmg, id.vars = "gene_id", variable.name = "sample",
               value.name = "value", variable.factor = FALSE)
g_long[, group := "Normal"]

long_dt <- rbind(t_long, g_long)
setnames(long_dt, "gene_id", "gene")
long_dt[, gene := factor(gene, levels = gene_order_fig)]
long_dt[, group := factor(group, levels = c("Normal", "Tumor"))]

cat(sprintf(" long_dt rows : %d (= 18 genes x (%d tumor + %d normal))\n",
            nrow(long_dt), N_TUMOR, N_NORMAL))
cat("\n")

# 5. Wilcoxon test + BH correction
cat(">>> STEP 5: Wilcoxon rank-sum + BH adjustment\n")

stat_dt <- long_dt[, .(
  n_tumor = sum(group == "Tumor" & !is.na(value)),
  n_normal = sum(group == "Normal" & !is.na(value)),
  median_tumor = median(value[group == "Tumor"], na.rm = TRUE),
  median_normal = median(value[group == "Normal"], na.rm = TRUE),
  p_value = tryCatch(
    wilcox.test(value[group == "Tumor"], value[group == "Normal"],
                exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
), by = gene]

stat_dt[, log2FC := median_tumor - median_normal] # log2
stat_dt[, p_adj_BH := p.adjust(p_value, method = "BH")]
stat_dt[, sig_label := dplyr::case_when(
  is.na(p_adj_BH) ~ "ns",
  p_adj_BH < 0.001 ~ "***",
  p_adj_BH < 0.01 ~ "**",
  p_adj_BH < 0.05 ~ "*",
  TRUE ~ "ns"
)]

# Restore gene as character for table output
stat_dt[, gene := as.character(gene)]
setcolorder(stat_dt, c("gene","n_tumor","n_normal",
                       "median_tumor","median_normal","log2FC",
                       "p_value","p_adj_BH","sig_label"))
# Order rows in stat_dt by panel order (gene_order_fig)
stat_dt[, gene := factor(gene, levels = gene_order_fig)]
setorder(stat_dt, gene)
stat_dt[, gene := as.character(gene)]

fwrite(stat_dt, TBL_OUT)
cat(sprintf(" [OK] %s\n", TBL_OUT))
cat(" All 18 gene stats:\n")
print(stat_dt)
cat(sprintf(" Significant @ BH<%.2f : %d / %d\n",
            ADJ_P_CUTOFF, sum(stat_dt$p_adj_BH < ADJ_P_CUTOFF, na.rm = TRUE),
            nrow(stat_dt)))
cat("\n")

# 6. Violin + boxplot per gene
cat(">>> STEP 6: Build violin plot (facet_wrap, 18 panels)\n")

# Restore factor for plotting
stat_dt_plot <- copy(stat_dt)
stat_dt_plot[, gene := factor(gene, levels = gene_order_fig)]

# y_pos for sig label: per-gene max value across both groups
ymax_dt <- long_dt[, .(y_max = max(value, na.rm = TRUE)), by = gene]
stat_dt_plot <- merge(stat_dt_plot, ymax_dt, by = "gene")
# y_pos slight offset above max
stat_dt_plot[, y_pos := y_max + 0.4]

group_colors <- c(Normal = "#4DAF4A", Tumor = "#E41A1C")

p <- ggplot(long_dt, aes(x = group, y = value, fill = group)) +
  geom_violin(width = 0.8, alpha = 0.6, color = "black", linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", color = "black",
               outlier.size = 0.4, outlier.alpha = 0.4, linewidth = 0.3) +
  geom_text(data = stat_dt_plot,
            aes(x = 1.5, y = y_pos, label = sig_label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  facet_wrap(~ gene, ncol = 6, nrow = 3, scales = "free_y") +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(
    x = NULL,
    y = "log2(RSEM + 1) expression"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "#EAEAEA", color = NA),
    axis.text.x = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 8),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

# 7. Save (PDF cairo + PNG ragg @600 dpi)
cat(">>> STEP 7: Save figure\n")
W <- 14
H <- 8

pdf_file <- file.path(FIGDIR, paste0(FIG_STEM, ".pdf"))
png_file <- file.path(FIGDIR, paste0(FIG_STEM, ".png"))

ggsave(pdf_file, plot = p, width = W, height = H, device = cairo_pdf)
cat(sprintf(" [OK] %s\n", pdf_file))

ggsave(png_file, plot = p, width = W, height = H, dpi = 600, device = ragg::agg_png)
cat(sprintf(" [OK] %s\n\n", png_file))

# 8. Summary Report
cat("========================================================================\n")
cat(" Stage 6 Tumor vs Normal Expression Summary Report \n")
cat("========================================================================\n")
cat(sprintf(" N_TUMOR (TCGA-ESCA primary) : %d\n", N_TUMOR))
cat(sprintf(" N_NORMAL (GTEx Esophagus - Mucosa) : %d\n", N_NORMAL))
cat(sprintf(" FMG analyzed : %d\n", length(fmg_genes)))
cat(sprintf(" Wilcoxon rank-sum (independent) : applied per gene\n"))
cat(sprintf(" Multiple testing correction : BH (Benjamini-Hochberg)\n"))
cat(sprintf(" Significance cutoff : adj.p < %.3f\n", ADJ_P_CUTOFF))
cat(sprintf(" Significant genes (BH adj.p<%.2f) : %d / %d\n\n",
            ADJ_P_CUTOFF, sum(stat_dt$p_adj_BH < ADJ_P_CUTOFF, na.rm = TRUE), nrow(stat_dt)))

up_n <- sum(stat_dt$log2FC > 0, na.rm = TRUE)
down_n <- sum(stat_dt$log2FC < 0, na.rm = TRUE)
cat(sprintf(" log2FC distribution: up (tumor>normal) = %d\n", up_n))
cat(sprintf(" down (tumor<normal)= %d\n\n", down_n))

cat(" Top 3 |log2FC| genes:\n")
top3_abs <- copy(stat_dt)
top3_abs[, abs_lfc := abs(log2FC)]
setorder(top3_abs, -abs_lfc)
print(head(top3_abs[, .(gene, median_tumor, median_normal, log2FC, p_adj_BH, sig_label)], 3))
cat("\n Top 3 up-regulated (largest +log2FC):\n")
top3_up <- copy(stat_dt); setorder(top3_up, -log2FC)
print(head(top3_up[, .(gene, median_tumor, median_normal, log2FC, p_adj_BH, sig_label)], 3))
cat("\n Top 3 down-regulated (largest -log2FC):\n")
top3_dn <- copy(stat_dt); setorder(top3_dn, log2FC)
print(head(top3_dn[, .(gene, median_tumor, median_normal, log2FC, p_adj_BH, sig_label)], 3))

cat("\n========================================================================\n")

# 9. sessionInfo
sess_file <- file.path(LOGDIR, "sessionInfo/sessionInfo_stage06.txt")
sink(sess_file)
cat(sprintf("Session info captured: %s\n\n", Sys.time()))
print(sessionInfo())
sink()
cat(sprintf("\n>>> sessionInfo saved: %s\n", sess_file))

cat(sprintf("\n=== Stage 6 Pipeline done (%s) ===\n", Sys.time()))

sink(type = "message")
sink()
close(log_con)

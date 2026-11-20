###############################################################################
## Stage 5: TCGA-ESCA CNV Gain/Loss Frequency (18 Primary FMG)
##
## - 4A: GISTIC2.0 thresholded (UCSC Xena TCGA_ESCA_gistic2_thresholded)
## [PMID: 21527027 -- Mermel 2011 Genome Biol]
## - 4B: threshold interpretation
## +2 high-level amplification
## +1 gain
## 0 diploid
## -1 heterozygous loss
## -2 deep deletion
## - 4C: Stacked barplot (gain frequency vs loss frequency per gene)
## [Liu et al. 2024 Figure 2B convention]
##
## Inputs:
## data/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes (24776 x 184)
## data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv
## tables/Tab1_2_FMG_list_primary.csv (18 Primary FMG)
##
## Sample filter:
## sample_type == "Primary Tumor" (clinical) intersect GISTIC sample columns
##
## Outputs (tables/):
## Tab5_TCGA_CNV_frequency.csv
## columns: gene, N_total,
## n_high_amp, high_amp_pct,
## n_gain, gain_pct,
## n_loss, loss_pct,
## n_deep_del, deep_del_pct,
## total_altered_pct
##
## Outputs (figures/, PDF + PNG @ 600 dpi):
## Fig2B_TCGA_CNV_gain_loss.{pdf,png}
## (Hybrid style horizontal stacked bars, 4-color, signed % x-axis)
##
## Reference PMID:
## - 21527027 : Mermel CH et al. 2011 Genome Biol -- GISTIC2.0
##
## Random seed: 1262118388 (Liu et al. 2024)
##
## Notes:
## - Hybrid stacked barplot style:
## * Horizontal bars (coord_flip)
## * 4-color (high_amp/gain/het_loss/deep_del)
## * gain side : gain (near 0) -> high_amp (at bar tip)
## * loss side : loss (near 0) -> deep_del (at bar tip, mirrored)
## * Signed % x-axis labels (-70%, ..., 0%, ..., 70%)
## * Vertical line at x=0
## - Color palette: ColorBrewer RdBu 4-class, color-blind-safe
## high_amp : "#CA0020" dark red
## gain : "#F4A582" light red
## loss : "#92C5DE" light blue
## deep_del : "#0571B0" dark blue
## - Gene order:
## * Figure y-axis (top->bottom) : Stage 4 Tab4 mutation prevalence, descending
## (TP53 at top, FGFR2 at bottom)
## (Liu 2024 Figure 2A/2B convention)
## * Table sort : total_altered_pct, descending
###############################################################################

# 0. Working directory + output
# run from the repository root (see README.md)

INPUT_GISTIC <- "data/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"
INPUT_CLI <- "data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv"
INPUT_FMG <- "tables/Tab1_2_FMG_list_primary.csv"
INPUT_TAB4 <- "tables/Tab4_TCGA_mutation_frequency.csv" # for gene order

TBLDIR <- "tables"
FIGDIR <- "figures"
LOGDIR <- "logs"

for (d in c(TBLDIR, FIGDIR, LOGDIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

FIG_STEM <- "Fig2B_TCGA_CNV_gain_loss"
TBL_OUT <- file.path(TBLDIR, "Tab5_TCGA_CNV_frequency.csv")

# 0.5. Reproducibility
set.seed(1262118388)

# 0.6. Logging start
log_file <- file.path(LOGDIR, "05_diagnostics.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat(sprintf("=== Stage 5 Pipeline start (%s) ===\n\n", Sys.time()))

# 1. load
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ragg)
})

# 2. Load input data
cat(">>> STEP 2: Load inputs\n")
gistic <- fread(INPUT_GISTIC)
cli <- fread(INPUT_CLI)
fmg <- fread(INPUT_FMG)
fmg_genes <- fmg$gene

cat(sprintf(" GISTIC dim : %d genes x %d cols\n", nrow(gistic), ncol(gistic)))
cat(sprintf(" Clinical rows : %d\n", nrow(cli)))
cat(sprintf(" FMG count : %d\n\n", length(fmg_genes)))

# 3. Sample filter: Primary Tumor
cat(">>> STEP 3: Sample filtering\n")
sample_cols_all <- grep("^TCGA-", colnames(gistic), value = TRUE)
cat(sprintf(" GISTIC sample columns : %d\n", length(sample_cols_all)))

primary_tumor_samples <- cli[sample_type == "Primary Tumor", sampleID]
cat(sprintf(" Clinical Primary Tumor : %d\n", length(primary_tumor_samples)))

samples_keep <- intersect(sample_cols_all, primary_tumor_samples)
N_TOTAL <- length(samples_keep)
cat(sprintf(" intersection (analysis) : %d\n", N_TOTAL))
n_dropped <- length(sample_cols_all) - N_TOTAL
if (n_dropped > 0) {
  cat(sprintf(" GISTIC samples dropped (not Primary Tumor): %d\n", n_dropped))
}
cat("\n")

# 4. The 18 FMG rows + value table
cat(">>> STEP 4: Extract 18 FMG CNV rows\n")
gene_col <- colnames(gistic)[1]
setnames(gistic, gene_col, "Gene_Symbol")

fmg_gistic <- gistic[Gene_Symbol %in% fmg_genes]
cat(sprintf(" FMG rows recovered : %d / %d\n", nrow(fmg_gistic), length(fmg_genes)))
missing_fmg <- setdiff(fmg_genes, fmg_gistic$Gene_Symbol)
if (length(missing_fmg) > 0) {
  cat(sprintf(" WARNING: missing in GISTIC: %s\n", paste(missing_fmg, collapse = ", ")))
}

mat_keep <- as.matrix(fmg_gistic[, ..samples_keep])
rownames(mat_keep) <- fmg_gistic$Gene_Symbol
cat(sprintf(" Analysis matrix : %d genes x %d samples\n\n", nrow(mat_keep), ncol(mat_keep)))

# 5. Frequency table
cat(">>> STEP 5: Compute CNV frequencies\n")
count_status <- function(row, status) sum(row == status, na.rm = TRUE)

freq_tbl <- data.table(
  gene = rownames(mat_keep),
  N_total = N_TOTAL,
  n_high_amp = apply(mat_keep, 1, count_status, status = 2L),
  n_gain = apply(mat_keep, 1, count_status, status = 1L),
  n_loss = apply(mat_keep, 1, count_status, status = -1L),
  n_deep_del = apply(mat_keep, 1, count_status, status = -2L)
)

freq_tbl[, high_amp_pct := round(100 * n_high_amp / N_total, 2)]
freq_tbl[, gain_pct := round(100 * n_gain / N_total, 2)]
freq_tbl[, loss_pct := round(100 * n_loss / N_total, 2)]
freq_tbl[, deep_del_pct := round(100 * n_deep_del / N_total, 2)]
freq_tbl[, total_altered_pct := high_amp_pct + gain_pct + loss_pct + deep_del_pct]

setcolorder(freq_tbl, c("gene","N_total",
                        "n_high_amp","high_amp_pct",
                        "n_gain","gain_pct",
                        "n_loss","loss_pct",
                        "n_deep_del","deep_del_pct",
                        "total_altered_pct"))
setorder(freq_tbl, -total_altered_pct)

fwrite(freq_tbl, TBL_OUT)
cat(sprintf(" [OK] %s\n", TBL_OUT))
cat(" Top 5 by total_altered_pct:\n")
print(head(freq_tbl, 5))
cat("\n")

# 6. Stacked barplot
cat(">>> STEP 6: Build stacked barplot\n")

# Gene order for figure x-axis: Stage 4 mutation prevalence (Tab4)
if (!file.exists(INPUT_TAB4)) {
  stop(sprintf("Tab4 not found at %s -- Stage 4 must run before Stage 5.", INPUT_TAB4))
}
tab4 <- fread(INPUT_TAB4)
gene_order_fig <- tab4$gene
cat(sprintf(" Gene order source : %s\n", INPUT_TAB4))
cat(sprintf(" Tab4 gene count : %d\n", length(gene_order_fig)))

# Sanity: gene set match
miss_in_tab4 <- setdiff(freq_tbl$gene, gene_order_fig)
miss_in_freq <- setdiff(gene_order_fig, freq_tbl$gene)
if (length(miss_in_tab4) > 0 || length(miss_in_freq) > 0) {
  cat(sprintf(" WARNING: gene set mismatch -- in CNV not Tab4: %s ; in Tab4 not CNV: %s\n",
              paste(miss_in_tab4, collapse=","), paste(miss_in_freq, collapse=",")))
}
cat(sprintf(" Figure x-axis order (Fig2A unify): %s\n", paste(gene_order_fig, collapse=", ")))

# Long format: gains positive, losses negative
# Stacking control via factor levels:
# ggplot geom_col(position = "stack") puts the first factor level at the base (near 0)
# and stacks the last level outward (at the tip).
# order:
# positive : gain (near 0) -> high_amp (tip)
# negative : loss (near 0) -> deep_del (tip, most negative)
# hence levels = c("gain","high_amp","loss","deep_del")
plot_df <- rbind(
  data.table(gene = freq_tbl$gene, alteration = "gain", percent = freq_tbl$gain_pct),
  data.table(gene = freq_tbl$gene, alteration = "high_amp", percent = freq_tbl$high_amp_pct),
  data.table(gene = freq_tbl$gene, alteration = "loss", percent = -freq_tbl$loss_pct),
  data.table(gene = freq_tbl$gene, alteration = "deep_del", percent = -freq_tbl$deep_del_pct)
)
# Gene factor: rev() so that TP53 is at the top of the horizontal plot (last factor level = top in coord_flip)
plot_df[, gene := factor(gene, levels = rev(gene_order_fig))]
plot_df[, alteration := factor(alteration,
                               levels = c("gain","high_amp","loss","deep_del"))]

cnv_colors <- c(
  high_amp = "#CA0020", # dark red
  gain = "#F4A582", # light red
  loss = "#92C5DE", # light blue
  deep_del = "#0571B0" # dark blue
)

# x axis range: nearest 10% above max abs
max_abs <- max(abs(plot_df$percent), na.rm = TRUE)
x_lim <- ceiling(max_abs / 10) * 10
x_breaks <- seq(-x_lim, x_lim, by = 10)
# Signed % labels: -70%, -60%, ..., 0%, ..., +70%
x_labels <- ifelse(
  x_breaks == 0, "0%",
  ifelse(x_breaks > 0, sprintf("+%d%%", x_breaks),
                          sprintf("%d%%", x_breaks))
)

p <- ggplot(plot_df, aes(x = gene, y = percent, fill = alteration)) +
  geom_col(width = 0.7, color = NA) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
  scale_fill_manual(
    values = cnv_colors,
    breaks = c("high_amp","gain","loss","deep_del"),
    labels = c("High-level amplification (+2)",
               "Gain (+1)",
               "Heterozygous loss (-1)",
               "Deep deletion (-2)"),
    name = NULL
  ) +
  scale_y_continuous(limits = c(-x_lim, x_lim),
                     breaks = x_breaks,
                     labels = x_labels,
                     expand = expansion(mult = c(0.02, 0.02))) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Frequency (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    legend.box = "horizontal"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# 7. Save (PDF cairo + PNG ragg @600 dpi)
cat(">>> STEP 7: Save figure\n")
# Hybrid horizontal layout, tall (8 in wide x 10 in tall)
W <- 8
H <- 10

pdf_file <- file.path(FIGDIR, paste0(FIG_STEM, ".pdf"))
png_file <- file.path(FIGDIR, paste0(FIG_STEM, ".png"))

ggsave(pdf_file, plot = p, width = W, height = H, device = cairo_pdf)
cat(sprintf(" [OK] %s\n", pdf_file))

ggsave(png_file, plot = p, width = W, height = H, dpi = 600, device = ragg::agg_png)
cat(sprintf(" [OK] %s\n\n", png_file))

# 8. Summary Report
cat("========================================================================\n")
cat(" Stage 5 TCGA-ESCA CNV Gain/Loss Summary Report \n")
cat("========================================================================\n")
cat(sprintf(" GISTIC sample columns : %d\n", length(sample_cols_all)))
cat(sprintf(" Primary Tumor analysis : %d\n", N_TOTAL))
cat(sprintf(" FMG analyzed : %d\n\n", length(fmg_genes)))

cat(" --- 18 FMG CNV frequency (sorted by total_altered_pct) ---\n")
print(freq_tbl)
cat("\n========================================================================\n")

# 9. sessionInfo
sess_file <- file.path(LOGDIR, "sessionInfo/sessionInfo_stage05.txt")
sink(sess_file)
cat(sprintf("Session info captured: %s\n\n", Sys.time()))
print(sessionInfo())
sink()
cat(sprintf("\n>>> sessionInfo saved: %s\n", sess_file))

cat(sprintf("\n=== Stage 5 Pipeline done (%s) ===\n", Sys.time()))

sink(type = "message")
sink()
close(log_con)

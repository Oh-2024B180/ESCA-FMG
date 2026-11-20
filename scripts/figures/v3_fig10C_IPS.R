###############################################################################
## Fig 10C remake (v3) -- IPS metrics boxplot, Fig10A style
##
## Figure settings (2026-05-10):
##   - jitter dots removed
##   - raw P numerical labels removed
##   - significance shown only by * stars
##   - larger fonts
##   - 92% data-range star positioning (Taselisib style)
##   - black panel + strip borders
##
## Inputs:
##   data/expr_ESCA_tumor.csv               (gene x sample log2(RSEM+1))
##   tables/Tab10_risk_scores_per_sample.csv (sampleID, RiskGroup)
##   tables/Tab14B_IPS_metrics_test.csv      (raw P per metric)
##   IOBR::IPS_calculation()                 (5 metrics: IPS, MHC, EC, SC, CP, AZ)
##
## Output:
##   figures/remake/Fig10C_IPS_metrics_boxplot_v3.{pdf,png}
##   figures/remake/v3_fig10C_diagnostics.log
##   cache/v3_fig10C_ips_long.rds            (per-sample IPS cache)
##
## Random seed: 1262118388 (CLAUDE.md §15)
###############################################################################

# 0. setup ---------------------------------------------------------------------
# run from the repository root (see README.md)

OUT_DIR  <- "figures/remake"
DATA_DIR <- "data"
TBL_DIR  <- "tables"
CACHE    <- "cache"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(CACHE))   dir.create(CACHE,   recursive = TRUE, showWarnings = FALSE)

set.seed(1262118388)

log_file <- file.path(OUT_DIR, "v3_fig10C_diagnostics.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
cat(sprintf("=== v3 Fig10C remake start (%s) ===\n\n", Sys.time()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(ragg)
})

COL_LOW  <- "#377EB8"
COL_HIGH <- "#E41A1C"

THEME_BASE <- theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank(),
    axis.text  = element_text(size = 11, colour = "black"),
    axis.title = element_text(size = 12),
    strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.6),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    panel.spacing = unit(0.6, "lines")
  )

sig_from_p <- function(p) {
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.05]   <- "*"
  out[!is.na(p) & p < 0.01]   <- "**"
  out[!is.na(p) & p < 0.001]  <- "***"
  out
}

# 1. load risk scores ---------------------------------------------------------
cat("[1] Load RS\n")
rs_dt <- fread(file.path(TBL_DIR, "Tab10_risk_scores_per_sample.csv"))
rs_dt[, RiskGroup := factor(RiskGroup, levels = c("Low", "High"))]
cat(sprintf("  RS samples: Low=%d, High=%d\n",
            sum(rs_dt$RiskGroup == "Low"), sum(rs_dt$RiskGroup == "High")))

# 2. compute IPS via IOBR (cached) --------------------------------------------
ips_cache <- file.path(CACHE, "v3_fig10C_ips_long.rds")
if (file.exists(ips_cache)) {
  cat("\n[2] Load IPS from cache\n")
  ips_dt <- readRDS(ips_cache)
} else {
  cat("\n[2] Compute IPS via IOBR::IPS_calculation()\n")
  suppressPackageStartupMessages(library(IOBR))

  expr_dt <- fread(file.path(DATA_DIR, "expr_ESCA_tumor.csv"))
  g_col <- colnames(expr_dt)[1]; setnames(expr_dt, g_col, "gene_id")
  expr_mat <- as.matrix(expr_dt[, -1])
  rownames(expr_mat) <- expr_dt$gene_id

  # subset to RS samples only
  expr_mat <- expr_mat[, intersect(colnames(expr_mat), rs_dt$sampleID), drop = FALSE]
  cat(sprintf("  expression matrix: %d genes x %d samples\n",
              nrow(expr_mat), ncol(expr_mat)))

  ips_res <- IPS_calculation(eset = expr_mat, plot = FALSE)
  ips_dt <- as.data.table(ips_res)
  if (!"sample" %in% colnames(ips_dt)) {
    if (!is.null(rownames(ips_res))) {
      ips_dt[, sample := rownames(ips_res)]
    } else {
      ips_dt[, sample := colnames(expr_mat)[seq_len(.N)]]
    }
  }
  setnames(ips_dt, "sample", "sampleID")
  saveRDS(ips_dt, ips_cache)
  cat(sprintf("  cached: %s\n", ips_cache))
}
cat("  IPS cols: "); cat(colnames(ips_dt), sep = ", "); cat("\n")

# 3. tidy + merge --------------------------------------------------------------
cat("\n[3] Long format + merge with RS\n")
test_dt <- fread(file.path(TBL_DIR, "Tab14B_IPS_metrics_test.csv"))
metric_levels <- test_dt$metric  # use the order from Tab14B (raw p ascending)
cat(sprintf("  metrics (order = Tab14B): %s\n",
            paste(metric_levels, collapse = ", ")))

# keep only metric columns + sampleID
keep_cols <- c("sampleID", intersect(metric_levels, colnames(ips_dt)))
ips_sub <- ips_dt[, ..keep_cols]
ips_long <- melt(ips_sub, id.vars = "sampleID",
                 variable.name = "metric", value.name = "value")
ips_long <- merge(ips_long, rs_dt[, .(sampleID, RiskGroup)], by = "sampleID")
ips_long[, metric := factor(metric, levels = metric_levels)]
ips_long[, value  := as.numeric(value)]
ips_long <- ips_long[!is.na(value)]

# 4. sig labels ----------------------------------------------------------------
cat("\n[4] Significance labels (raw P)\n")
sig_10C <- test_dt[, .(metric, raw_p = p_value)]
sig_10C[, sig := sig_from_p(raw_p)]
sig_10C[, metric := factor(metric, levels = metric_levels)]

y_top <- ips_long[, .(ymax = max(value, na.rm = TRUE),
                      ymin = min(value, na.rm = TRUE)),
                  by = metric]
sig_10C <- merge(sig_10C, y_top, by = "metric")
sig_10C[, ylabel := ymin + (ymax - ymin) * 0.92]
print(sig_10C)

# 5. plot ----------------------------------------------------------------------
cat("\n[5] Plot\n")
fig10C <- ggplot(ips_long, aes(x = RiskGroup, y = value, fill = RiskGroup)) +
  geom_violin(alpha = 0.6, colour = "black", linewidth = 0.4, trim = TRUE) +
  geom_boxplot(width = 0.18, fill = "white", colour = "black",
               outlier.size = 0.4, linewidth = 0.4) +
  geom_text(data = sig_10C,
            aes(x = 1.5, y = ylabel, label = sig),
            inherit.aes = FALSE, size = 6, fontface = "bold") +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(Low = COL_LOW, High = COL_HIGH),
                    name = "Risk group") +
  labs(x = NULL, y = "IPS score") +
  THEME_BASE

# save
pdf_path <- file.path(OUT_DIR, "Fig10C_IPS_metrics_boxplot_v3.pdf")
png_path <- file.path(OUT_DIR, "Fig10C_IPS_metrics_boxplot_v3.png")
ggsave(pdf_path, fig10C, width = 9, height = 6, device = cairo_pdf)
agg_png(png_path, width = 9, height = 6, units = "in", res = 600)
  print(fig10C)
invisible(dev.off())
cat(sprintf("  saved: %s\n", basename(pdf_path)))

cat(sprintf("\n=== v3 Fig10C remake done (%s) ===\n", Sys.time()))
sink(type = "message"); sink(); close(log_con)

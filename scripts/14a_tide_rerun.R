###############################################################################
## Stage 14A re-run, part 2 of 2 (2026-10-19): risk-group tests and the TIDE panel on the corrected TIDE output
##
## Part 1 (scripts/14a_tide_rerun.py) re-ran TIDEpy 1.3.9 on a singly logged, gene-centred input
## (data/esca_tide_output_hgnc_centered.txt). This part repeats, unchanged, the code that used the April output:
##   - tests: stage 14A of scripts/14_treatment_response_prediction.R ([14A.2]-[14A.4]: Wilcoxon for the 11
##     continuous metrics, chi-square or Fisher for Responder / No.benefits / CTL.flag); only the input file differs
##   - panel: the Fig10A block of figures/remake/v3_fig8_10_remake.R (manuscript Figure 11a), same size and style
## Outputs tables/Tab14A_TIDE_metrics_test.csv, tables/Tab14A_responder_chisq.csv,
##         figures/remake/Fig10A_TIDE_metrics_boxplot_v3.{pdf,png}, logs/14a_tide_rerun_R.log
###############################################################################
# run from the repository root (see README.md)
INPUT_TIDE <- "data/esca_tide_output_hgnc_centered.txt"
INPUT_RS <- "tables/Tab10_risk_scores_per_sample.csv"
INPUT_CLI <- "data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv"
TBLDIR <- "tables"
OUT_DIR <- "figures/remake"
P_CUT <- 0.05
set.seed(1262118388)

log_con <- file("logs/14a_tide_rerun_R.log", open = "wt")
sink(log_con, split = TRUE)
cat(sprintf("=== Stage 14A re-run, part 2 (%s) ===\n", Sys.time()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ragg)
})
cat(sprintf(" R %s, data.table %s, ggplot2 %s, ragg %s\n", paste(R.Version()[c("major", "minor")], collapse = "."),
            packageVersion("data.table"), packageVersion("ggplot2"), packageVersion("ragg")))

# ---------------------------------------------------------------- stage 14, [Common] (unchanged)
rs_dt <- fread(INPUT_RS)
rs_dt[, RiskGroup := factor(RiskGroup, levels = c("Low","High"))]
cli_dt <- fread(INPUT_CLI)
cli_dt[, histology_simple := fcase(
  histological_type == "Esophagus Squamous Cell Carcinoma", "ESCC",
  histological_type == "Esophagus Adenocarcinoma, NOS", "EAC",
  default = NA_character_
)]
rs_cli <- merge(rs_dt[, .(sampleID, RiskGroup, RiskScore)],
                cli_dt[, .(sampleID, histology_simple)],
                by = "sampleID", all.x = TRUE)

wilcox_lh <- function(values, group) {
  ok <- !is.na(values) & is.finite(values) & !is.na(group)
  values <- values[ok]; group <- group[ok]
  if (length(unique(group)) < 2) {
    return(list(W = NA_real_, p = NA_real_,
                med_low = NA_real_, med_high = NA_real_,
                n_low = sum(group == "Low"), n_high = sum(group == "High")))
  }
  wt <- wilcox.test(values ~ group, exact = FALSE)
  list(W = unname(wt$statistic), p = wt$p.value,
       med_low = median(values[group == "Low"], na.rm = TRUE),
       med_high = median(values[group == "High"], na.rm = TRUE),
       n_low = sum(group == "Low"),
       n_high = sum(group == "High"))
}

# ---------------------------------------------------------------- stage 14A (unchanged apart from INPUT_TIDE)
tide_dt <- fread(INPUT_TIDE)
cn1 <- colnames(tide_dt)[1]
setnames(tide_dt, cn1, "sampleID")
setnames(tide_dt, gsub(" ", ".", colnames(tide_dt), fixed = TRUE))
TIDE_METRICS <- c("TIDE","Dysfunction","Exclusion","MDSC","CAF","TAM.M2",
                  "IFNG","MSI.Score","CD274","CD8","CTL")
TIDE_BOOLS <- c("Responder","No.benefits","CTL.flag")
stopifnot(all(TIDE_METRICS %in% colnames(tide_dt)), all(TIDE_BOOLS %in% colnames(tide_dt)))
tide_dt[, (TIDE_BOOLS) := lapply(.SD, function(x) {
  out <- toupper(as.character(x))
  ifelse(out %in% c("TRUE","T","YES","1"), TRUE,
  ifelse(out %in% c("FALSE","F","NO","0"), FALSE, NA))
}), .SDcols = TIDE_BOOLS]
tide_rs <- merge(rs_cli, tide_dt, by = "sampleID")
N_14A <- nrow(tide_rs)
cat(sprintf(" N (TIDE x RS): %d\n", N_14A))
stopifnot(N_14A == 181)

test_a_rows <- lapply(TIDE_METRICS, function(m) {
  res <- wilcox_lh(tide_rs[[m]], tide_rs$RiskGroup)
  data.table(metric = m, n_low = res$n_low, n_high = res$n_high,
             median_low = res$med_low, median_high = res$med_high,
             W = res$W, p_value = res$p)
})
test_a <- rbindlist(test_a_rows)
test_a[, p_adj_BH := p.adjust(p_value, method = "BH")]
test_a[, sig := fifelse(p_value < 0.001, "***",
                  fifelse(p_value < 0.01, "**",
                  fifelse(p_value < 0.05, "*", "")))]
setorder(test_a, p_value)
fwrite(test_a, file.path(TBLDIR, "Tab14A_TIDE_metrics_test.csv"))
print(test_a)
cat(sprintf(" raw p<%.2f: %d / %d metrics\n", P_CUT, sum(test_a$p_value < P_CUT, na.rm = TRUE), length(TIDE_METRICS)))

chisq_rows <- lapply(TIDE_BOOLS, function(b) {
  tab <- table(RiskGroup = tide_rs$RiskGroup, Flag = tide_rs[[b]])
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(data.table(flag = b, n_total = sum(tab),
                      method = NA, statistic = NA, p_value = NA,
                      n_low_T = NA, n_low_F = NA,
                      n_high_T = NA, n_high_F = NA))
  }
  exp_min <- min(chisq.test(tab, simulate.p.value = FALSE)$expected)
  use_fisher <- exp_min < 5
  res <- if (use_fisher) {
    f <- fisher.test(tab); list(stat = NA_real_, p = f$p.value, m = "fisher")
  } else {
    c2 <- chisq.test(tab, correct = FALSE)
    list(stat = unname(c2$statistic), p = c2$p.value, m = "chisq")
  }
  ncol_t <- if ("TRUE" %in% colnames(tab)) tab[, "TRUE", drop = FALSE] else NA
  ncol_f <- if ("FALSE" %in% colnames(tab)) tab[, "FALSE", drop = FALSE] else NA
  data.table(
    flag = b,
    n_total = sum(tab),
    method = res$m,
    statistic = res$stat,
    p_value = res$p,
    n_low_T = if (!is.na(ncol_t[1])) tab["Low", "TRUE"] else 0L,
    n_low_F = if (!is.na(ncol_f[1])) tab["Low", "FALSE"] else 0L,
    n_high_T = if (!is.na(ncol_t[1])) tab["High", "TRUE"] else 0L,
    n_high_F = if (!is.na(ncol_f[1])) tab["High", "FALSE"] else 0L
  )
})
chisq_dt <- rbindlist(chisq_rows)
chisq_dt[, p_adj_BH := p.adjust(p_value, method = "BH")]
chisq_dt[, sig := fifelse(p_value < 0.05, "*", "")]
fwrite(chisq_dt, file.path(TBLDIR, "Tab14A_responder_chisq.csv"))
print(chisq_dt)

# ---------------------------------------------------------------- Figure 11a (v3_fig8_10_remake.R, Fig10A block)
COL_LOW  <- "#377EB8"
COL_HIGH <- "#E41A1C"
THEME_BASE <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    axis.text  = element_text(colour = "black", size = 11),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10)
  )
sig_from_p <- function(p) {
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.05]   <- "*"
  out[!is.na(p) & p < 0.01]   <- "**"
  out[!is.na(p) & p < 0.001]  <- "***"
  out
}
save_figure <- function(plot_obj, stem, w_in = 8, h_in = 8) {
  pdf_path <- file.path(OUT_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(OUT_DIR, paste0(stem, ".png"))
  ggsave(pdf_path, plot_obj, width = w_in, height = h_in, device = cairo_pdf)
  agg_png(png_path, width = w_in, height = h_in, units = "in", res = 600)
    print(plot_obj)
  invisible(dev.off())
  cat(sprintf("  saved: %s (%.1f x %.1f in)\n", stem, w_in, h_in))
}

tide_raw <- fread(INPUT_TIDE, header = TRUE)
setnames(tide_raw, 1, "sampleID")
tide_test <- fread(file.path(TBLDIR, "Tab14A_TIDE_metrics_test.csv"))
tide_metrics <- tide_test$metric
metric_to_col <- c(
  "MSI.Score"   = "MSI Score",
  "Dysfunction" = "Dysfunction",
  "Exclusion"   = "Exclusion",
  "CAF"         = "CAF",
  "MDSC"        = "MDSC",
  "CD274"       = "CD274",
  "CD8"         = "CD8",
  "TAM.M2"      = "TAM M2",
  "TIDE"        = "TIDE",
  "IFNG"        = "IFNG",
  "CTL"         = "CTL"
)
tide_long <- rbindlist(lapply(tide_metrics, function(m) {
  col <- metric_to_col[m]
  if (is.na(col) || !col %in% colnames(tide_raw)) return(NULL)
  d <- data.table(sampleID = tide_raw$sampleID,
                  metric = m,
                  value = suppressWarnings(as.numeric(tide_raw[[col]])))
  d[!is.na(value)]
}))
tide_long <- merge(tide_long, rs_dt[, .(sampleID, RiskGroup)], by = "sampleID")
tide_long[, metric := factor(metric, levels = tide_metrics)]
sig_10A <- tide_test[, .(metric, raw_p = p_value)]
sig_10A[, sig := sig_from_p(raw_p)]
sig_10A[, metric := factor(metric, levels = tide_metrics)]
y_top <- tide_long[, .(ymax = max(value, na.rm = TRUE),
                       ymin = min(value, na.rm = TRUE)), by = metric]
sig_10A <- merge(sig_10A, y_top, by = "metric")
sig_10A[, ylabel := ymin + (ymax - ymin) * 0.92]

fig10A <- ggplot(tide_long, aes(x = RiskGroup, y = value, fill = RiskGroup)) +
  geom_boxplot(width = 0.55, outlier.size = 0.6, alpha = 0.85, colour = "black") +
  geom_text(data = sig_10A,
            aes(x = 1.5, y = ylabel, label = sig),
            inherit.aes = FALSE, size = 6, fontface = "bold") +
  facet_wrap(~ metric, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(Low = COL_LOW, High = COL_HIGH),
                    name = "Risk group") +
  labs(x = NULL, y = "Score") +
  THEME_BASE +
  theme(
    legend.position = "top",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.6),
    strip.text = element_text(size = 11, face = "bold"),
    panel.spacing = unit(0.6, "lines")
  )
save_figure(fig10A, "Fig10A_TIDE_metrics_boxplot_v3", w_in = 10, h_in = 8)
cat(sprintf("=== end (%s) ===\n", Sys.time()))
sink()
close(log_con)

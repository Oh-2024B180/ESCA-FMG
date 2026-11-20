###############################################################################
## Stage 14: Treatment Response Prediction
##
## NOTE: 14A below ran on data/esca_tide_output_hgnc_self.txt, which is NOT a TIDE web output: it was made
## on 2026-04-19 with TIDEpy 1.3.9 locally ("-c Other --ignore_norm") from an input logged twice (log2(x+1) of
## values already log2(norm_count+1)). 14A was re-run with a singly logged, gene-centred input by
## scripts/14a_tide_rerun.py and scripts/14a_tide_rerun.R, which now write the Tab14A_* tables and the Figure 11a
## panel; this file is kept as it ran in May (14B IPS is unchanged).
## NOTE: 14C below passes GDSC2_Res (natural-log IC50) to calcPhenotype() without exp(); the
## oncoPredict vignette exponentiates it first. 14C was re-run that way (scripts/14b_oncopredict_predictions.R, run B)
## and scripts/14c_oncopredict_vignette.R now writes the Tab14C_* tables and the Figure 11c panel.
##
## input:
## - data/esca_tide_output_hgnc_self.txt (14A TIDE output, n=181; see NOTE above)
## - data/expr_ESCA_tumor.csv (14B/14C input)
## - data/GDSC2_Expr (RMA Normalized and Log Transformed).rds (14C training)
## - data/GDSC2_Res.rds (14C training IC50)
## - data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv (for subtype stratification)
## - tables/Tab10_risk_scores_per_sample.csv (RiskGroup join key)
##
## output:
## - tables/Tab14A_TIDE_metrics_test.csv
## - tables/Tab14A_responder_chisq.csv
## - tables/Tab14B_IPS_metrics_test.csv
## - tables/Tab14C_oncoPredict_IC50_full.csv
## - tables/Tab14C_oncoPredict_drug_test.csv
## - tables/Tab14C_oncoPredict_top30_drugs.csv
## - figures/Fig10A_TIDE_metrics_boxplot.{pdf,png}
## - figures/Fig10B_TIDE_responder_proportion.{pdf,png}
## - figures/Fig10C_IPS_metrics_boxplot.{pdf,png}
## - figures/Fig10D_oncoPredict_top_drugs_boxplot.{pdf,png}
## - logs/14_treatment_response_diagnostics.txt
##
## cutoff:
## - 14A: Wilcoxon Low-RS vs High-RS, raw p < 0.05 (per Liu 2024 §2.13); Responder chi-square (Fisher fallback if expected < 5)
## - 14B: IOBR::IPS_calculation 5 metrics (IPS, MHC, EC, SC, CP); Wilcoxon + raw p < 0.05
## - 14C: oncoPredict::calcPhenotype with GDSC2 training; Wilcoxon + raw p < 0.05; top 30 by raw p
## - BH-adjusted p values are retained as supplementary columns in all output CSVs.
##
## PMID:
## - 30127393: Jiang et al. 2018 Nat Med TIDE
## - 28052254: Charoentong et al. 2017 Cell Rep IPS/TCIA framework
## (the TCIA portal was not used; IPS from IOBR::IPS_calculation)
## - 34276676: Zeng et al. 2021 Front Immunol IOBR R package
## - 34260682: Maeser et al. 2021 Brief Bioinform oncoPredict
##
## Random seed: 1262118388
##
## Notes:
## - 14C oncoPredict::calcPhenotype run time: ~5-10 min (GDSC2, 198 drugs x 181 samples).
## memory use ~4 GB.
###############################################################################

# run from the repository root (see README.md)

INPUT_TIDE <- "data/esca_tide_output_hgnc_self.txt"
INPUT_EXPR <- "data/expr_ESCA_tumor.csv"
INPUT_RS <- "tables/Tab10_risk_scores_per_sample.csv"
INPUT_CLI <- "data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv"
INPUT_GDSC_EXPR <- "data/GDSC2_Expr (RMA Normalized and Log Transformed).rds"
INPUT_GDSC_RES <- "data/GDSC2_Res.rds"

TBLDIR <- "tables"
FIGDIR <- "figures"
LOGDIR <- "logs"
CACHE <- "cache"
for (d in c(TBLDIR, FIGDIR, LOGDIR, CACHE)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

P_CUT <- 0.05  # raw p significance cutoff (Liu 2024 §2.13)
TOP_N_DRUG <- 30 # Stage 14C top drugs by raw p

RISK_COL <- c(Low = "#377EB8", High = "#E41A1C")

set.seed(1262118388)
t_start <- Sys.time()

log_file <- file.path(LOGDIR, "14_diagnostics.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE); sink(log_con, type = "message")
cat(sprintf("=== Stage 14 start (%s) ===\n\n", Sys.time()))

cat("[0] Loading packages\n")
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ragg)
  library(IOBR)
})

cat(sprintf(" R version : %s\n", paste(R.Version()[c("major","minor")], collapse = ".")))
cat(sprintf(" IOBR : %s\n", packageVersion("IOBR")))

# Common: load RS + sample alignment
cat("\n[Common] Load RS + clinical (subtype) for all sub-stages\n")
rs_dt <- fread(INPUT_RS)
rs_dt[, RiskGroup := factor(RiskGroup, levels = c("Low","High"))]
samples_rs <- rs_dt$sampleID
cat(sprintf(" RS samples (Stage 10): %d (Low=%d, High=%d)\n",
            nrow(rs_dt),
            sum(rs_dt$RiskGroup == "Low"),
            sum(rs_dt$RiskGroup == "High")))

cli_dt <- fread(INPUT_CLI)
cli_dt[, histology_simple := fcase(
  histological_type == "Esophagus Squamous Cell Carcinoma", "ESCC",
  histological_type == "Esophagus Adenocarcinoma, NOS", "EAC",
  default = NA_character_
)]
rs_cli <- merge(rs_dt[, .(sampleID, RiskGroup, RiskScore)],
                cli_dt[, .(sampleID, histology_simple)],
                by = "sampleID", all.x = TRUE)
cat(sprintf(" histology_simple: ESCC=%d, EAC=%d, NA=%d\n",
            sum(rs_cli$histology_simple == "ESCC", na.rm = TRUE),
            sum(rs_cli$histology_simple == "EAC", na.rm = TRUE),
            sum(is.na(rs_cli$histology_simple))))

# Helper: Wilcoxon Low vs High with NA-safe handling
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

# Helper: ggplot boxplot panel with Wilcoxon p
make_box_panel <- function(plot_dt, metric_levels, ncol = 4,
                           title_text = NULL, subtitle_text = NULL,
                           jitter = TRUE) {
  plot_dt[, metric := factor(metric, levels = metric_levels)]
  ann_dt <- plot_dt[, .(p = wilcox.test(value ~ RiskGroup, exact = FALSE)$p.value,
                        y = max(value, na.rm = TRUE) +
                            0.08 * (max(value, na.rm = TRUE) -
                                     min(value, na.rm = TRUE) + 1e-6)),
                    by = metric]
  ann_dt[, lbl := sprintf("p = %.2g", p)]

  p <- ggplot(plot_dt, aes(x = RiskGroup, y = value, fill = RiskGroup)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, color = "black", linewidth = 0.4)
  if (jitter) {
    p <- p + geom_jitter(width = 0.15, size = 0.7, alpha = 0.45, color = "grey25")
  }
  p +
    geom_text(data = ann_dt, aes(x = 1.5, y = y, label = lbl),
              inherit.aes = FALSE, size = 3.0, fontface = "plain") +
    facet_wrap(~ metric, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = RISK_COL, name = "RiskGroup") +
    labs(x = NULL, y = "Score") +
    theme_bw(base_size = 11) +
    theme(legend.position = "top",
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(face = "bold"))
}

###############################################################################
# 14A: TIDE immunotherapy response prediction (precomputed TIDE output; see the NOTE in the header)
###############################################################################
cat("\n========================================================================\n")
cat(" Stage 14A: TIDE \n")
cat("========================================================================\n")

if (file.exists(INPUT_TIDE)) {   # the April TIDE output is not distributed; 14A was superseded by scripts/14a_tide_rerun.*
cat("[14A.1] Load TIDE output\n")
tide_dt <- fread(INPUT_TIDE)
# fread leaves the first (unnamed) column as 'V1'; read.delim renames it 'X'.
# Rename whichever shows up. Also normalize space dot in remaining columns.
cn1 <- colnames(tide_dt)[1]
setnames(tide_dt, cn1, "sampleID")
setnames(tide_dt, gsub(" ", ".", colnames(tide_dt), fixed = TRUE))
cat(sprintf(" TIDE rows: %d, cols: %d\n", nrow(tide_dt), ncol(tide_dt)))
cat(sprintf(" cols (normalized): %s\n", paste(colnames(tide_dt), collapse = ", ")))

# Continuous metrics (numeric)
TIDE_METRICS <- c("TIDE","Dysfunction","Exclusion","MDSC","CAF","TAM.M2",
                  "IFNG","MSI.Score","CD274","CD8","CTL")
# Boolean / categorical responder flags
TIDE_BOOLS <- c("Responder","No.benefits","CTL.flag")

stopifnot(all(TIDE_METRICS %in% colnames(tide_dt)))
stopifnot(all(TIDE_BOOLS %in% colnames(tide_dt)))

# Coerce booleans (TIDE writes "True"/"False" strings)
tide_dt[, (TIDE_BOOLS) := lapply(.SD, function(x) {
  out <- toupper(as.character(x))
  ifelse(out %in% c("TRUE","T","YES","1"), TRUE,
  ifelse(out %in% c("FALSE","F","NO","0"), FALSE, NA))
}), .SDcols = TIDE_BOOLS]

cat("[14A.2] Join TIDE with RiskGroup (181-sample 3-way intersect)\n")
tide_rs <- merge(rs_cli, tide_dt, by = "sampleID")
N_14A <- nrow(tide_rs)
cat(sprintf(" N (TIDE RS): %d\n", N_14A))
stopifnot(N_14A > 50)

cat("[14A.3] Wilcoxon test for 11 continuous metrics\n")
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
cat(" Tab14A_TIDE_metrics_test.csv\n")
print(test_a)
n_sig_a <- sum(test_a$p_value < P_CUT, na.rm = TRUE)
cat(sprintf(" raw p<%.2f: %d / %d metrics\n", P_CUT, n_sig_a, length(TIDE_METRICS)))

cat("\n[14A.4] Responder/No.benefits/CTL.flag chi-square (Fisher fallback)\n")
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
  # Ensure both columns exist for cleaner reporting
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
cat(" Tab14A_responder_chisq.csv\n")
print(chisq_dt)

cat("\n[14A.5] Render Fig10A TIDE 11-metric boxplot\n")
plot_a <- melt(tide_rs[, c("sampleID","RiskGroup", TIDE_METRICS), with = FALSE],
               id.vars = c("sampleID","RiskGroup"),
               variable.name = "metric", value.name = "value")
plot_a <- plot_a[!is.na(value) & is.finite(value)]

p10A <- make_box_panel(plot_a, metric_levels = test_a$metric,
                        ncol = 4,
                        title_text = "TIDE metrics by RiskGroup (Stage 14A)",
                        subtitle_text = sprintf(
                          "n = %d (Low = %d, High = %d); Wilcoxon raw p annotated; %d / %d metric raw p<%.2f",
                          N_14A,
                          sum(tide_rs$RiskGroup == "Low"),
                          sum(tide_rs$RiskGroup == "High"),
                          n_sig_a, length(TIDE_METRICS), P_CUT))

W10A <- 12; H10A <- 9
ggsave(file.path(FIGDIR, "Fig10A_TIDE_metrics_boxplot.pdf"),
       plot = p10A, width = W10A, height = H10A, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig10A_TIDE_metrics_boxplot.png"),
       plot = p10A, width = W10A, height = H10A, dpi = 600, device = ragg::agg_png)
cat(" [OK] Fig10A_TIDE_metrics_boxplot.{pdf,png}\n")

cat("\n[14A.6] Render Fig10B Responder/No.benefits proportion stacked bar\n")
prop_rows <- lapply(TIDE_BOOLS, function(b) {
  d <- tide_rs[!is.na(get(b)),
               .(n = .N), by = .(RiskGroup, flag_value = get(b))]
  d[, total := sum(n), by = RiskGroup]
  d[, prop := n / total]
  d[, flag := b]
  d[, flag_value := factor(flag_value, levels = c(FALSE, TRUE),
                            labels = c("No","Yes"))]
  d
})
prop_dt <- rbindlist(prop_rows)
chisq_lookup <- setNames(chisq_dt$p_value, chisq_dt$flag)
prop_dt[, flag_label := sprintf("%s (p = %.2g)",
                                  flag, chisq_lookup[flag])]
flag_levels <- sapply(TIDE_BOOLS, function(b)
  sprintf("%s (p = %.2g)", b, chisq_lookup[b]))
prop_dt[, flag_label := factor(flag_label, levels = flag_levels)]

p10B <- ggplot(prop_dt, aes(x = RiskGroup, y = prop, fill = flag_value)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%d", n)),
            position = position_stack(vjust = 0.5), size = 3, color = "white",
            fontface = "bold") +
  scale_fill_manual(values = c("No" = "grey70", "Yes" = "#E41A1C"),
                    name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.04))) +
  facet_wrap(~ flag_label, ncol = 3) +
  labs(x = NULL, y = "Proportion of samples") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold"))

W10B <- 9; H10B <- 5
ggsave(file.path(FIGDIR, "Fig10B_TIDE_responder_proportion.pdf"),
       plot = p10B, width = W10B, height = H10B, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig10B_TIDE_responder_proportion.png"),
       plot = p10B, width = W10B, height = H10B, dpi = 600, device = ragg::agg_png)
cat(" [OK] Fig10B_TIDE_responder_proportion.{pdf,png}\n")
} else {
  cat("[14A] skipped:", INPUT_TIDE, "not found (14A was superseded by scripts/14a_tide_rerun.*)\n")
  n_sig_a <- NA_integer_; TIDE_METRICS <- character(0)
  chisq_dt <- data.table(flag = character(), method = character(), p_value = numeric(), p_adj_BH = numeric(), sig = character())
}

###############################################################################
# 14B: IPS via IOBR::IPS_calculation
###############################################################################
cat("\n========================================================================\n")
cat(" Stage 14B: IPS \n")
cat("========================================================================\n")

cat("[14B.1] Load expression matrix\n")
expr_dt <- fread(INPUT_EXPR)
gene_col <- colnames(expr_dt)[1]
expr_genes <- expr_dt[[gene_col]]
expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- expr_genes
# de-dup gene names if any
dup_idx <- duplicated(rownames(expr_mat))
if (any(dup_idx)) {
  cat(sprintf(" Dropping %d duplicate gene rows\n", sum(dup_idx)))
  expr_mat <- expr_mat[!dup_idx, , drop = FALSE]
}
cat(sprintf(" expr_mat: %d genes x %d samples\n", nrow(expr_mat), ncol(expr_mat)))

cat("[14B.2] IOBR::IPS_calculation()\n")
ips_res <- IPS_calculation(eset = expr_mat, plot = FALSE)
cat(sprintf(" IPS result: %d rows x %d cols\n", nrow(ips_res), ncol(ips_res)))
cat(sprintf(" cols: %s\n", paste(colnames(ips_res), collapse = ", ")))
print(head(ips_res))

# Locate sample identifier column
ips_dt <- as.data.table(ips_res)
if (!"sample" %in% colnames(ips_dt)) {
  if (!is.null(rownames(ips_res))) {
    ips_dt[, sample := rownames(ips_res)]
  } else {
    ips_dt[, sample := colnames(expr_mat)[seq_len(.N)]]
  }
}
# IOBR returns: ID, MHC, EC, SC, CP, AZ, IPS (depending on version)
IPS_METRICS <- intersect(c("IPS","MHC","EC","SC","CP","AZ"),
                          colnames(ips_dt))
cat(sprintf(" Detected IPS metrics: %s\n",
            paste(IPS_METRICS, collapse = ", ")))
stopifnot(length(IPS_METRICS) >= 4) # at least IPS + 3 components

cat("[14B.3] Join with RS + Wilcoxon Low vs High\n")
# IOBR may name the sample column differently; accept either 'sample' or 'ID'
sample_col_ips <- if ("sample" %in% colnames(ips_dt)) "sample" else "ID"
ips_dt[, sampleID := get(sample_col_ips)]
ips_rs <- merge(rs_cli, ips_dt[, c("sampleID", IPS_METRICS), with = FALSE],
                by = "sampleID")
N_14B <- nrow(ips_rs)
cat(sprintf(" N (IPS RS): %d\n", N_14B))

test_b_rows <- lapply(IPS_METRICS, function(m) {
  res <- wilcox_lh(ips_rs[[m]], ips_rs$RiskGroup)
  data.table(metric = m, n_low = res$n_low, n_high = res$n_high,
             median_low = res$med_low, median_high = res$med_high,
             W = res$W, p_value = res$p)
})
test_b <- rbindlist(test_b_rows)
test_b[, p_adj_BH := p.adjust(p_value, method = "BH")]
test_b[, sig := fifelse(p_value < 0.001, "***",
                  fifelse(p_value < 0.01, "**",
                  fifelse(p_value < 0.05, "*", "")))]
setorder(test_b, p_value)
fwrite(test_b, file.path(TBLDIR, "Tab14B_IPS_metrics_test.csv"))
cat(" Tab14B_IPS_metrics_test.csv\n")
print(test_b)
n_sig_b <- sum(test_b$p_value < P_CUT, na.rm = TRUE)
cat(sprintf(" raw p<%.2f: %d / %d metrics\n", P_CUT, n_sig_b, length(IPS_METRICS)))

cat("\n[14B.4] Render Fig10C IPS metrics boxplot\n")
plot_b <- melt(ips_rs[, c("sampleID","RiskGroup", IPS_METRICS), with = FALSE],
               id.vars = c("sampleID","RiskGroup"),
               variable.name = "metric", value.name = "value")
plot_b <- plot_b[!is.na(value) & is.finite(value)]

p10C <- make_box_panel(plot_b, metric_levels = test_b$metric,
                        ncol = 3,
                        title_text = "IPS metrics by RiskGroup (Stage 14B)",
                        subtitle_text = sprintf(
                          "n = %d (Low = %d, High = %d); IOBR::IPS_calculation; Wilcoxon raw p annotated; %d / %d raw p<%.2f",
                          N_14B,
                          sum(ips_rs$RiskGroup == "Low"),
                          sum(ips_rs$RiskGroup == "High"),
                          n_sig_b, length(IPS_METRICS), P_CUT))

W10C <- 10; H10C <- 7
ggsave(file.path(FIGDIR, "Fig10C_IPS_metrics_boxplot.pdf"),
       plot = p10C, width = W10C, height = H10C, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig10C_IPS_metrics_boxplot.png"),
       plot = p10C, width = W10C, height = H10C, dpi = 600, device = ragg::agg_png)
cat(" [OK] Fig10C_IPS_metrics_boxplot.{pdf,png}\n")

###############################################################################
# 14C: oncoPredict drug sensitivity (GDSC2)
###############################################################################
cat("\n========================================================================\n")
cat(" Stage 14C: oncoPredict \n")
cat("========================================================================\n")

oncoP_ok <- requireNamespace("oncoPredict", quietly = TRUE)
if (!oncoP_ok) {
  cat("[WARN] oncoPredict not installed. Skipping 14C.\n")
  cat(" To enable: install.packages('oncoPredict')\n")
} else {
  suppressPackageStartupMessages(library(oncoPredict))
  cat(sprintf(" oncoPredict : %s\n", packageVersion("oncoPredict")))

  cat("[14C.1] Load GDSC2 training data\n")
  gdsc_expr <- readRDS(INPUT_GDSC_EXPR)
  gdsc_res <- readRDS(INPUT_GDSC_RES)
  cat(sprintf(" GDSC2 expr : %d genes x %d cell lines\n",
              nrow(gdsc_expr), ncol(gdsc_expr)))
  cat(sprintf(" GDSC2 IC50 : %d cell lines x %d drugs\n",
              nrow(gdsc_res), ncol(gdsc_res)))

  cat("[14C.2] Prepare TCGA test expression (LASSO samples 181)\n")
  expr_181_idx <- match(samples_rs, colnames(expr_mat))
  expr_181_idx <- expr_181_idx[!is.na(expr_181_idx)]
  expr_test <- expr_mat[, expr_181_idx, drop = FALSE]
  cat(sprintf(" expr_test : %d genes x %d samples\n",
              nrow(expr_test), ncol(expr_test)))

  CACHE_IC50 <- file.path(CACHE, "cache_oncoPredict_IC50_stage14.rds")
  if (file.exists(CACHE_IC50)) {
    cat(sprintf("[14C.3] Loading cached IC50 prediction: %s\n", CACHE_IC50))
    ic50_dt <- readRDS(CACHE_IC50)
  } else {
    cat("[14C.3] Running oncoPredict::calcPhenotype() in 4 batches (GDSC2)\n")
    cat(" Batch strategy reduces peak memory; per-batch CSV is saved so a\n")
    cat(" crash in one batch loses only that batch's drugs.\n")
    onco_outdir <- file.path(CACHE, "oncoPredict_output")
    if (!dir.exists(onco_outdir)) dir.create(onco_outdir, recursive = TRUE)

    n_drugs <- ncol(gdsc_res)
    BATCHES <- 4 # ~50 drugs each; keeps peak memory bounded
    batch_size <- ceiling(n_drugs / BATCHES)
    batch_idx_list <- split(seq_len(n_drugs),
                            ceiling(seq_len(n_drugs) / batch_size))

    # Move into output dir so calcPhenotype's mandatory file outputs land there
    old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)

    pred_list <- vector("list", length(batch_idx_list))
    for (b in seq_along(batch_idx_list)) {
      cols_b <- batch_idx_list[[b]]
      batch_dir <- file.path(onco_outdir, sprintf("batch%02d", b))
      if (!dir.exists(batch_dir)) dir.create(batch_dir, recursive = TRUE)
      pred_csv_b <- file.path(batch_dir, "calcPhenotype_Output",
                              "DrugPredictions.csv")
      if (file.exists(pred_csv_b)) {
        cat(sprintf(" Batch %d/%d: cached (%d drugs) skip\n",
                    b, length(batch_idx_list), length(cols_b)))
        pred_list[[b]] <- fread(pred_csv_b)
        next
      }
      cat(sprintf(" Batch %d/%d: training on drugs %d-%d (%d drugs)\n",
                  b, length(batch_idx_list),
                  min(cols_b), max(cols_b), length(cols_b)))
      setwd(batch_dir)
      tryCatch({
        calcPhenotype(
          trainingExprData = gdsc_expr,
          trainingPtype = gdsc_res[, cols_b, drop = FALSE],
          testExprData = expr_test,
          batchCorrect = "eb",
          powerTransformPhenotype = TRUE,
          removeLowVaryingGenes = 0.4, # heavier filter for memory
          minNumSamples = 10,
          printOutput = TRUE,
          removeLowVaringGenesFrom = "rawData"
        )
      }, finally = setwd(old_wd))
      gc(reset = TRUE, full = TRUE)

      if (!file.exists(pred_csv_b)) {
        stop(sprintf("Batch %d failed: %s missing", b, pred_csv_b))
      }
      pred_list[[b]] <- fread(pred_csv_b)
      cat(sprintf(" [OK] batch %d wrote %s (%d cols)\n",
                  b, pred_csv_b, ncol(pred_list[[b]])))
    }

    # Merge batches by first column (sampleID); column union via merge keeps order
    ic50_dt <- pred_list[[1]]
    setnames(ic50_dt, 1, "sampleID")
    if (length(pred_list) > 1) {
      for (b in 2:length(pred_list)) {
        nxt <- pred_list[[b]]
        setnames(nxt, 1, "sampleID")
        ic50_dt <- merge(ic50_dt, nxt, by = "sampleID", sort = FALSE)
      }
    }
    saveRDS(ic50_dt, CACHE_IC50)
    cat(sprintf(" Cached: %s (%d samples %d drugs)\n",
                CACHE_IC50, nrow(ic50_dt), ncol(ic50_dt) - 1L))
  }

  cat(sprintf(" IC50 matrix: %d samples x %d drugs\n",
              nrow(ic50_dt), ncol(ic50_dt) - 1L))
  fwrite(ic50_dt, file.path(TBLDIR, "Tab14C_oncoPredict_IC50_full.csv"))

  cat("[14C.4] Wilcoxon Low vs High per drug\n")
  drug_cols <- setdiff(colnames(ic50_dt), "sampleID")
  ic50_join <- merge(rs_cli, ic50_dt, by = "sampleID")
  rg_vec <- ic50_join$RiskGroup

  test_c_rows <- lapply(drug_cols, function(dr) {
    v <- as.numeric(ic50_join[[dr]])
    res <- wilcox_lh(v, rg_vec)
    data.table(drug = dr, n_low = res$n_low, n_high = res$n_high,
               median_low = res$med_low, median_high = res$med_high,
               W = res$W, p_value = res$p)
  })
  test_c <- rbindlist(test_c_rows)
  test_c[, p_adj_BH := p.adjust(p_value, method = "BH")]
  test_c[, log_FC_HighMinusLow := median_high - median_low]
  test_c[, direction := fifelse(log_FC_HighMinusLow < 0, "High_more_sensitive",
                          fifelse(log_FC_HighMinusLow > 0, "High_more_resistant", "tie"))]
  setorder(test_c, p_value)
  fwrite(test_c, file.path(TBLDIR, "Tab14C_oncoPredict_drug_test.csv"))
  cat(sprintf(" Tab14C_oncoPredict_drug_test.csv (%d drugs)\n", nrow(test_c)))

  n_sig_c <- sum(test_c$p_value < P_CUT, na.rm = TRUE)
  cat(sprintf(" raw p<%.2f: %d / %d drugs\n", P_CUT, n_sig_c, nrow(test_c)))

  top_c <- head(test_c, TOP_N_DRUG)
  fwrite(top_c, file.path(TBLDIR, "Tab14C_oncoPredict_top30_drugs.csv"))
  cat(" Tab14C_oncoPredict_top30_drugs.csv (top 30 by raw p)\n")
  print(head(top_c[, .(drug, n_low, n_high, median_low, median_high,
                        log_FC_HighMinusLow, p_value, p_adj_BH, direction)], 10))

  cat("\n[14C.5] Render Fig10D oncoPredict top 30 drugs boxplot\n")
  drugs_show <- top_c$drug
  plot_c <- melt(ic50_join[, c("sampleID","RiskGroup", drugs_show), with = FALSE],
                 id.vars = c("sampleID","RiskGroup"),
                 variable.name = "drug", value.name = "ic50")
  plot_c <- plot_c[!is.na(ic50) & is.finite(ic50)]

  # Annotate per-drug raw p
  pann <- top_c[, .(drug,
                    p_lbl = sprintf("p = %.2g", p_value),
                    BH = p_adj_BH)]
  plot_c[, drug := factor(drug, levels = drugs_show)]
  pann[, drug := factor(drug, levels = drugs_show)]
  ymax_dt <- plot_c[, .(y = max(ic50, na.rm = TRUE) +
                            0.10 * (max(ic50, na.rm = TRUE) -
                                     min(ic50, na.rm = TRUE) + 1e-6)),
                    by = drug]
  pann <- merge(pann, ymax_dt, by = "drug")

  p10D <- ggplot(plot_c, aes(x = RiskGroup, y = ic50, fill = RiskGroup)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, color = "black", linewidth = 0.35) +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.40, color = "grey25") +
    geom_text(data = pann, aes(x = 1.5, y = y, label = p_lbl),
              inherit.aes = FALSE, size = 2.6) +
    facet_wrap(~ drug, scales = "free_y", ncol = 6) +
    scale_fill_manual(values = RISK_COL, name = "RiskGroup") +
    labs(x = NULL, y = "Predicted log(IC50)") +
    theme_bw(base_size = 9) +
    theme(legend.position = "top",
          strip.text = element_text(face = "bold", size = 7),
          axis.text.x = element_text(face = "bold", size = 7))

  W10D <- 14
  H10D <- ceiling(TOP_N_DRUG / 6) * 1.8 + 1.5
  ggsave(file.path(FIGDIR, "Fig10D_oncoPredict_top_drugs_boxplot.pdf"),
         plot = p10D, width = W10D, height = H10D, device = cairo_pdf)
  ggsave(file.path(FIGDIR, "Fig10D_oncoPredict_top_drugs_boxplot.png"),
         plot = p10D, width = W10D, height = H10D, dpi = 600,
         device = ragg::agg_png)
  cat(" [OK] Fig10D_oncoPredict_top_drugs_boxplot.{pdf,png}\n")
}

# Summary
t_end <- Sys.time(); elapsed_total <- as.numeric(t_end - t_start, units = "mins")
cat("\n========================================================================\n")
cat(" Stage 14 Treatment Response Summary \n")
cat("========================================================================\n")
cat(sprintf(" 14A TIDE raw p<%.2f metrics : %d / %d\n", P_CUT, n_sig_a, length(TIDE_METRICS)))
cat(sprintf(" 14A Responder chi-square :\n"))
print(chisq_dt[, .(flag, method, p_value, p_adj_BH, sig)])
cat(sprintf("\n 14B IPS raw p<%.2f metrics : %d / %d\n", P_CUT, n_sig_b, length(IPS_METRICS)))
print(test_b[, .(metric, median_low, median_high, p_value, p_adj_BH, sig)])
if (oncoP_ok) {
  cat(sprintf("\n 14C oncoPredict raw p<%.2f drugs: %d / %d\n", P_CUT, n_sig_c, nrow(test_c)))
  cat(sprintf(" Top 5 by raw p:\n"))
  print(head(test_c[, .(drug, p_value, p_adj_BH, direction)], 5))
} else {
  cat("\n 14C oncoPredict: SKIPPED (package not installed)\n")
}
cat(sprintf("\n Total elapsed: %.2f min\n", elapsed_total))
cat("========================================================================\n")

sess_file <- file.path(LOGDIR, "sessionInfo/sessionInfo_stage14.txt")
sink(sess_file); cat(sprintf("Session info: %s\n\n", Sys.time())); print(sessionInfo()); sink()
cat(sprintf("\n>>> sessionInfo: %s\n", sess_file))
cat(sprintf("\n=== Stage 14 done (%s) ===\n", Sys.time()))

sink(type = "message"); sink(); close(log_con)

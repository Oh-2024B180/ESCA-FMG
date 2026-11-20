###############################################################################
## Stage 11: External Validation GSE53625 (n=179 ESCC) + AnnoProbe
##
## Model: 7-gene protein-coding LASSO-Cox (Stage 10)
## FOLR3, ZBTB7C, RPRML, HEPACAM2, CLCA1, CACNG8, LIPC
##
## 11A: GSE53625 ESCC (GPL18109); GSE180791 excluded
## 11B: gene-level z-score normalization (cross-platform)
## 11C: KM (log-rank), time-dependent ROC, Concordance index
##
## Probe annotation:
## - probe SYMBOL from **AnnoProbe::idmap("GPL18109")** (not biomaRt ENST)
## AnnoProbe maps probes by sequence alignment and covers all 71,584 probes;
## the biomaRt ENST mapping covered 4.3% (3,054 probes) and none of the 7 genes.
## - Clinical extraction is column-name-aware
##
## Settings:
## 1. probe annotation : AnnoProbe::idmap
## 2. GSE180791 : excluded
## 3. Survival units : months -> days (x 30.4375)
## 4. Sample filter : tumor only (179 ESCC)
## 5. Median cutoff : (a) TCGA training first, (b) fallback if imbalanced
##
## Random seed: 1262118388
###############################################################################

# run from the repository root (see README.md)

INPUT_GSE <- "data/GSE53625_series_matrix.txt.gz"
INPUT_LASSO <- "tables/Tab10_LASSO_selected_genes.csv"
INPUT_RS_TCGA <- "tables/Tab10_risk_scores_per_sample.csv"
INPUT_TCGA_KM <- "tables/Tab10_KM_logrank_results.csv"
INPUT_TCGA_RC <- "tables/Tab10_timeROC_AUC.csv"

TBLDIR <- "tables"; FIGDIR <- "figures"; LOGDIR <- "logs"
CACHE <- "cache"
for (d in c(TBLDIR, FIGDIR, LOGDIR, CACHE)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

get_tcga_median <- function() median(fread(INPUT_RS_TCGA)$RiskScore)

MONTHS_TO_DAYS <- 365.25 / 12
ROC_TIMES_DAYS <- c(365, 1095, 1825)
ROC_FALLBACK <- c(365, 730, 1095)
AT_RISK_MIN <- 10
MIN_GROUP_FRAC <- 0.10
RISK_COL <- c(Low = "#377EB8", High = "#E41A1C")

set.seed(1262118388)
t_start <- Sys.time()

log_file <- file.path(LOGDIR, "11_diagnostics.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE); sink(log_con, type = "message")
cat(sprintf("=== Stage 11 (AnnoProbe path) start (%s) ===\n\n", Sys.time()))

cat("[1/14] Loading packages\n")
suppressPackageStartupMessages({
  library(data.table); library(GEOquery); library(Biobase)
  library(AnnoProbe)
  library(survival); library(timeROC); library(survminer)
  library(ggplot2); library(ragg); library(cowplot)
})
cat(sprintf(" AnnoProbe : %s\n", packageVersion("AnnoProbe")))

cat("\n[2/14] Loading Stage 10 model\n")
lasso_dt <- fread(INPUT_LASSO)
rs_tcga_dt <- fread(INPUT_RS_TCGA)
km_tcga_dt <- fread(INPUT_TCGA_KM)
TCGA_RS_MEDIAN <- get_tcga_median()
cat(sprintf(" LASSO genes: %d (%s)\n", nrow(lasso_dt),
            paste(lasso_dt$gene, collapse = ", ")))
cat(sprintf(" TCGA RS median: %.4f\n", TCGA_RS_MEDIAN))

cat("\n[3/14] GEOquery: loading GSE53625\n")
eset_list <- tryCatch(getGEO(filename = INPUT_GSE, AnnotGPL = FALSE,
                              destdir = CACHE, GSEMatrix = TRUE, getGPL = FALSE),
                       error = function(e) { cat(sprintf(" [warn] %s\n", conditionMessage(e))); NULL })
if (is.null(eset_list)) stop("GEOquery failed.")
eset <- if (inherits(eset_list, "list")) eset_list[[1]] else eset_list
expr_raw <- exprs(eset); pheno <- pData(eset)
cat(sprintf(" ExpressionSet: %d probes x %d samples\n", nrow(expr_raw), ncol(expr_raw)))
cat(sprintf(" expr value range: [%.2f, %.2f] (median=%.2f) Agilent log2 scale\n",
            min(expr_raw, na.rm = TRUE), max(expr_raw, na.rm = TRUE),
            median(expr_raw, na.rm = TRUE)))

cat("\n[3.5/14] AnnoProbe::idmap('GPL18109') for probe SYMBOL\n")
ids <- AnnoProbe::idmap("GPL18109", type = "pipe", destdir = CACHE)
cat(sprintf(" AnnoProbe annotation: %d rows x %d cols\n", nrow(ids), ncol(ids)))
cat(sprintf(" columns: %s\n", paste(colnames(ids), collapse = ", ")))
sym_col <- grep("symbol|SYMBOL|gene|GENE", colnames(ids),
                value = TRUE, ignore.case = TRUE)[1]
probe_col <- colnames(ids)[1]
cat(sprintf(" probe column: '%s' symbol column: '%s'\n", probe_col, sym_col))
cat(sprintf(" unique probes : %d\n", length(unique(ids[[probe_col]]))))
cat(sprintf(" unique SYMBOLs: %d\n", length(unique(ids[[sym_col]]))))

ids_clean <- ids[!is.na(ids[[sym_col]]) & ids[[sym_col]] != "", ]
cat(sprintf(" after dropping empty SYMBOL: %d rows\n", nrow(ids_clean)))

cat("\n[4/14] Filter tumor samples\n")
# direct column access (Li 2014 / GEOquery 2.74+ format)
tissue_col <- grep("tissue", colnames(pheno), value = TRUE, ignore.case = TRUE)
cat(sprintf(" tissue columns: %s\n", paste(tissue_col, collapse = ", ")))
tissue_col1 <- tissue_col[grepl(":ch1$", tissue_col)][1]
if (is.na(tissue_col1)) tissue_col1 <- tissue_col[1]
cat(sprintf(" using: %s\n", tissue_col1))
tumor_idx <- grep("cancer tissue", pheno[[tissue_col1]], ignore.case = TRUE)
tumor_samples <- rownames(pheno)[tumor_idx]
cat(sprintf(" tumor samples: %d / %d\n", length(tumor_samples), nrow(pheno)))

cat("\n[5/14] Extract clinical OS event + OS time\n")
death_col <- grep("death at fu", colnames(pheno), value = TRUE, ignore.case = TRUE)[1]
surv_col <- grep("survival time", colnames(pheno), value = TRUE, ignore.case = TRUE)[1]
cat(sprintf(" death column: %s\n", death_col))
cat(sprintf(" surv column: %s\n", surv_col))

pdat_t <- pheno[tumor_samples, ]
death_raw <- as.character(pdat_t[[death_col]])
surv_raw <- as.character(pdat_t[[surv_col]])
os_event <- ifelse(tolower(trimws(death_raw)) == "yes", 1L,
              ifelse(tolower(trimws(death_raw)) == "no", 0L, NA_integer_))
os_months <- suppressWarnings(as.numeric(surv_raw))
os_days <- os_months * MONTHS_TO_DAYS

clin_dt <- data.table(sampleID = tumor_samples, OS = os_event,
                      OS_months = os_months, OS.time = os_days)
clin_dt <- clin_dt[!is.na(OS) & !is.na(OS.time) & OS.time > 0]
N_val <- nrow(clin_dt); samples_keep <- clin_dt$sampleID
cat(sprintf(" N_validation: %d (events=%d, %.1f%%)\n",
            N_val, sum(clin_dt$OS), 100 * mean(clin_dt$OS)))
cat(sprintf(" OS.time (days): median=%.0f, max=%.0f\n",
            median(clin_dt$OS.time), max(clin_dt$OS.time)))

cat("\n[6/14] 7 LASSO gene coverage check (AnnoProbe mapping)\n")
lasso_genes <- lasso_dt$gene
cov_dt <- data.table(gene = lasso_genes, coefficient = lasso_dt$coefficient)
cov_dt[, n_probes := vapply(gene, function(g)
       sum(ids_clean[[sym_col]] == g, na.rm = TRUE), integer(1))]
cov_dt[, in_GSE53625 := n_probes > 0]
cov_dt[, status := fifelse(in_GSE53625, "MAPPED", "MISSING")]
fwrite(cov_dt, file.path(TBLDIR, "Tab11_GSE53625_LASSO_gene_coverage.csv"))
print(cov_dt)
n_mapped <- sum(cov_dt$in_GSE53625)
cat(sprintf("\n Coverage: %d / %d (%.1f%%)\n",
            n_mapped, length(lasso_genes), 100 * n_mapped / length(lasso_genes)))
if (n_mapped == 0) stop("0 LASSO genes mapped cannot compute Risk Score.")

cat("\n[7/14] Probe Gene aggregation (mean rule, LASSO genes only)\n")
mapped_genes <- cov_dt[in_GSE53625 == TRUE, gene]
expr_lasso <- matrix(NA_real_, nrow = length(mapped_genes),
                     ncol = length(samples_keep),
                     dimnames = list(mapped_genes, samples_keep))
for (g in mapped_genes) {
  probe_ids_g <- as.character(ids_clean[[probe_col]][ids_clean[[sym_col]] == g])
  probe_ids_g <- intersect(probe_ids_g, rownames(expr_raw))
  if (length(probe_ids_g) == 1) {
    expr_lasso[g, ] <- expr_raw[probe_ids_g, samples_keep]
  } else if (length(probe_ids_g) > 1) {
    expr_lasso[g, ] <- colMeans(expr_raw[probe_ids_g, samples_keep, drop = FALSE],
                                 na.rm = TRUE)
  }
}
cat(sprintf(" LASSO gene expression matrix: %d %d\n",
            nrow(expr_lasso), ncol(expr_lasso)))
cat(" per-gene NA count:\n"); print(rowSums(is.na(expr_lasso)))

cat("\n[8/14] Z-score normalization per gene\n")
expr_lasso_z <- t(scale(t(expr_lasso), center = TRUE, scale = TRUE))
cat(sprintf(" z-scored matrix: %d %d\n", nrow(expr_lasso_z), ncol(expr_lasso_z)))

cat("\n[9/14] Risk Score computation\n")
coef_lookup <- setNames(lasso_dt$coefficient, lasso_dt$gene)
matched_coefs <- coef_lookup[mapped_genes]
risk_score <- as.numeric(t(expr_lasso_z) %*% matched_coefs)
names(risk_score) <- samples_keep
cat(" RS dist:\n"); print(summary(risk_score))

cat("\n[10/14] Median cutoff strategy\n")
group_a <- ifelse(risk_score <= TCGA_RS_MEDIAN, "Low", "High")
frac_a <- min(sum(group_a == "Low"), sum(group_a == "High")) / length(group_a)
cat(sprintf(" (a) TCGA median %.4f: Low=%d High=%d min_frac=%.1f%%\n",
            TCGA_RS_MEDIAN, sum(group_a == "Low"), sum(group_a == "High"), 100 * frac_a))
if (frac_a < MIN_GROUP_FRAC) {
  cat(sprintf(" [FALLBACK] min_frac < 10%% validation median\n"))
  cutoff_used <- median(risk_score); cutoff_source <- "validation_median"
} else {
  cat(" [USE TCGA] balance OK\n")
  cutoff_used <- TCGA_RS_MEDIAN; cutoff_source <- "TCGA_training_median"
}
risk_group <- factor(ifelse(risk_score <= cutoff_used, "Low", "High"),
                     levels = c("Low","High"))
cat(sprintf(" cutoff: %.4f (%s) Low/High: %d/%d\n", cutoff_used, cutoff_source,
            sum(risk_group == "Low"), sum(risk_group == "High")))

rs_dt <- data.table(sampleID = samples_keep, OS = clin_dt$OS, OS.time = clin_dt$OS.time,
                    OS_months = clin_dt$OS_months, RiskScore = risk_score,
                    RiskGroup = as.character(risk_group),
                    cutoff_used = cutoff_used, cutoff_src = cutoff_source)
fwrite(rs_dt, file.path(TBLDIR, "Tab11_GSE53625_RS_per_sample.csv"))

cat("\n[11/14] KM + log-rank + HR\n")
km_data <- data.frame(time = clin_dt$OS.time, event = clin_dt$OS, RiskGroup = risk_group)
fit_km <- survfit(Surv(time, event) ~ RiskGroup, data = km_data)
sd_lr <- survdiff(Surv(time, event) ~ RiskGroup, data = km_data)
lr_p <- 1 - pchisq(sd_lr$chisq, length(sd_lr$n) - 1)
fit_cox <- coxph(Surv(time, event) ~ RiskGroup, data = km_data)
hrs <- summary(fit_cox)
HR <- hrs$conf.int[1, "exp(coef)"]
HR_low <- hrs$conf.int[1, "lower .95"]; HR_high <- hrs$conf.int[1, "upper .95"]
cox_p <- hrs$coefficients[1, "Pr(>|z|)"]
med_os <- summary(fit_km)$table
cat(sprintf(" log-rank p: %.4g\n", lr_p))
cat(sprintf(" HR (High vs Low): %.3f (95%% CI %.3f %.3f, Cox p=%.4g)\n",
            HR, HR_low, HR_high, cox_p))
fwrite(data.table(cohort = "GSE53625", n = N_val,
                  n_low = sum(risk_group == "Low"), n_high = sum(risk_group == "High"),
                  median_OS_low_days = med_os["RiskGroup=Low", "median"],
                  median_OS_high_days = med_os["RiskGroup=High", "median"],
                  HR = HR, HR_low_95 = HR_low, HR_high_95 = HR_high,
                  cox_p_value = cox_p, logrank_p = lr_p,
                  cutoff_used = cutoff_used, cutoff_source = cutoff_source),
       file.path(TBLDIR, "Tab11_GSE53625_KM_logrank.csv"))

cat("\n[12/14] timeROC + Concordance index\n")
at_risk_5yr <- sum(clin_dt$OS.time >= 1825)
cat(sprintf(" at-risk @ 5-yr: %d\n", at_risk_5yr))
if (at_risk_5yr < AT_RISK_MIN) {
  roc_times_used <- ROC_FALLBACK; roc_labels <- c("1-yr","2-yr","3-yr")
  cat(" [WARN] fallback to 1/2/3-yr\n")
} else {
  roc_times_used <- ROC_TIMES_DAYS; roc_labels <- c("1-yr","3-yr","5-yr")
}
roc_res <- timeROC(T = clin_dt$OS.time, delta = clin_dt$OS, marker = risk_score,
                   cause = 1, times = roc_times_used, iid = TRUE)
auc_vec <- as.numeric(roc_res$AUC)
# 2026-10-06: vect_sd_1 already is the standard error of the AUC; the earlier sqrt() made the 95% CI
# far too wide. Corrected table: helpers/11_fix_cindex_timeROC_tables.R
sd_vec <- tryCatch(as.numeric(roc_res$inference$vect_sd_1),
                      error = function(e) rep(NA_real_, length(auc_vec)))
auc_low <- pmax(0, auc_vec - 1.96 * sd_vec)
auc_high <- pmin(1, auc_vec + 1.96 * sd_vec)
roc_tbl <- data.table(cohort = "GSE53625", time_label = roc_labels,
                      time_days = roc_times_used, AUC = auc_vec,
                      AUC_low_95 = auc_low, AUC_high_95 = auc_high, SE = sd_vec)
fwrite(roc_tbl, file.path(TBLDIR, "Tab11_GSE53625_timeROC_AUC.csv"))
print(roc_tbl)

# 2026-10-06: the C-index of the FIXED score, higher score = shorter survival. The earlier code read
# summary(coxph(Surv ~ risk_score))$concordance from a model refitted inside GSE53625; a refit
# re-estimates the sign of the score, so the inverted association seen here (HR 0.85, AUC below 0.5)
# still gave a value above 0.5 (0.517). The fixed score has C = 0.483 (95% CI 0.424-0.543).
cidx_fixed <- function(time, event, score) {
  k <- concordance(Surv(time, event) ~ score, reverse = TRUE)
  c(C = unname(k$concordance), se = sqrt(k$var))
}
cidx_v <- cidx_fixed(clin_dt$OS.time, clin_dt$OS, risk_score)
cidx_t <- cidx_fixed(rs_tcga_dt$OS.time, rs_tcga_dt$OS, rs_tcga_dt$RiskScore)
c_v <- cidx_v["C"]; c_v_se <- cidx_v["se"]
c_t <- cidx_t["C"]; c_t_se <- cidx_t["se"]
cidx_out <- data.table(cohort = c("TCGA-ESCA training","GSE53625 validation"),
                       n = c(nrow(rs_tcga_dt), N_val),
                       c_index = c(c_t, c_v),
                       c_low_95 = c(c_t - 1.96 * c_t_se, c_v - 1.96 * c_v_se),
                       c_high_95 = c(c_t + 1.96 * c_t_se, c_v + 1.96 * c_v_se),
                       c_se = c(c_t_se, c_v_se))
fwrite(cidx_out, file.path(TBLDIR, "Tab11_GSE53625_concordance_index.csv"))
cat(sprintf(" C-index: TCGA=%.3f GSE53625=%.3f\n", c_t, c_v))

cat("\n[13/14] Generating figures Fig7A-7E\n")

# Fig7A KM
# HR/CI/p annotation moved from title to upper-right INSIDE the plot panel.
ss <- ggsurvplot(fit_km, data = km_data, pval = FALSE, risk.table = TRUE, conf.int = FALSE,
                 palette = c(unname(RISK_COL["Low"]), unname(RISK_COL["High"])),
                 legend.labs = c("Low risk","High risk"), legend.title = "Risk group",
                 xlab = "Time (days)", ylab = "Overall survival",
                 title = NULL,
                 ggtheme = theme_bw(base_size = 11), risk.table.fontsize = 3.5)
hr_lbl_7A <- sprintf("HR = %.2f (95%% CI %.2f - %.2f)\nlog-rank p = %.3g",
                     HR, HR_low, HR_high, lr_p)
ss$plot <- ss$plot +
  annotate("label", x = Inf, y = Inf, label = hr_lbl_7A,
           hjust = 1.05, vjust = 1.2, size = 3.6, fontface = "plain",
           label.size = 0.3, fill = "white", label.r = unit(0, "pt"))
combined_7A <- cowplot::plot_grid(ss$plot + theme(legend.position = "top"),
                                   ss$table + theme(plot.title = element_blank()),
                                   ncol = 1, align = "v", axis = "lr",
                                   rel_heights = c(3, 1))
ggsave(file.path(FIGDIR, "Fig7A_GSE53625_KM_by_risk_group.pdf"),
       plot = combined_7A, width = 8, height = 7, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig7A_GSE53625_KM_by_risk_group.png"),
       plot = combined_7A, width = 8, height = 7, dpi = 600, device = ragg::agg_png)
cat(" [OK] Fig7A_GSE53625_KM_by_risk_group.{pdf,png}\n")

# Fig7B timeROC
roc_dfs <- lapply(seq_along(roc_times_used), function(k) {
  data.frame(FPR = roc_res$FP[, k], TPR = roc_res$TP[, k],
             time_label = sprintf("%s (AUC=%.3f, 95%% CI %.3f-%.3f)",
                                  roc_labels[k], auc_vec[k], auc_low[k], auc_high[k]))
})
roc_df <- do.call(rbind, roc_dfs)
palette_used <- c("#377EB8","#4DAF4A","#E41A1C")
names(palette_used) <- unique(roc_df$time_label)
p7B <- ggplot(roc_df, aes(x = FPR, y = TPR, color = time_label)) +
  geom_line(linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = palette_used, name = NULL) + coord_equal() +
  labs(x = "False positive rate", y = "True positive rate") +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.65, 0.20),
        legend.background = element_rect(fill = "white", color = "grey70"))
ggsave(file.path(FIGDIR, "Fig7B_GSE53625_timeROC.pdf"),
       plot = p7B, width = 6.5, height = 6, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig7B_GSE53625_timeROC.png"),
       plot = p7B, width = 6.5, height = 6, dpi = 600, device = ragg::agg_png)
cat(" [OK] Fig7B_GSE53625_timeROC.{pdf,png}\n")

# Fig7C 3-panel
ord <- order(risk_score)
rs_sorted <- risk_score[ord]; samples_sort <- samples_keep[ord]
rg_sorted <- as.character(risk_group)[ord]
os_event_s <- clin_dt$OS[ord]; os_time_s <- clin_dt$OS.time[ord]
N <- length(ord)
p1 <- ggplot(data.frame(rank = 1:N, RS = rs_sorted,
                         RiskGroup = factor(rg_sorted, levels = c("Low","High"))),
             aes(x = rank, y = RS, color = RiskGroup)) +
  geom_point(size = 0.8) +
  geom_hline(yintercept = cutoff_used, linetype = "dashed", color = "grey40") +
  annotate("text", x = N * 0.05, y = cutoff_used,
           label = sprintf("cutoff = %.4f (%s)", cutoff_used,
                           ifelse(cutoff_source == "TCGA_training_median",
                                  "TCGA training", "validation median")),
           vjust = -0.5, size = 3, color = "grey40") +
  scale_color_manual(values = c(Low = unname(RISK_COL["Low"]), High = unname(RISK_COL["High"])),
                     name = "Risk group") +
  labs(x = NULL, y = "Risk Score") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right", axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), panel.grid.minor = element_blank())
p2 <- ggplot(data.frame(rank = 1:N, OS_time = os_time_s,
                         Event = factor(ifelse(os_event_s == 1, "Death","Censored"),
                                         levels = c("Censored","Death"))),
             aes(x = rank, y = OS_time, color = Event)) +
  geom_point(size = 0.9, alpha = 0.8) +
  scale_color_manual(values = c(Censored = "grey30", Death = "#E41A1C"), name = "OS event") +
  labs(x = NULL, y = "Time (days)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right", axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), panel.grid.minor = element_blank())

heat_mat <- expr_lasso_z[, ord, drop = FALSE]
heat_dt2 <- as.data.table(heat_mat, keep.rownames = "gene")
heat_long <- melt(heat_dt2, id.vars = "gene", variable.name = "sample", value.name = "z")
heat_long[, sample := factor(sample, levels = samples_sort)]
heat_long[, gene := factor(gene, levels = rev(mapped_genes))]
p3 <- ggplot(heat_long, aes(x = sample, y = gene, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#0571B0", mid = "white", high = "#CA0020",
                        midpoint = 0, name = "z-score",
                        limits = c(-3, 3), oob = scales::squish) +
  labs(x = sprintf("GSE53625 samples sorted by Risk Score (low high, n=%d)", N), y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8), panel.grid = element_blank(),
        legend.position = "right")
combined_7C <- cowplot::plot_grid(p1, p2, p3, ncol = 1, align = "v", axis = "lr",
                                   rel_heights = c(1.1, 1.1,
                                                    max(1.5, 0.18 * length(mapped_genes) + 0.6)))
H_7C <- 3 + 0.18 * length(mapped_genes) + 1
ggsave(file.path(FIGDIR, "Fig7C_GSE53625_RS_distribution_3panel.pdf"),
       plot = combined_7C, width = 11, height = H_7C, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(FIGDIR, "Fig7C_GSE53625_RS_distribution_3panel.png"),
       plot = combined_7C, width = 11, height = H_7C, dpi = 600,
       device = ragg::agg_png, limitsize = FALSE)
cat(" [OK] Fig7C_GSE53625_RS_distribution_3panel.{pdf,png}\n")

# Fig7D HR forest
forest_dt <- data.table(
  cohort = c("TCGA-ESCA training (n=181)", sprintf("GSE53625 validation (n=%d)", N_val)),
  HR = c(km_tcga_dt$HR, HR),
  HR_low = c(km_tcga_dt$HR_low_95, HR_low),
  HR_high = c(km_tcga_dt$HR_high_95, HR_high),
  p_value = c(km_tcga_dt$logrank_p, lr_p))
forest_dt[, cohort := factor(cohort, levels = rev(cohort))]
forest_dt[, label := sprintf("HR=%.2f (95%% CI %.2f-%.2f), p=%.2g",
                              HR, HR_low, HR_high, p_value)]
p7D <- ggplot(forest_dt, aes(y = cohort, x = HR)) +
  geom_errorbarh(aes(xmin = HR_low, xmax = HR_high), height = 0.2, color = "black") +
  geom_point(size = 4, color = "#E41A1C") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = label), hjust = 0, nudge_x = 0.05, size = 3.5) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10),
                limits = c(0.4, max(forest_dt$HR_high) * 3.5)) +
  labs(x = "Hazard ratio (High vs Low risk, log scale)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "bold"))
ggsave(file.path(FIGDIR, "Fig7D_TCGA_vs_GSE53625_HR_comparison.pdf"),
       plot = p7D, width = 10, height = 4, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig7D_TCGA_vs_GSE53625_HR_comparison.png"),
       plot = p7D, width = 10, height = 4, dpi = 600, device = ragg::agg_png)
cat(" [OK] Fig7D_TCGA_vs_GSE53625_HR_comparison.{pdf,png}\n")

# Fig7E c-index bar
p7E <- ggplot(cidx_out, aes(x = cohort, y = c_index, fill = cohort)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = c_low_95, ymax = c_high_95),
                width = 0.2, color = "black", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", c_index)),
            vjust = -1.5, size = 4, fontface = "bold") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.10))) +
  scale_fill_manual(values = c("TCGA-ESCA training" = "#377EB8",
                                "GSE53625 validation" = "#E41A1C"), guide = "none") +
  labs(x = NULL, y = "Concordance index") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(face = "bold"))
ggsave(file.path(FIGDIR, "Fig7E_TCGA_vs_GSE53625_cindex_comparison.pdf"),
       plot = p7E, width = 7, height = 5, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig7E_TCGA_vs_GSE53625_cindex_comparison.png"),
       plot = p7E, width = 7, height = 5, dpi = 600, device = ragg::agg_png)
cat(" [OK] Fig7E_TCGA_vs_GSE53625_cindex_comparison.{pdf,png}\n")

# Summary
t_end <- Sys.time(); elapsed <- as.numeric(t_end - t_start, units = "mins")
cat("\n========================================================================\n")
cat(" Stage 11 (AnnoProbe) External Validation Summary\n")
cat("========================================================================\n")
cat(sprintf(" Cohort : GSE53625 ESCC (Li 2014)\n"))
cat(sprintf(" Annotation : AnnoProbe::idmap('GPL18109') (sequence-based)\n"))
cat(sprintf(" N (valid OS) : %d\n", N_val))
cat(sprintf(" Events : %d (%.1f%%)\n", sum(clin_dt$OS), 100 * mean(clin_dt$OS)))
cat(sprintf(" 7 LASSO mapped: %d / 7\n", n_mapped))
cat(sprintf(" cutoff source : %s (%.4f)\n", cutoff_source, cutoff_used))
cat(sprintf(" Low / High n : %d / %d\n",
            sum(risk_group == "Low"), sum(risk_group == "High")))
cat(sprintf(" log-rank p : %.4g\n", lr_p))
cat(sprintf(" HR : %.3f (95%% CI %.3f %.3f)\n", HR, HR_low, HR_high))
cat(sprintf(" C-idx valid. : %.3f (TCGA: %.3f)\n", c_v, c_t))
cat("\n AUC table:\n"); print(roc_tbl)
cat(sprintf("\n Total elapsed: %.2f min\n", elapsed))
cat("========================================================================\n")

sess_file <- file.path(LOGDIR, "sessionInfo/sessionInfo_stage11.txt")
sink(sess_file); cat(sprintf("Session info: %s\n\n", Sys.time())); print(sessionInfo()); sink()
cat(sprintf("\n>>> sessionInfo: %s\n", sess_file))
cat(sprintf("\n=== Stage 11 (AnnoProbe) done (%s) ===\n", Sys.time()))
sink(type = "message"); sink(); close(log_con)

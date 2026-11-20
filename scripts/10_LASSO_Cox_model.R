###############################################################################
## Stage 10: LASSO-Cox Prognostic Risk Model (protein-coding only)
##
## - 10A: Univariate Cox p < 0.01
## - 10B: model-friendly gene SYMBOL (alphanumeric, no '-' start, no special)
## - 10C: glmnet family="cox", alpha=1, 10-fold CV, lambda.min
## - 10D: RS = sum(coefficient x expression)
## - 10E: Median RS cutoff
## - 10F: KM + timeROC + RS distribution
##
## The stage 9 candidate pool is protein-coding only,
## so the genes satisfy 10B (lncRNAs such as RP11-* are excluded).
##
## Fig6A/6B: to avoid title-axis overlap,
## par(mar = c(5, 4, 6, 4)) and mtext title at line = 4.
##
## Inputs:
## data/expr_ESCA_tumor.csv
## data/survival_ESCA_survival.csv
## data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv
## tables/Tab9_DEG_significant_C1_vs_C2.csv (protein-coding only)
##
## Outputs:
## tables/Tab10_*.csv (univariate Cox, LASSO selected, RS, KM, timeROC)
## figures/Fig6A-Fig6E.{pdf,png}
## logs/10_lasso_cox_diagnostics.txt + sessionInfo_stage10.txt
## cache/cache_LASSO_cv_glmnet_stage10.rds + cache_LASSO_glmnet_stage10.rds
##
## Random seed: 1262118388 (CLAUDE.md §4.3)
###############################################################################

# run from the repository root (see README.md)

INPUT_EXPR <- "data/expr_ESCA_tumor.csv"
INPUT_SURV <- "data/survival_ESCA_survival.csv"
INPUT_CLI <- "data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv"
INPUT_DEG <- "tables/Tab9_DEG_significant_C1_vs_C2.csv"

TBLDIR <- "tables"
FIGDIR <- "figures"
LOGDIR <- "logs"
CACHE <- "cache"
for (d in c(TBLDIR, FIGDIR, LOGDIR, CACHE)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

COX_P_CUT <- 0.01
LASSO_NFOLDS <- 10
ROC_TIMES <- c(365, 1095, 1825)
ROC_FALLBACK <- c(365, 730, 1095)
AT_RISK_MIN <- 10

RISK_COL <- c(Low = "#377EB8", High = "#E41A1C")

set.seed(1262118388)
t_start <- Sys.time()

log_file <- file.path(LOGDIR, "10_diagnostics.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE); sink(log_con, type = "message")
cat(sprintf("=== Stage 10 (PC-only candidate pool) start (%s) ===\n\n", Sys.time()))

cat("[1/14] Loading packages\n")
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(glmnet)
  library(timeROC)
  library(survminer)
  library(ggplot2)
  library(ragg)
  library(cowplot)
})

cat("\n[2/14] Loading inputs and 3-way sample intersect\n")
expr_dt <- fread(INPUT_EXPR)
surv_dt <- fread(INPUT_SURV)
cli_dt <- fread(INPUT_CLI)
deg_dt <- fread(INPUT_DEG)

cat(sprintf(" expr cols (samples) : %d\n", ncol(expr_dt) - 1))
cat(sprintf(" survival rows : %d\n", nrow(surv_dt)))
cat(sprintf(" candidate DEG (Stage 9 PC-only) : %d\n", nrow(deg_dt)))

key_g <- colnames(expr_dt)[1]; setnames(expr_dt, key_g, "gene_id")
expr_samples <- setdiff(colnames(expr_dt), "gene_id")
primary_tumor_samples <- cli_dt[sample_type == "Primary Tumor", sampleID]
surv_sample_col <- intersect(c("sample","sampleID","Sample"), colnames(surv_dt))[1]
setnames(surv_dt, surv_sample_col, "sampleID")

samples_3way <- Reduce(intersect, list(expr_samples, primary_tumor_samples, surv_dt$sampleID))
cat(sprintf(" 3-way intersect : %d\n", length(samples_3way)))

surv_keep <- surv_dt[sampleID %in% samples_3way][order(match(sampleID, samples_3way))]
surv_keep <- surv_keep[!is.na(OS) & !is.na(OS.time) & OS.time > 0]
samples_keep <- surv_keep$sampleID
N_train <- length(samples_keep)
cat(sprintf(" After OS/OS.time NA filter : %d (N_train)\n", N_train))

cat("\n[3/14] Candidate gene matrix\n")
candidate_genes <- deg_dt$gene
n_in_expr <- sum(candidate_genes %in% expr_dt$gene_id)
cat(sprintf(" candidate genes detected in expr: %d / %d (%.1f%%)\n",
            n_in_expr, length(candidate_genes),
            100 * n_in_expr / length(candidate_genes)))

cand_dt <- expr_dt[gene_id %in% candidate_genes, c("gene_id", samples_keep), with = FALSE]
cand_mat <- as.matrix(cand_dt[, ..samples_keep])
rownames(cand_mat) <- cand_dt$gene_id
cat(sprintf(" candidate matrix dim : %d genes %d samples\n\n",
            nrow(cand_mat), ncol(cand_mat)))

os_event <- as.integer(surv_keep$OS)
os_time <- as.numeric(surv_keep$OS.time)

cat("[4/14] Univariate Cox screening\n")
t_uni <- Sys.time()
uni_rows <- vector("list", nrow(cand_mat))
for (i in seq_len(nrow(cand_mat))) {
  gene <- rownames(cand_mat)[i]
  x_g <- as.numeric(cand_mat[i, ])
  res <- tryCatch({
    df <- data.frame(time = os_time, event = os_event, x = x_g)
    fit <- coxph(Surv(time, event) ~ x, data = df)
    s <- summary(fit)
    data.table(gene = gene, n = s$n, n_event = s$nevent,
               coef = s$coefficients[1, "coef"], HR = s$coefficients[1, "exp(coef)"],
               HR_low = s$conf.int[1, "lower .95"], HR_high = s$conf.int[1, "upper .95"],
               z = s$coefficients[1, "z"], p_value = s$coefficients[1, "Pr(>|z|)"],
               logrank_p = s$sctest["pvalue"])
  }, error = function(e) data.table(gene = gene, n = NA, n_event = NA,
               coef = NA_real_, HR = NA_real_, HR_low = NA_real_, HR_high = NA_real_,
               z = NA_real_, p_value = NA_real_, logrank_p = NA_real_))
  uni_rows[[i]] <- res
}
uni_dt <- rbindlist(uni_rows)
elapsed_uni <- as.numeric(Sys.time() - t_uni, units = "mins")
cat(sprintf(" univariate Cox done in %.2f min\n", elapsed_uni))
uni_dt[, p_adj_BH := p.adjust(p_value, method = "BH")]
fwrite(uni_dt[order(p_value)], file.path(TBLDIR, "Tab10_univariate_Cox_results.csv"))
uni_p001 <- uni_dt[!is.na(p_value) & p_value < COX_P_CUT][order(p_value)]
fwrite(uni_p001, file.path(TBLDIR, "Tab10_univariate_Cox_p001_genes.csv"))
cat(sprintf(" Tab10_univariate_Cox_results.csv (%d genes)\n", nrow(uni_dt)))
cat(sprintf(" Tab10_univariate_Cox_p001_genes.csv (%d pass p<%.2f)\n\n",
            nrow(uni_p001), COX_P_CUT))

cat("[5/14] LASSO-Cox\n")
screened_genes <- uni_p001$gene
x_lasso <- t(cand_mat[screened_genes, , drop = FALSE])
y_lasso <- Surv(time = os_time, event = os_event)
cat(sprintf(" LASSO design: %d samples %d genes\n", nrow(x_lasso), ncol(x_lasso)))

set.seed(1262118388)
cv_fit <- cv.glmnet(x_lasso, y_lasso, family = "cox", alpha = 1, nfolds = LASSO_NFOLDS,
                    standardize = TRUE, type.measure = "deviance")
cat(sprintf(" lambda.min : %.6f\n", cv_fit$lambda.min))
cat(sprintf(" lambda.1se : %.6f\n", cv_fit$lambda.1se))

# Save cache for cosmetic re-renders
saveRDS(cv_fit, file.path(CACHE, "cache_LASSO_cv_glmnet_stage10.rds"))
saveRDS(cv_fit$glmnet.fit, file.path(CACHE, "cache_LASSO_glmnet_stage10.rds"))

coefs_full <- coef(cv_fit, s = "lambda.min")
coef_vec <- as.numeric(coefs_full); names(coef_vec) <- rownames(coefs_full)
sel_idx <- which(coef_vec != 0)
selected_genes <- names(coef_vec)[sel_idx]
selected_coefs <- coef_vec[sel_idx]
n_sel <- length(selected_genes)
cat(sprintf(" LASSO selected : %d genes (Liu 2024 = 7 for AML)\n", n_sel))

sel_dt <- data.table(gene = selected_genes, coefficient = selected_coefs,
                     HR = exp(selected_coefs),
                     direction = ifelse(selected_coefs > 0, "risk_up", "protective"))
sel_dt <- sel_dt[order(-abs(coefficient))]
fwrite(sel_dt, file.path(TBLDIR, "Tab10_LASSO_selected_genes.csv"))

formula_str <- paste(sprintf("(%+.4f %s)", selected_coefs, selected_genes), collapse = " ")
cat(sprintf(" Risk Score formula:\n RS = %s\n\n", formula_str))

cat("[6/14] Risk Score + median cutoff\n")
x_sel <- t(cand_mat[selected_genes, , drop = FALSE])
risk_score <- as.numeric(x_sel %*% selected_coefs)
names(risk_score) <- samples_keep
med_cutoff <- median(risk_score)
risk_group <- factor(ifelse(risk_score <= med_cutoff, "Low", "High"), levels = c("Low","High"))
rs_dt <- data.table(sampleID = samples_keep, OS = os_event, OS.time = os_time,
                    RiskScore = risk_score, RiskGroup = as.character(risk_group))
fwrite(rs_dt, file.path(TBLDIR, "Tab10_risk_scores_per_sample.csv"))
cat(sprintf(" median RS cutoff : %.4f\n", med_cutoff))
cat(sprintf(" Low / High n : %d / %d\n\n",
            sum(risk_group == "Low"), sum(risk_group == "High")))

cat("[7/14] KM + log-rank + HR\n")
km_data <- data.frame(time = os_time, event = os_event, RiskGroup = risk_group)
fit_km <- survfit(Surv(time, event) ~ RiskGroup, data = km_data)
sd_lr <- survdiff(Surv(time, event) ~ RiskGroup, data = km_data)
lr_p <- 1 - pchisq(sd_lr$chisq, length(sd_lr$n) - 1)
fit_cox_rg <- coxph(Surv(time, event) ~ RiskGroup, data = km_data)
hr_summary <- summary(fit_cox_rg)
HR <- hr_summary$conf.int[1, "exp(coef)"]
HR_low <- hr_summary$conf.int[1, "lower .95"]
HR_high <- hr_summary$conf.int[1, "upper .95"]
cox_p <- hr_summary$coefficients[1, "Pr(>|z|)"]
med_os <- summary(fit_km)$table
med_os_low <- med_os["RiskGroup=Low", "median"]
med_os_high <- med_os["RiskGroup=High", "median"]
cat(sprintf(" log-rank p : %.4g\n", lr_p))
cat(sprintf(" HR (High vs Low): %.3f (95%% CI %.3f %.3f, Cox p=%.4g)\n",
            HR, HR_low, HR_high, cox_p))
fwrite(data.table(n_low = sum(risk_group == "Low"), n_high = sum(risk_group == "High"),
                  median_OS_low_days = med_os_low, median_OS_high_days = med_os_high,
                  HR = HR, HR_low_95 = HR_low, HR_high_95 = HR_high,
                  cox_p_value = cox_p, logrank_p = lr_p),
       file.path(TBLDIR, "Tab10_KM_logrank_results.csv"))

cat("\n[8/14] Time-dependent ROC\n")
at_risk_5yr <- sum(os_time >= 1825)
cat(sprintf(" at-risk @ 5-yr : %d\n", at_risk_5yr))
if (at_risk_5yr < AT_RISK_MIN) {
  cat(sprintf(" [WARN] fallback to 1/2/3-yr (5-yr at-risk %d < %d)\n", at_risk_5yr, AT_RISK_MIN))
  roc_times_used <- ROC_FALLBACK
  roc_labels <- c("1-yr","2-yr","3-yr")
} else {
  roc_times_used <- ROC_TIMES; roc_labels <- c("1-yr","3-yr","5-yr")
}
roc_res <- timeROC(T = os_time, delta = os_event, marker = risk_score, cause = 1,
                   times = roc_times_used, iid = TRUE)
auc_vec <- as.numeric(roc_res$AUC)
# 2026-10-06: vect_sd_1 already is the standard error of the AUC; the earlier sqrt() made the 95% CI
# far too wide (0.337-1.000 for an AUC of 0.755). Corrected table: helpers/11_fix_cindex_timeROC_tables.R
sd_vec <- tryCatch(as.numeric(roc_res$inference$vect_sd_1),
                    error = function(e) rep(NA_real_, length(auc_vec)))
auc_low <- pmax(0, auc_vec - 1.96 * sd_vec)
auc_high <- pmin(1, auc_vec + 1.96 * sd_vec)
roc_tbl <- data.table(time_label = roc_labels, time_days = roc_times_used,
                      AUC = auc_vec, AUC_low_95 = auc_low,
                      AUC_high_95 = auc_high, SE = sd_vec)
fwrite(roc_tbl, file.path(TBLDIR, "Tab10_timeROC_AUC.csv"))
print(roc_tbl)

# Fig 6A LASSO coefficient path (with par(mar) cosmetic fix)
cat("\n[9/14] Fig6A LASSO coefficient path\n")
draw_6A <- function() {
  # Lambda lines also removed (lambda values preserved in log).
  par(mar = c(5, 4, 2, 2))
  par(oma = c(0, 0, 0, 0))
  plot(cv_fit$glmnet.fit, xvar = "lambda", label = FALSE, main = "")
}
cairo_pdf(file.path(FIGDIR, "Fig6A_LASSO_coefficient_path.pdf"), width = 7, height = 5.5)
draw_6A(); dev.off()
ragg::agg_png(file.path(FIGDIR, "Fig6A_LASSO_coefficient_path.png"), width = 7, height = 5.5,
              units = "in", res = 600)
draw_6A(); dev.off()
cat(" [OK] Fig6A_LASSO_coefficient_path.{pdf,png}\n")

# Fig 6B CV deviance (with par(mar) fix)
cat("\n[10/14] Fig6B CV deviance plot\n")
draw_6B <- function() {
  par(mar = c(5, 4, 2, 2))
  par(oma = c(0, 0, 0, 0))
  plot(cv_fit, main = "")
}
cairo_pdf(file.path(FIGDIR, "Fig6B_CV_deviance_plot.pdf"), width = 7, height = 5.5)
draw_6B(); dev.off()
ragg::agg_png(file.path(FIGDIR, "Fig6B_CV_deviance_plot.png"), width = 7, height = 5.5,
              units = "in", res = 600)
draw_6B(); dev.off()
cat(" [OK] Fig6B_CV_deviance_plot.{pdf,png}\n")

# Fig 6C KM
cat("\n[11/14] Fig6C KM by risk group\n")
# HR/CI/p annotation moved from title to upper-right INSIDE the plot panel.
ss <- ggsurvplot(fit_km, data = km_data, pval = FALSE, risk.table = TRUE, conf.int = FALSE,
                 palette = c(unname(RISK_COL["Low"]), unname(RISK_COL["High"])),
                 legend.labs = c("Low risk","High risk"), legend.title = "Risk group",
                 xlab = "Time (days)", ylab = "Overall survival",
                 title = NULL,
                 ggtheme = theme_bw(base_size = 11), risk.table.fontsize = 3.5)
hr_lbl_6C <- sprintf("HR = %.2f (95%% CI %.2f - %.2f)\nlog-rank p = %.3g",
                     HR, HR_low, HR_high, lr_p)
ss$plot <- ss$plot +
  annotate("label", x = Inf, y = Inf, label = hr_lbl_6C,
           hjust = 1.05, vjust = 1.2, size = 3.6, fontface = "plain",
           label.size = 0.3, fill = "white", label.r = unit(0, "pt"))
combined_6C <- cowplot::plot_grid(ss$plot + theme(legend.position = "top"),
                                   ss$table + theme(plot.title = element_blank()),
                                   ncol = 1, align = "v", axis = "lr",
                                   rel_heights = c(3, 1))
ggsave(file.path(FIGDIR, "Fig6C_KM_by_risk_group.pdf"),
       plot = combined_6C, width = 8, height = 7, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig6C_KM_by_risk_group.png"),
       plot = combined_6C, width = 8, height = 7, dpi = 600, device = ragg::agg_png)
cat(" [OK] Fig6C_KM_by_risk_group.{pdf,png}\n")

# Fig 6D timeROC
cat("\n[12/14] Fig6D time-dependent ROC\n")
roc_dfs <- lapply(seq_along(roc_times_used), function(k) {
  data.frame(FPR = roc_res$FP[, k], TPR = roc_res$TP[, k],
             time_label = sprintf("%s (AUC=%.3f, 95%% CI %.3f-%.3f)",
                                  roc_labels[k], auc_vec[k], auc_low[k], auc_high[k]))
})
roc_df <- do.call(rbind, roc_dfs)
palette_used <- c("#377EB8","#4DAF4A","#E41A1C")
names(palette_used) <- unique(roc_df$time_label)
p6d <- ggplot(roc_df, aes(x = FPR, y = TPR, color = time_label)) +
  geom_line(linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = palette_used, name = NULL) + coord_equal() +
  labs(x = "False positive rate", y = "True positive rate") +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.65, 0.20),
        legend.background = element_rect(fill = "white", color = "grey70"))
ggsave(file.path(FIGDIR, "Fig6D_timeROC.pdf"), plot = p6d,
       width = 6.5, height = 6, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig6D_timeROC.png"), plot = p6d,
       width = 6.5, height = 6, dpi = 600, device = ragg::agg_png)
cat(" [OK] Fig6D_timeROC.{pdf,png}\n")

# Fig 6E 3-panel risk distribution
cat("\n[13/14] Fig6E 3-panel risk distribution\n")
ord <- order(risk_score)
rs_sorted <- risk_score[ord]; samples_sort <- samples_keep[ord]
rg_sorted <- as.character(risk_group)[ord]
os_event_s <- os_event[ord]; os_time_s <- os_time[ord]
N <- length(ord)

p1 <- ggplot(data.frame(rank = 1:N, RS = rs_sorted,
                         RiskGroup = factor(rg_sorted, levels = c("Low","High"))),
             aes(x = rank, y = RS, color = RiskGroup)) +
  geom_point(size = 0.8) +
  geom_hline(yintercept = med_cutoff, linetype = "dashed", color = "grey40") +
  annotate("text", x = N * 0.05, y = med_cutoff,
           label = sprintf("median = %.4f", med_cutoff), vjust = -0.5, size = 3, color = "grey40") +
  scale_color_manual(values = c(Low = unname(RISK_COL["Low"]), High = unname(RISK_COL["High"])),
                     name = "Risk group") +
  labs(x = NULL, y = "Risk Score") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right", axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), panel.grid.minor = element_blank())

p2 <- ggplot(data.frame(rank = 1:N, OS_time = os_time_s,
                         Event = factor(ifelse(os_event_s == 1, "Death", "Censored"),
                                         levels = c("Censored","Death"))),
             aes(x = rank, y = OS_time, color = Event)) +
  geom_point(size = 0.9, alpha = 0.8) +
  scale_color_manual(values = c(Censored = "grey30", Death = "#E41A1C"), name = "OS event") +
  labs(x = NULL, y = "Time (days)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right", axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), panel.grid.minor = element_blank())

heat_mat <- cand_mat[selected_genes, samples_sort, drop = FALSE]
heat_z <- t(scale(t(heat_mat), center = TRUE, scale = TRUE))
heat_dt <- as.data.table(heat_z, keep.rownames = "gene")
heat_long <- melt(heat_dt, id.vars = "gene", variable.name = "sample", value.name = "z")
heat_long[, sample := factor(sample, levels = samples_sort)]
heat_long[, gene := factor(gene, levels = rev(selected_genes))]

p3 <- ggplot(heat_long, aes(x = sample, y = gene, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#0571B0", mid = "white", high = "#CA0020",
                        midpoint = 0, name = "z-score",
                        limits = c(-3, 3), oob = scales::squish) +
  labs(x = sprintf("Samples sorted by Risk Score (low high, n=%d)", N), y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 7), panel.grid = element_blank(),
        legend.position = "right")

combined_6E <- cowplot::plot_grid(p1, p2, p3, ncol = 1, align = "v", axis = "lr",
                                   rel_heights = c(1.1, 1.1,
                                                    max(1.5, 0.18 * length(selected_genes) + 0.6)))
H_6E <- 3 + 0.18 * length(selected_genes) + 1
ggsave(file.path(FIGDIR, "Fig6E_risk_distribution_3panel.pdf"),
       plot = combined_6E, width = 11, height = H_6E, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(FIGDIR, "Fig6E_risk_distribution_3panel.png"),
       plot = combined_6E, width = 11, height = H_6E, dpi = 600,
       device = ragg::agg_png, limitsize = FALSE)
cat(" [OK] Fig6E_risk_distribution_3panel.{pdf,png}\n")

# Summary
t_end <- Sys.time()
elapsed_total <- as.numeric(t_end - t_start, units = "mins")
cat("\n========================================================================\n")
cat(" Stage 10 LASSO-Cox Summary Report\n")
cat("========================================================================\n")
cat(sprintf(" Sample alignment : %d (expr Primary Tumor survival)\n", N_train))
cat(sprintf(" Candidate pool : %d (Stage 9 PC-only sig DEG)\n", nrow(deg_dt)))
cat(sprintf(" Univariate Cox p<0.01 : %d genes\n", nrow(uni_p001)))
cat(sprintf(" LASSO selected : %d genes (all protein_coding by upstream filter)\n", n_sel))
cat(sprintf(" Median RS cutoff : %.4f\n", med_cutoff))
cat(sprintf(" Low / High n : %d / %d\n",
            sum(risk_group == "Low"), sum(risk_group == "High")))
cat(sprintf(" log-rank p : %.4g\n", lr_p))
cat(sprintf(" HR (High vs Low) : %.3f (95%% CI %.3f %.3f)\n", HR, HR_low, HR_high))
cat("\n AUC table:\n")
print(roc_tbl)
cat(sprintf("\n ROC times used: %s (5-yr fallback: %s)\n",
            paste(roc_labels, collapse = ", "),
            ifelse(identical(roc_times_used, ROC_FALLBACK), "YES", "NO")))
cat(sprintf(" at-risk @ 5-yr: %d\n", at_risk_5yr))
cat(sprintf("\n Risk Score formula:\n RS = %s\n", formula_str))
cat("\n Selected genes:\n"); print(sel_dt)
cat(sprintf("\n Total elapsed: %.2f min\n", elapsed_total))
cat("========================================================================\n")

sess_file <- file.path(LOGDIR, "sessionInfo/sessionInfo_stage10.txt")
sink(sess_file); cat(sprintf("Session info captured: %s\n\n", Sys.time())); print(sessionInfo()); sink()
cat(sprintf(">>> sessionInfo: %s\n", sess_file))

cat(sprintf("\n=== Stage 10 done (%s, %.2f min) ===\n", Sys.time(), elapsed_total))
sink(type = "message"); sink(); close(log_con)

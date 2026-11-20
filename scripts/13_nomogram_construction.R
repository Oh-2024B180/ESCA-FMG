###############################################################################
## Stage 13: Nomogram (TCGA-ESCA, RS + clinical variables)
##
## 13A: Variable selection
## - Candidates: RiskScore + age + gender + pathologic_stage + grade
## + residual_tumor + histology
## - Multivariate Cox with all seven candidates (complete cases)
## - Final nomogram variables: **RS + histology + stage + residual**
## 13B: rms::cph + nomogram + calibrate (Harrell 2015 textbook)
## 13C: Calibration bootstrap B=1000 at 1, 2 and 3 years (5 years only if at least AT_RISK_MIN patients remain at risk)
## (5-yr at-risk = 4 < 10, so 1/2/3-yr only)
## 13D: C-index comparison: nomogram / RS alone / clinical only
##
## Inputs:
## data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv
## tables/Tab10_risk_scores_per_sample.csv
##
## Outputs:
## Tab13_univariate_Cox_clinical.csv (per-variable univariate)
## Tab13_multivariate_Cox.csv (full 7-variable multivariate)
## Tab13_nomogram_summary.csv (final 4-variable model coefs)
## Tab13_C_index_comparison.csv (nomogram / RS / clinical)
## Tab13_nomogram_timeROC_AUC.csv (1/2/3-year nomogram AUC)
## Fig9A_nomogram.{pdf,png}
## Fig9B_calibration.{pdf,png}
## Fig9C_C_index_comparison.{pdf,png}
## logs/13_nomogram_diagnostics.txt + sessionInfo_stage13.txt
##
## Random seed: 1262118388 (CLAUDE.md §4.3)
###############################################################################

# run from the repository root (see README.md)

INPUT_CLI <- "data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv"
INPUT_RS <- "tables/Tab10_risk_scores_per_sample.csv"

TBLDIR <- "tables"; FIGDIR <- "figures"; LOGDIR <- "logs"
for (d in c(TBLDIR, FIGDIR, LOGDIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

CALIB_B <- 1000
# 1-/2-/3-yr (5-yr fallback per low at-risk; 2-yr added 2026-05-07 for
# consistency with Results body "1-, 2- and 3-year OS probability" framing and Fig9D
# timeROC. Verified at-risk @ 2-yr = 36 (>= AT_RISK_MIN); 3-yr = 14; 5-yr = 4.)
CALIB_TIMES_DAYS <- c(365, 730, 1095)
AT_RISK_MIN <- 10

set.seed(1262118388)
t_start <- Sys.time()

log_file <- file.path(LOGDIR, "13_diagnostics.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE); sink(log_con, type = "message")
cat(sprintf("=== Stage 13 Nomogram start (%s) ===\n\n", Sys.time()))

cat("[1/12] Loading packages\n")
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(rms)
  library(ggplot2)
  library(ragg)
  library(cowplot)
})
cat(sprintf(" rms : %s\n", packageVersion("rms")))
cat(sprintf(" Hmisc : %s\n", packageVersion("Hmisc")))

cat("\n[2/12] Loading clinical + Stage 10 RS\n")
cli_dt <- fread(INPUT_CLI)
rs_dt <- fread(INPUT_RS)
cli_dt <- cli_dt[sampleID %in% rs_dt$sampleID]
cat(sprintf(" clinical rows matching RS: %d\n", nrow(cli_dt)))

# 3. Build merged pheno + consolidate variables
cat("\n[3/12] Variable consolidation\n")

# Stage: collapse sub-stages (IA/IB I, IIA/IIB II, etc.)
collapse_stage <- function(s) {
  s[s == "" | grepl("^\\[", s)] <- NA
  out <- rep(NA_character_, length(s))
  out[grepl("Stage I[ABC]$", s)] <- "I"
  out[grepl("Stage II[ABC]$", s)] <- "II"
  out[grepl("Stage III[ABC]$", s)] <- "III"
  out[grepl("Stage IV[ABC]$", s)] <- "IV"
  out[s == "Stage I"] <- "I"
  out[s == "Stage II"] <- "II"
  out[s == "Stage III"] <- "III"
  out[s == "Stage IV"] <- "IV"
  out
}

# Residual: R0 vs R1+R2 (collapse R1+R2 due to small N), exclude RX
collapse_residual <- function(r) {
  r[r == "" | grepl("^\\[", r)] <- NA
  out <- rep(NA_character_, length(r))
  out[r == "R0"] <- "R0"
  out[r %in% c("R1","R2")] <- "R1+R2"
  out[r == "RX"] <- NA
  out
}

# Grade: G1/G2/G3, exclude GX
collapse_grade <- function(g) {
  g[g == "" | grepl("^\\[", g)] <- NA
  out <- rep(NA_character_, length(g))
  out[g %in% c("G1","G2","G3")] <- g[g %in% c("G1","G2","G3")]
  out[g == "GX"] <- NA
  out
}

# Histology
collapse_histology <- function(h) {
  ifelse(h == "Esophagus Squamous Cell Carcinoma", "ESCC",
  ifelse(h == "Esophagus Adenocarcinoma, NOS", "EAC", NA))
}

cli_dt[, stage_collapsed := collapse_stage(pathologic_stage)]
cli_dt[, residual_collapsed := collapse_residual(residual_tumor)]
cli_dt[, grade_collapsed := collapse_grade(neoplasm_histologic_grade)]
cli_dt[, histology_collapsed := collapse_histology(histological_type)]
cli_dt[, age := age_at_initial_pathologic_diagnosis]
cli_dt[, sex := gender]

# Merge with RS
ph <- merge(rs_dt[, .(sampleID, OS, OS.time, RiskScore)],
            cli_dt[, .(sampleID, age, sex, stage_collapsed, residual_collapsed,
                       grade_collapsed, histology_collapsed)],
            by = "sampleID")
cat(sprintf(" merged rows: %d\n", nrow(ph)))

# Variable distribution after collapse
cat("\n Stage: "); print(table(ph$stage_collapsed, useNA = "always"))
cat("\n Residual: "); print(table(ph$residual_collapsed, useNA = "always"))
cat("\n Grade: "); print(table(ph$grade_collapsed, useNA = "always"))
cat("\n Histology: "); print(table(ph$histology_collapsed, useNA = "always"))

# 4. Univariate Cox per candidate variable
cat("\n[4/12] Univariate Cox per candidate (n varies due to NA)\n")

uni_run <- function(var, ph) {
  df <- ph[!is.na(get(var)) & !is.na(OS) & !is.na(OS.time) & OS.time > 0]
  if (is.character(df[[var]])) df[, (var) := factor(get(var))]
  fit <- tryCatch(coxph(as.formula(sprintf("Surv(OS.time, OS) ~ %s", var)),
                         data = df), error = function(e) NULL)
  if (is.null(fit)) return(data.table(variable = var, level = NA, n = nrow(df),
                                       HR = NA, HR_low = NA, HR_high = NA, p_value = NA))
  s <- summary(fit)
  rows <- lapply(seq_len(nrow(s$coefficients)), function(i) {
    data.table(variable = var, level = rownames(s$coefficients)[i],
               n = s$n, HR = s$coefficients[i, "exp(coef)"],
               HR_low = s$conf.int[i, "lower .95"],
               HR_high = s$conf.int[i, "upper .95"],
               p_value = s$coefficients[i, "Pr(>|z|)"])
  })
  rbindlist(rows)
}

vars_to_test <- c("RiskScore", "age", "sex", "stage_collapsed",
                  "residual_collapsed", "grade_collapsed", "histology_collapsed")
uni_dt <- rbindlist(lapply(vars_to_test, uni_run, ph = ph))
fwrite(uni_dt, file.path(TBLDIR, "Tab13_univariate_Cox_clinical.csv"))
print(uni_dt)

# 5. Multivariate Cox (full 7-variable, complete-case)
cat("\n[5/12] Multivariate Cox (full 7-variable model)\n")
ph_full <- ph[complete.cases(ph[, .(RiskScore, age, sex, stage_collapsed,
                                     residual_collapsed, grade_collapsed,
                                     histology_collapsed)])]
ph_full[, sex := factor(sex, levels = c("MALE","FEMALE"))]
ph_full[, stage_collapsed := factor(stage_collapsed, levels = c("I","II","III","IV"))]
ph_full[, residual_collapsed := factor(residual_collapsed, levels = c("R0","R1+R2"))]
ph_full[, grade_collapsed := factor(grade_collapsed, levels = c("G1","G2","G3"))]
ph_full[, histology_collapsed := factor(histology_collapsed, levels = c("ESCC","EAC"))]
cat(sprintf(" complete-case n: %d\n", nrow(ph_full)))
cat(sprintf(" complete-case deaths: %d\n", sum(ph_full$OS)))
print(ph_full[, .(n = .N, deaths = sum(OS)), by = sex])          # patients and deaths by sex (Results 3.8)

mv_full <- coxph(Surv(OS.time, OS) ~ RiskScore + age + sex + stage_collapsed +
                  residual_collapsed + grade_collapsed + histology_collapsed,
                  data = ph_full)
mv_summary <- summary(mv_full)
mv_dt <- data.table(
  term = rownames(mv_summary$coefficients),
  HR = mv_summary$coefficients[, "exp(coef)"],
  HR_low = mv_summary$conf.int[, "lower .95"],
  HR_high = mv_summary$conf.int[, "upper .95"],
  p_value = mv_summary$coefficients[, "Pr(>|z|)"]
)
fwrite(mv_dt, file.path(TBLDIR, "Tab13_multivariate_Cox.csv"))
cat(" Multivariate Cox results:\n"); print(mv_dt)

# Identify p<0.05 variables
sig_terms <- mv_dt[p_value < 0.05, term]
cat(sprintf("\n Significant (p<0.05) terms in multivariate: %d\n", length(sig_terms)))
print(sig_terms)

# 6. Final nomogram model (four variables)
cat("\n[6/12] Final nomogram model RS + histology + stage + residual\n")
ph_nom <- ph[complete.cases(ph[, .(RiskScore, stage_collapsed,
                                     residual_collapsed, histology_collapsed)])]
ph_nom[, stage_collapsed := factor(stage_collapsed, levels = c("I","II","III","IV"))]
ph_nom[, residual_collapsed := factor(residual_collapsed, levels = c("R0","R1+R2"))]
ph_nom[, histology_collapsed := factor(histology_collapsed, levels = c("ESCC","EAC"))]
cat(sprintf(" N (nomogram CC): %d\n", nrow(ph_nom)))

# rms::cph requires datadist set globally
ph_nom_df <- as.data.frame(ph_nom)
dd <- datadist(ph_nom_df); options(datadist = "dd")

nom_fit <- cph(Surv(OS.time, OS) ~ RiskScore + histology_collapsed +
                stage_collapsed + residual_collapsed,
               data = ph_nom_df, x = TRUE, y = TRUE, surv = TRUE,
               time.inc = 365)

cat(" cph fit summary:\n")
print(nom_fit)

# Build summary table
nom_coefs <- summary(nom_fit)
coef_vec <- coef(nom_fit)
se_vec <- sqrt(diag(vcov(nom_fit)))
z_vec <- as.numeric(coef_vec) / se_vec
nom_summary <- data.table(
  term = names(coef_vec),
  coef = as.numeric(coef_vec),
  HR = as.numeric(exp(coef_vec)),
  se = se_vec,
  z = z_vec,
  p_value = 2 * pnorm(-abs(z_vec))
)
fwrite(nom_summary, file.path(TBLDIR, "Tab13_nomogram_summary.csv"))
cat(" nomogram model coefficients:\n"); print(nom_summary)

# 7. Build nomogram + Fig8A
cat("\n[7/12] Fig8A nomogram\n")
# Survival functions for nomogram
surv_obj <- Survival(nom_fit)
surv_1yr <- function(x) surv_obj(365, x)
surv_2yr <- function(x) surv_obj(730, x)
surv_3yr <- function(x) surv_obj(1095, x)

nom <- nomogram(nom_fit,
                fun = list(surv_1yr, surv_2yr, surv_3yr),
                funlabel = c("1-year survival", "2-year survival",
                             "3-year survival"),
                lp = TRUE,
                fun.at = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 0.95))

W <- 11; H <- 8
draw_nom <- function() {
  par(mar = c(3, 4, 1, 1))
  plot(nom, xfrac = 0.30, cex.axis = 0.85, cex.var = 0.95,
       lmgp = 0.2)
}
cairo_pdf(file.path(FIGDIR, "Fig9A_nomogram.pdf"), width = W, height = H)
draw_nom(); dev.off()
ragg::agg_png(file.path(FIGDIR, "Fig9A_nomogram.png"), width = W, height = H,
              units = "in", res = 600)
draw_nom(); dev.off()
cat(" [OK] Fig9A_nomogram.{pdf,png}\n")

# 8. Calibration plots (1-/2-/3-yr; 5-yr fallback if at-risk >= 10)
cat("\n[8/12] Fig9B calibration (1-/2-/3-yr, B=1000 bootstrap)\n")

# At-risk check at each candidate time (only 5-yr can be skipped per fallback)
for (t_d in c(365, 730, 1095, 1825)) {
  cat(sprintf("   at-risk @ %.0f-yr: %d\n", t_d / 365,
              sum(ph_nom_df$OS.time >= t_d)))
}
at_risk_5yr <- sum(ph_nom_df$OS.time >= 1825)
cat(sprintf(" 5-yr fallback rule: skip 5-yr calibration if < %d\n",
            AT_RISK_MIN))
if (at_risk_5yr < AT_RISK_MIN) {
  times_calib <- c(365, 730, 1095)
  labels_calib <- c("1-year","2-year","3-year")
} else {
  times_calib <- c(365, 730, 1095, 1825)
  labels_calib <- c("1-year","2-year","3-year","5-year")
}
cat(sprintf(" calibration times: %s\n", paste(labels_calib, collapse = ", ")))

# Need to refit cph with explicit time.inc per calibration time
draw_calib <- function() {
  par(mfrow = c(1, length(times_calib)),
      mar = c(4, 4, 4, 1), oma = c(0, 0, 2, 0))
  for (i in seq_along(times_calib)) {
    t_d <- times_calib[i]
    fit_t <- cph(Surv(OS.time, OS) ~ RiskScore + histology_collapsed +
                  stage_collapsed + residual_collapsed,
                  data = ph_nom_df, x = TRUE, y = TRUE, surv = TRUE,
                  time.inc = t_d)
    cal <- calibrate(fit_t, cmethod = "KM", method = "boot", u = t_d, m = 50,
                     B = CALIB_B)
    # Per-panel main label retained as a small horizon (1-year/3-year etc.) -
    # this acts as a sub-panel identifier, not a figure title.
    plot(cal, xlim = c(0, 1), ylim = c(0, 1),
         xlab = "Nomogram predicted survival",
         ylab = "Observed survival",
         main = labels_calib[i])
    abline(0, 1, lty = 3, col = "grey50")
  }
}

W2 <- 6 * length(times_calib); H2 <- 6
cairo_pdf(file.path(FIGDIR, "Fig9B_calibration.pdf"), width = W2, height = H2)
draw_calib(); dev.off()
ragg::agg_png(file.path(FIGDIR, "Fig9B_calibration.png"), width = W2, height = H2,
              units = "in", res = 600)
draw_calib(); dev.off()
cat(" [OK] Fig9B_calibration.{pdf,png}\n")

# 9. C-index comparison: nomogram / RS / clinical only
cat("\n[9/12] C-index comparison: nomogram vs RS alone vs clinical only\n")

ph_cm <- ph_nom_df # complete-case
fit_nom_c <- coxph(Surv(OS.time, OS) ~ RiskScore + histology_collapsed +
                    stage_collapsed + residual_collapsed, data = ph_cm)
fit_rs_c <- coxph(Surv(OS.time, OS) ~ RiskScore, data = ph_cm)
fit_cli_c <- coxph(Surv(OS.time, OS) ~ histology_collapsed +
                    stage_collapsed + residual_collapsed, data = ph_cm)

c_extract <- function(fit, label) {
  s <- summary(fit)
  c_idx <- s$concordance["C"]
  c_se <- s$concordance["se(C)"]
  data.table(model = label, c_index = c_idx,
             c_low_95 = c_idx - 1.96 * c_se,
             c_high_95 = c_idx + 1.96 * c_se, c_se = c_se)
}
cidx_dt <- rbindlist(list(
  c_extract(fit_nom_c, "Nomogram (RS + histology + stage + residual)"),
  c_extract(fit_rs_c, "RS alone"),
  c_extract(fit_cli_c, "Clinical only (histology + stage + residual)")
))
fwrite(cidx_dt, file.path(TBLDIR, "Tab13_C_index_comparison.csv"))
cat(" c-index comparison:\n"); print(cidx_dt)

# Fig8C bar
cidx_dt[, model := factor(model, levels = model)]
p8C <- ggplot(cidx_dt, aes(x = model, y = c_index, fill = model)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = c_low_95, ymax = c_high_95),
                width = 0.2, color = "black", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", c_index)),
            vjust = -1.5, size = 4, fontface = "bold") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.10))) +
  scale_fill_manual(values = c("#E41A1C","#377EB8","#4DAF4A"), guide = "none") +
  labs(x = NULL, y = "Concordance index") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(face = "bold", size = 8))

ggsave(file.path(FIGDIR, "Fig9C_C_index_comparison.pdf"),
       plot = p8C, width = 9, height = 5, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig9C_C_index_comparison.png"),
       plot = p8C, width = 9, height = 5, dpi = 600, device = ragg::agg_png)
cat(" [OK] Fig9C_C_index_comparison.{pdf,png}\n")

# 10. Time-dependent ROC AUC for nomogram (1/2/3-year)
# Reference: Blanche P 2013 Stat Med — timeROC (PMID 24027076)
cat("\n[10/12] Time-dependent ROC AUC (nomogram linear predictor)\n")
suppressPackageStartupMessages(library(timeROC))
nom_lp <- predict(nom_fit, type = "lp")
tr_nom <- timeROC(T = ph_nom_df$OS.time, delta = ph_nom_df$OS,
                  marker = nom_lp, cause = 1,
                  times = c(365, 730, 1095), iid = TRUE)
auc_vec <- tr_nom$AUC
se_vec_tr <- tr_nom$inference$vect_sd_1
nom_roc_dt <- data.table(
  time_label = c("1-yr", "2-yr", "3-yr"),
  time_days = c(365, 730, 1095),
  AUC = auc_vec,
  AUC_low_95 = pmax(0, auc_vec - 1.96 * se_vec_tr),
  AUC_high_95 = pmin(1, auc_vec + 1.96 * se_vec_tr),
  SE = se_vec_tr
)
fwrite(nom_roc_dt, file.path(TBLDIR, "Tab13_nomogram_timeROC_AUC.csv"))
cat(" Nomogram timeROC AUC (1/2/3-yr):\n"); print(nom_roc_dt)

# 11. Summary
t_end <- Sys.time(); elapsed <- as.numeric(t_end - t_start, units = "mins")
cat("\n========================================================================\n")
cat(" Stage 13 Nomogram Summary Report\n")
cat("========================================================================\n")
cat(sprintf(" Multivariate Cox CC (full 7 vars) : %d\n", nrow(ph_full)))
cat(sprintf(" Nomogram model CC (4 vars) : %d\n", nrow(ph_nom_df)))
cat(sprintf(" Calibration bootstrap B : %d\n", CALIB_B))
cat(sprintf(" Calibration time points : %s\n",
            paste(labels_calib, collapse = ", ")))
cat(sprintf(" 5-yr fallback : %s\n",
            ifelse(at_risk_5yr < AT_RISK_MIN, "YES (at-risk < 10)", "NO")))

cat("\n --- Nomogram model coefficients ---\n")
print(nom_summary)
cat("\n --- C-index comparison ---\n")
print(cidx_dt)
cat("\n --- Nomogram time-dependent AUC (1/2/3-yr) ---\n")
print(nom_roc_dt)
cat(sprintf("\n Total elapsed: %.2f min\n", elapsed))
cat("========================================================================\n")

# sessionInfo
sess_file <- file.path(LOGDIR, "sessionInfo/sessionInfo_stage13.txt")
sink(sess_file); cat(sprintf("Session info: %s\n\n", Sys.time())); print(sessionInfo()); sink()
cat(sprintf("\n>>> sessionInfo: %s\n", sess_file))
cat(sprintf("\n=== Stage 13 Nomogram done (%s) ===\n", Sys.time()))
sink(type = "message"); sink(); close(log_con)

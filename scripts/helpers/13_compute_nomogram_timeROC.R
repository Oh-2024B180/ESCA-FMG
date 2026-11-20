###############################################################################
## Stage 13 helper — Compute nomogram time-dependent AUC (1/2/3-year)
##
## Inputs:
##   - data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv
##   - tables/Tab10_risk_scores_per_sample.csv
##
## Outputs:
##   - tables/Tab13_nomogram_timeROC_AUC.csv
##   - logs/fix_logs/13_compute_nomogram_timeROC.log
##
## Random seed: 1262118388 (CLAUDE.md §4.3)
###############################################################################

# run from the repository root (see README.md)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(rms)
  library(timeROC)
})

set.seed(1262118388)

LOGDIR <- "logs/fix_logs"
if (!dir.exists(LOGDIR)) dir.create(LOGDIR, recursive = TRUE)

log_file <- file.path(LOGDIR, "13_compute_nomogram_timeROC.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE); sink(log_con, type = "message")

cat(sprintf("=== Stage 13 nomogram timeROC start (%s) ===\n\n", Sys.time()))

# Reproduce nomogram cohort (same logic as scripts/13_nomogram_construction.R)
INPUT_CLI <- "data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv"
INPUT_RS <- "tables/Tab10_risk_scores_per_sample.csv"

cli_dt <- fread(INPUT_CLI)
rs_dt <- fread(INPUT_RS)
cli_dt <- cli_dt[sampleID %in% rs_dt$sampleID]

collapse_stage <- function(s) {
  s[s == "" | grepl("^\\[", s)] <- NA
  out <- rep(NA_character_, length(s))
  out[grepl("Stage I[ABC]$", s)] <- "I"
  out[grepl("Stage II[ABC]$", s)] <- "II"
  out[grepl("Stage III[ABC]$", s)] <- "III"
  out[grepl("Stage IV[ABC]$", s)] <- "IV"
  out[s == "Stage I"] <- "I"; out[s == "Stage II"] <- "II"
  out[s == "Stage III"] <- "III"; out[s == "Stage IV"] <- "IV"
  out
}
collapse_residual <- function(r) {
  r[r == "" | grepl("^\\[", r)] <- NA
  out <- rep(NA_character_, length(r))
  out[r == "R0"] <- "R0"; out[r %in% c("R1","R2")] <- "R1+R2"
  out[r == "RX"] <- NA
  out
}
collapse_histology <- function(h) {
  ifelse(h == "Esophagus Squamous Cell Carcinoma", "ESCC",
         ifelse(h == "Esophagus Adenocarcinoma, NOS", "EAC", NA))
}

cli_dt[, stage_collapsed := collapse_stage(pathologic_stage)]
cli_dt[, residual_collapsed := collapse_residual(residual_tumor)]
cli_dt[, histology_collapsed := collapse_histology(histological_type)]

ph <- merge(rs_dt[, .(sampleID, OS, OS.time, RiskScore)],
            cli_dt[, .(sampleID, stage_collapsed, residual_collapsed,
                       histology_collapsed)],
            by = "sampleID")

ph_nom <- ph[complete.cases(ph[, .(RiskScore, stage_collapsed,
                                    residual_collapsed, histology_collapsed)])]
ph_nom[, stage_collapsed := factor(stage_collapsed, levels = c("I","II","III","IV"))]
ph_nom[, residual_collapsed := factor(residual_collapsed, levels = c("R0","R1+R2"))]
ph_nom[, histology_collapsed := factor(histology_collapsed, levels = c("ESCC","EAC"))]

cat(sprintf("Nomogram cohort N (CC): %d (expected 145)\n", nrow(ph_nom)))
stopifnot(nrow(ph_nom) == 145)

ph_nom_df <- as.data.frame(ph_nom)
dd <- datadist(ph_nom_df); options(datadist = "dd")

nom_fit <- cph(Surv(OS.time, OS) ~ RiskScore + histology_collapsed +
                 stage_collapsed + residual_collapsed,
               data = ph_nom_df, x = TRUE, y = TRUE, surv = TRUE,
               time.inc = 365)

# Linear predictor (centered) — use raw lp from cph
lp <- predict(nom_fit, type = "lp")

# timeROC at 1/2/3-year
tr <- timeROC(T = ph_nom_df$OS.time, delta = ph_nom_df$OS,
              marker = lp, cause = 1,
              times = c(365, 730, 1095),
              iid = TRUE)

# Build output table with 95% CI
auc_vals <- tr$AUC
se_vals <- tr$inference$vect_sd_1
ci_low <- pmax(0, auc_vals - 1.96 * se_vals)
ci_high <- pmin(1, auc_vals + 1.96 * se_vals)

out_dt <- data.table(
  time_label = c("1-yr", "2-yr", "3-yr"),
  time_days = c(365, 730, 1095),
  AUC = auc_vals,
  AUC_low_95 = ci_low,
  AUC_high_95 = ci_high,
  SE = se_vals
)

cat("\nNomogram timeROC AUC results:\n")
print(out_dt)

OUT_CSV <- "tables/Tab13_nomogram_timeROC_AUC.csv"
fwrite(out_dt, OUT_CSV)
cat(sprintf("\n[OK] %s\n", OUT_CSV))

cat(sprintf("\n=== End (%s) ===\n", Sys.time()))
sink(NULL); sink(NULL, type = "message")
close(log_con)

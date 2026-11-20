###############################################################################
## Helper: corrected discrimination tables for Stage 10 and Stage 11
##
## Inputs : tables/Tab10_risk_scores_per_sample.csv, tables/Tab11_GSE53625_RS_per_sample.csv
##
## Outputs: tables/Tab10_timeROC_AUC.csv, tables/Tab11_GSE53625_timeROC_AUC.csv,
##          tables/Tab11_GSE53625_concordance_index.csv
###############################################################################
suppressPackageStartupMessages({library(survival); library(timeROC); library(data.table)})
ROOT <- ""   # repository root = working directory (see README.md)
TBL <- paste0(ROOT, "tables/")
tr <- fread(paste0(TBL, "Tab10_risk_scores_per_sample.csv"))
va <- fread(paste0(TBL, "Tab11_GSE53625_RS_per_sample.csv"))

auc_table <- function(time, event, marker, times, labels) {
  r <- timeROC(T = time, delta = event, marker = marker, cause = 1, times = times, iid = TRUE)
  se <- as.numeric(r$inference$vect_sd_1)                      # standard error of the AUC, not a variance
  auc <- as.numeric(r$AUC)
  data.table(time_label = labels, time_days = times, AUC = auc,
             AUC_low_95 = pmax(0, auc - 1.96 * se), AUC_high_95 = pmin(1, auc + 1.96 * se), SE = se)
}
cidx <- function(time, event, score) {                         # fixed score, higher score = shorter survival
  k <- concordance(Surv(time, event) ~ score, reverse = TRUE)
  c(C = unname(k$concordance), se = sqrt(k$var))
}

t10 <- auc_table(tr$OS.time, tr$OS, tr$RiskScore, c(365, 730, 1095), c("1-yr", "2-yr", "3-yr"))
t11 <- cbind(cohort = "GSE53625", auc_table(va$OS.time, va$OS, va$RiskScore, c(365, 1095, 1825), c("1-yr", "3-yr", "5-yr")))
ct <- cidx(tr$OS.time, tr$OS, tr$RiskScore); cv <- cidx(va$OS.time, va$OS, va$RiskScore)
cc <- data.table(cohort = c("TCGA-ESCA training", "GSE53625 validation"), n = c(nrow(tr), nrow(va)),
                 c_index = c(ct["C"], cv["C"]),
                 c_low_95 = c(ct["C"] - 1.96 * ct["se"], cv["C"] - 1.96 * cv["se"]),
                 c_high_95 = c(ct["C"] + 1.96 * ct["se"], cv["C"] + 1.96 * cv["se"]),
                 c_se = c(ct["se"], cv["se"]))
fwrite(t10, paste0(TBL, "Tab10_timeROC_AUC.csv"))
fwrite(t11, paste0(TBL, "Tab11_GSE53625_timeROC_AUC.csv"))
fwrite(cc, paste0(TBL, "Tab11_GSE53625_concordance_index.csv"))
options(digits = 6)
print(t10); print(t11); print(cc)

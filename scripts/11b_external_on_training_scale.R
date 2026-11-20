###############################################################################
## Stage 11b (2026-10-19): external evaluation in GSE53625 with the score on the training scale
##
## Why: Stage 11 computed the GSE53625 score as sum(beta * z), z = gene-wise z-scores within GSE53625, and split it at
## the TCGA median of the log2(norm_count + 1)-scale score (-0.0026). Two scale mismatches: (1) the cutoff belongs to
## the log2 scale, not to the z scale; (2) the TCGA coefficients were fitted on genes with different spreads
## (TCGA SD 1.7-3.4), so applying them to unit-variance z-scores changes the relative weight of the seven genes.
## Fix: place each GSE53625 gene on its TCGA distribution,
##   x' = mean_TCGA + SD_TCGA * z_GSE53625,
## so the TCGA coefficients and the TCGA median cutoff apply unchanged. Everything else is as in Stage 11:
## same 179 tumours with OS, same probe-to-gene averaging (AnnoProbe GPL18109 annotation, read from the cached file,
## no download), same timeROC times (1, 3, 5 years), same interval methods as Table 4 and Figure 8b
## (C-index: survival::concordance, higher score = shorter survival, C +/- 1.96 SE; AUC: timeROC iid, AUC +/- 1.96 SD).
## Inputs  tables/Tab10_LASSO_selected_genes.csv, tables/Tab10_risk_scores_per_sample.csv, data/expr_ESCA_tumor.csv,
##         data/GSE53625_series_matrix.txt.gz, cache/GPL18109_pipe.rda, tables/Tab11_GSE53625_RS_per_sample.csv
## Outputs tables/Tab11b_GSE53625_RS_per_sample_training_scale.csv, tables/Tab11b_GSE53625_metrics_training_scale.csv,
##         logs/11b_external_on_training_scale.log
## The manuscript is not changed by this script.
###############################################################################
# run from the repository root (see README.md)
log_con <- file("logs/11b_external_on_training_scale.log", open = "wt")
sink(log_con, split = TRUE)
suppressPackageStartupMessages({
  library(data.table); library(survival); library(timeROC); library(GEOquery); library(Biobase)
})
cat(sprintf("=== Stage 11b (%s): R %s, survival %s, timeROC %s, GEOquery %s ===\n", Sys.time(),
            paste(R.Version()[c("major", "minor")], collapse = "."), packageVersion("survival"),
            packageVersion("timeROC"), packageVersion("GEOquery")))

# ---------------------------------------------------------------- training model and scale
lasso <- fread("tables/Tab10_LASSO_selected_genes.csv")
rs_t <- fread("tables/Tab10_risk_scores_per_sample.csv")
ex <- fread("data/expr_ESCA_tumor.csv"); setnames(ex, 1, "gene")
Xt <- as.matrix(ex[match(lasso$gene, gene), rs_t$sampleID, with = FALSE]); rownames(Xt) <- lasso$gene
stopifnot(max(abs(as.numeric(t(Xt) %*% lasso$coefficient) - rs_t$RiskScore)) < 1e-10)   # Tab10 score reproduced
mu_t <- rowMeans(Xt); sd_t <- apply(Xt, 1, sd)
cut_t <- median(rs_t$RiskScore)
cat(sprintf("TCGA: %d tumours, median score (cutoff) %.4f\n", ncol(Xt), cut_t))
print(data.table(gene = lasso$gene, beta = lasso$coefficient, mean_TCGA = round(mu_t, 3), SD_TCGA = round(sd_t, 3)))

# ---------------------------------------------------------------- GSE53625 as in Stage 11
pub <- fread("tables/Tab11_GSE53625_RS_per_sample.csv")
eset <- getGEO(filename = "data/GSE53625_series_matrix.txt.gz", AnnotGPL = FALSE, destdir = "cache",
               GSEMatrix = TRUE, getGPL = FALSE)
eset <- if (inherits(eset, "list")) eset[[1]] else eset
expr_raw <- exprs(eset)
e <- new.env(); load("cache/GPL18109_pipe.rda", envir = e); ids <- get("GPL18109_pipe", e)
ids <- ids[!is.na(ids$symbol) & ids$symbol != "", ]
samples <- pub$sampleID
Xg <- t(sapply(lasso$gene, function(g) {
  p <- intersect(as.character(ids$probe_id[ids$symbol == g]), rownames(expr_raw))
  if (length(p) == 1) expr_raw[p, samples] else colMeans(expr_raw[p, samples, drop = FALSE], na.rm = TRUE)
}))
Zg <- t(scale(t(Xg)))
rs_z <- as.numeric(t(Zg) %*% lasso$coefficient)
d0 <- max(abs(rs_z - pub$RiskScore))
cat(sprintf("\nGSE53625: %d tumours, %d deaths; Stage 11 score rebuilt, max |difference| %.1e\n", length(samples),
            sum(pub$OS), d0))
stopifnot(d0 < 1e-10)

# ---------------------------------------------------------------- score on the training scale
Xmap <- mu_t + sd_t * Zg
rs_new <- as.numeric(t(Xmap) %*% lasso$coefficient)
grp_new <- factor(ifelse(rs_new <= cut_t, "Low", "High"), levels = c("Low", "High"))
grp_old <- factor(ifelse(rs_z <= cut_t, "Low", "High"), levels = c("Low", "High"))
stopifnot(all(grp_old == pub$RiskGroup))
cat(sprintf("new score: median %.4f (TCGA %.4f); Spearman with the Stage 11 score %.3f; %d of %d patients change group\n",
            median(rs_new), cut_t, cor(rs_new, rs_z, method = "spearman"), sum(grp_new != grp_old), length(grp_new)))

metrics <- function(score, grp, label) {
  k <- concordance(Surv(pub$OS.time, pub$OS) ~ score, reverse = TRUE); se <- sqrt(k$var)
  r <- timeROC(T = pub$OS.time, delta = pub$OS, marker = score, cause = 1, times = c(365, 1095, 1825), iid = TRUE)
  sdv <- r$inference$vect_sd_1
  s <- summary(coxph(Surv(OS.time, OS) ~ grp, data = data.table(pub, grp = grp)))
  lr <- survdiff(Surv(pub$OS.time, pub$OS) ~ grp)
  km <- summary(survfit(Surv(pub$OS.time, pub$OS) ~ grp))$table
  data.table(score = label, n = length(score), deaths = sum(pub$OS),
             n_low = sum(grp == "Low"), n_high = sum(grp == "High"),
             deaths_low = sum(pub$OS[grp == "Low"]), deaths_high = sum(pub$OS[grp == "High"]),
             median_OS_low_days = km["grp=Low", "median"], median_OS_high_days = km["grp=High", "median"],
             HR_high_vs_low = s$conf.int[1, 1], HR_low95 = s$conf.int[1, 3], HR_high95 = s$conf.int[1, 4],
             logrank_p = pchisq(lr$chisq, 1, lower.tail = FALSE),
             C = k$concordance, C_low95 = k$concordance - 1.96 * se, C_high95 = k$concordance + 1.96 * se,
             AUC_1y = r$AUC[1], AUC_1y_low = r$AUC[1] - 1.96 * sdv[1], AUC_1y_high = r$AUC[1] + 1.96 * sdv[1],
             AUC_3y = r$AUC[2], AUC_3y_low = r$AUC[2] - 1.96 * sdv[2], AUC_3y_high = r$AUC[2] + 1.96 * sdv[2],
             AUC_5y = r$AUC[3], AUC_5y_low = r$AUC[3] - 1.96 * sdv[3], AUC_5y_high = r$AUC[3] + 1.96 * sdv[3])
}
m <- rbind(metrics(rs_z, grp_old, "Stage 11: beta * z, TCGA log2 median"),
           metrics(rs_new, grp_new, "Stage 11b: training scale, TCGA log2 median"))
fwrite(m, "tables/Tab11b_GSE53625_metrics_training_scale.csv")
fwrite(data.table(sampleID = pub$sampleID, OS = pub$OS, OS.time = pub$OS.time, OS_months = pub$OS_months,
                  RiskScore = rs_new, RiskGroup = as.character(grp_new), cutoff_used = cut_t,
                  cutoff_src = "TCGA_training_median", RiskScore_stage11 = rs_z, RiskGroup_stage11 = as.character(grp_old)),
       "tables/Tab11b_GSE53625_RS_per_sample_training_scale.csv")

cat("\n")
for (i in 1:2) with(m[i], cat(sprintf(
  "%s\n  High/Low %d/%d (deaths %d/%d) | HR %.2f (%.2f-%.2f), log-rank p = %.3f | C = %.3f (%.3f-%.3f)\n  AUC 1 y %.3f (%.3f-%.3f), 3 y %.3f (%.3f-%.3f), 5 y %.3f (%.3f-%.3f)\n",
  score, n_high, n_low, deaths_high, deaths_low, HR_high_vs_low, HR_low95, HR_high95, logrank_p, C, C_low95, C_high95,
  AUC_1y, AUC_1y_low, AUC_1y_high, AUC_3y, AUC_3y_low, AUC_3y_high, AUC_5y, AUC_5y_low, AUC_5y_high)))
cat(sprintf("\n=== end (%s) ===\n", Sys.time()))
sink(); close(log_con)

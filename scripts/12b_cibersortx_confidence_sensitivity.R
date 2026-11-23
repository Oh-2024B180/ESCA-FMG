###############################################################################
## Stage 12b: risk score - CIBERSORTx correlations by CIBERSORTx confidence class
##
## Inputs  data/cibersortx_LM22_TCGA_ESCA_Bmode.csv (CIBERSORTx B-mode output, column "P-value"),
##         tables/Tab8B_CIBERSORT_results.csv (fractions per sample), tables/Tab10_risk_scores_per_sample.csv,
##         tables/Tab12_RS_immune_spearman.csv
##
## Output  printed table only
###############################################################################
suppressPackageStartupMessages(library(data.table))
cx <- fread("data/cibersortx_LM22_TCGA_ESCA_Bmode.csv")
fr <- fread("tables/Tab8B_CIBERSORT_results.csv")
rs <- fread("tables/Tab10_risk_scores_per_sample.csv")
sp <- fread("tables/Tab12_RS_immune_spearman.csv")[feature_category == "CIBERSORT"]
wx <- fread("tables/Tab12_RS_immune_wilcoxon.csv")[feature_category == "CIBERSORT"]
cx[, conf := fifelse(`P-value` < 0.05, "confident", fifelse(`P-value` < 0.10, "borderline", "low"))]
print(cx[, .N, by = conf])
cells <- setdiff(names(fr), c("sampleID", "cluster", "histology", "P-value", "Correlation", "RMSE"))
d <- merge(fr, rs[, .(sampleID, RiskScore, RiskGroup)], by = "sampleID")
run <- function(dd) rbindlist(lapply(cells, function(f) {
  ct <- suppressWarnings(cor.test(dd$RiskScore, dd[[f]], method = "spearman", exact = FALSE))
  wt <- tryCatch(wilcox.test(dd[RiskGroup == "Low"][[f]], dd[RiskGroup == "High"][[f]], exact = FALSE)$p.value,
                 error = function(e) NA_real_)
  data.table(feature = f, rho = unname(ct$estimate), p_spearman = ct$p.value, p_wilcoxon = wt)
}))
keep <- cx[conf == "confident", Mixture]
conf <- run(d[sampleID %in% keep])
all <- run(d)
chk <- merge(conf, merge(sp[, .(feature, rho_tab = rho, ps_tab = p_value)], wx[, .(feature, pw_tab = p_value)],
                         by = "feature"), by = "feature")
stopifnot(max(abs(chk$rho - chk$rho_tab), abs(chk$p_spearman - chk$ps_tab),
              abs(chk$p_wilcoxon - chk$pw_tab), na.rm = TRUE) < 1e-12)          # confident run = stage 12
cat(sprintf("samples: confident %d (stage 12), all %d (left out: %s)\n", sum(d$sampleID %in% keep), nrow(d),
            paste(setdiff(d$sampleID, keep), collapse = ", ")))
for (v in c("p_spearman", "p_wilcoxon")) {
  cat(sprintf("%s < 0.05, confident samples: %s\n", v, paste(conf[get(v) < 0.05, feature], collapse = ", ")))
  cat(sprintf("%s < 0.05, all samples      : %s\n", v, paste(all[get(v) < 0.05, feature], collapse = ", ")))
}
print(merge(conf, all, by = "feature", suffixes = c("_confident", "_all"))[order(p_spearman_confident)][1:6])

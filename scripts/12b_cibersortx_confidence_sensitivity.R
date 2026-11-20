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
cx[, conf := fifelse(`P-value` < 0.05, "confident", fifelse(`P-value` < 0.10, "borderline", "low"))]
print(cx[, .N, by = conf])
cells <- setdiff(names(fr), c("sampleID", "cluster", "histology"))
d <- merge(fr, rs[, .(sampleID, RiskScore)], by = "sampleID")
run <- function(dd) rbindlist(lapply(cells, function(f) {
  ct <- suppressWarnings(cor.test(dd$RiskScore, dd[[f]], method = "spearman", exact = FALSE))
  data.table(feature = f, rho = unname(ct$estimate), p = ct$p.value)
}))
all <- run(d)
chk <- merge(all, sp[, .(feature, rho_tab = rho, p_tab = p_value)], by = "feature")
stopifnot(max(abs(chk$rho - chk$rho_tab), na.rm = TRUE) < 1e-12)
keep <- cx[conf == "confident", Mixture]
conf <- run(d[sampleID %in% keep])
cat(sprintf("samples: all %d, confident %d (dropped: %s)\n", nrow(d), sum(d$sampleID %in% keep),
            paste(setdiff(d$sampleID, keep), collapse = ", ")))
cat("p < 0.05, all samples      :", paste(all[p < 0.05, feature], collapse = ", "), "\n")
cat("p < 0.05, confident samples:", paste(conf[p < 0.05, feature], collapse = ", "), "\n")
print(merge(all, conf, by = "feature", suffixes = c("_all", "_confident"))[order(p_all)][1:6])

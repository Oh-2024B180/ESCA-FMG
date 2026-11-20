###############################################################################
## Stage 12c: discrimination of the risk score within each histology
##
## Inputs  tables/Tab10_risk_scores_per_sample.csv (RiskScore, RiskGroup, OS, OS.time),
##         tables/Tab7_cluster_assignment.csv (histology)
##
## Output  tables/Tab12c_within_histology_discrimination.csv
###############################################################################
suppressPackageStartupMessages({ library(survival); library(data.table) })
rs <- fread("tables/Tab10_risk_scores_per_sample.csv")
cl <- fread("tables/Tab7_cluster_assignment.csv")
d <- merge(rs[, .(sampleID, RiskScore, RiskGroup, OS, OS.time)], cl[, .(sampleID, histology)], by = "sampleID")
d[, H := fifelse(grepl("SQUAMOUS|ESCC", toupper(histology)), "ESCC",
         fifelse(grepl("ADENO|EAC", toupper(histology)), "EAC", NA_character_))]
stopifnot(!anyNA(d$H), nrow(d) == 181)

res <- rbindlist(lapply(c("ESCC", "EAC"), function(h) {
  x <- d[H == h]
  x[, high := as.integer(RiskGroup == "High")]
  k <- concordance(Surv(OS.time, OS) ~ RiskScore, data = x, reverse = TRUE)
  se <- sqrt(k$var)
  u <- summary(coxph(Surv(OS.time, OS) ~ RiskScore, data = x))
  g <- summary(coxph(Surv(OS.time, OS) ~ high, data = x))
  lr <- survdiff(Surv(OS.time, OS) ~ high, data = x)
  data.table(histology = h, n = nrow(x), events = sum(x$OS),
             C = k$concordance, C_low95 = k$concordance - 1.96 * se, C_high95 = k$concordance + 1.96 * se,
             HR_per_unit = u$conf.int[1], HR_per_unit_low95 = u$conf.int[3], HR_per_unit_high95 = u$conf.int[4],
             HR_per_unit_p = u$coefficients[1, 5],
             n_high = sum(x$high), n_low = sum(1 - x$high),
             HR_group = g$conf.int[1], HR_group_low95 = g$conf.int[3], HR_group_high95 = g$conf.int[4],
             logrank_p = pchisq(lr$chisq, df = 1, lower.tail = FALSE))
}))
fwrite(res, "tables/Tab12c_within_histology_discrimination.csv")
for (i in seq_len(nrow(res))) with(res[i], cat(sprintf(
  "%-4s n = %d (%d events) | C %.3f (%.3f-%.3f) | HR per unit %.2f (%.2f-%.2f), p = %.2g | high/low %d/%d | HR %.2f (%.2f-%.2f) | log-rank p = %.2g\n",
  histology, n, events, C, C_low95, C_high95, HR_per_unit, HR_per_unit_low95, HR_per_unit_high95, HR_per_unit_p,
  n_high, n_low, HR_group, HR_group_low95, HR_group_high95, logrank_p)))
cat("[OK] tables/Tab12c_within_histology_discrimination.csv\n")

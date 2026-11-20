###############################################################################
## Stage 14C after the oncoPredict check (2026-10-19)
##
## The GDSC2 training response (GDSC2_Res.rds) holds natural-log IC50 values. The oncoPredict vignette
## (vignettes/calcPhenotype.Rmd) exponentiates them before calcPhenotype(), which then power-transforms the response,
## "to not have both transformations"; stage 14C had passed the log values directly. The predictions used here are
## run B of scripts/14b_oncopredict_predictions.R (oncoPredict 1.2, R 4.4.1, 63 min; trainingPtype = exp(GDSC2_Res),
## everything else as in stage 14C: batchCorrect "eb", powerTransformPhenotype TRUE, removeLowVaryingGenes 0.4,
## minNumSamples 10, removeLowVaringGenesFrom "rawData"), saved as tables/Tab14Cb_oncoPredict_IC50_vignette.csv
## (calcPhenotype returns them on the IC50 scale). Run A of the same script reproduced the stage 14C table exactly.
## This script puts the predictions on the natural-log scale used before, repeats the stage 14C tests unchanged and
## redraws Figure 11c with the Fig10D block of figures/remake/v3_fig8_10_remake.R (unchanged).
## Outputs tables/Tab14C_oncoPredict_IC50_full.csv, tables/Tab14C_oncoPredict_drug_test.csv,
##         tables/Tab14C_oncoPredict_top30_drugs.csv, figures/remake/Fig10D_oncoPredict_top10_boxplot_v3.{pdf,png},
##         logs/14c_oncopredict_vignette.log
###############################################################################
# run from the repository root (see README.md)
TBLDIR <- "tables"; OUT_DIR <- "figures/remake"; P_CUT <- 0.05; TOP_N_DRUG <- 30
log_con <- file("logs/14c_oncopredict_vignette.log", open = "wt")
sink(log_con, split = TRUE)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ragg) })
cat(sprintf("=== Stage 14C, vignette predictions (%s): R %s ===\n", Sys.time(), paste(R.Version()[c("major", "minor")], collapse = ".")))

rs_dt <- fread("tables/Tab10_risk_scores_per_sample.csv")
rs_dt[, RiskGroup := factor(RiskGroup, levels = c("Low", "High"))]
pred <- fread("tables/Tab14Cb_oncoPredict_IC50_vignette.csv")
drug_cols <- setdiff(colnames(pred), "sampleID")
stopifnot(nrow(pred) == 181, all(pred$sampleID %in% rs_dt$sampleID), all(as.matrix(pred[, ..drug_cols]) > 0, na.rm = TRUE))
ic50_dt <- copy(pred)
ic50_dt[, (drug_cols) := lapply(.SD, log), .SDcols = drug_cols]                 # natural log, as in the published table
cat(sprintf("predictions: %d samples x %d drugs, missing values %d\n", nrow(ic50_dt), length(drug_cols),
            sum(is.na(as.matrix(ic50_dt[, ..drug_cols])))))
fwrite(ic50_dt, file.path(TBLDIR, "Tab14C_oncoPredict_IC50_full.csv"))

# ---------------------------------------------------------------- stage 14C tests (unchanged)
wilcox_lh <- function(values, group) {
  ok <- !is.na(values) & is.finite(values) & !is.na(group)
  values <- values[ok]; group <- group[ok]
  wt <- wilcox.test(values ~ group, exact = FALSE)
  list(W = unname(wt$statistic), p = wt$p.value,
       med_low = median(values[group == "Low"], na.rm = TRUE), med_high = median(values[group == "High"], na.rm = TRUE),
       n_low = sum(group == "Low"), n_high = sum(group == "High"))
}
ic50_join <- merge(rs_dt[, .(sampleID, RiskGroup)], ic50_dt, by = "sampleID")
test_c <- rbindlist(lapply(drug_cols, function(dr) {
  res <- wilcox_lh(as.numeric(ic50_join[[dr]]), ic50_join$RiskGroup)
  data.table(drug = dr, n_low = res$n_low, n_high = res$n_high, median_low = res$med_low, median_high = res$med_high,
             W = res$W, p_value = res$p)
}))
test_c[, p_adj_BH := p.adjust(p_value, method = "BH")]
test_c[, log_FC_HighMinusLow := median_high - median_low]
test_c[, direction := fifelse(log_FC_HighMinusLow < 0, "High_more_sensitive",
                        fifelse(log_FC_HighMinusLow > 0, "High_more_resistant", "tie"))]
setorder(test_c, p_value)
fwrite(test_c, file.path(TBLDIR, "Tab14C_oncoPredict_drug_test.csv"))
fwrite(head(test_c, TOP_N_DRUG), file.path(TBLDIR, "Tab14C_oncoPredict_top30_drugs.csv"))
cat(sprintf("raw p < %.2f: %d / %d drugs (%.1f%%); smallest BH-adjusted p %.3f\n", P_CUT, sum(test_c$p_value < P_CUT),
            nrow(test_c), 100 * mean(test_c$p_value < P_CUT), min(test_c$p_adj_BH)))
print(test_c[1:10, .(drug, p_value = signif(p_value, 3), median_low = round(median_low, 3), median_high = round(median_high, 3),
                     direction)])

# ---------------------------------------------------------------- Figure 11c (v3_fig8_10_remake.R, Fig10D block)
COL_LOW <- "#377EB8"; COL_HIGH <- "#E41A1C"
THEME_BASE <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black", size = 11),
        axis.title = element_text(size = 12), strip.text = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))
sig_from_p <- function(p) {
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.05] <- "*"; out[!is.na(p) & p < 0.01] <- "**"; out[!is.na(p) & p < 0.001] <- "***"
  out
}
save_figure <- function(plot_obj, stem, w_in = 8, h_in = 8) {
  ggsave(file.path(OUT_DIR, paste0(stem, ".pdf")), plot_obj, width = w_in, height = h_in, device = cairo_pdf)
  agg_png(file.path(OUT_DIR, paste0(stem, ".png")), width = w_in, height = h_in, units = "in", res = 600)
  print(plot_obj)
  invisible(dev.off())
  cat(sprintf("  saved: %s (%.1f x %.1f in)\n", stem, w_in, h_in))
}
drug_test <- fread(file.path(TBLDIR, "Tab14C_oncoPredict_top30_drugs.csv"))
setorder(drug_test, p_value)
drug_top10 <- drug_test[seq_len(min(10, .N))]
ic50_plot <- merge(fread(file.path(TBLDIR, "Tab14C_oncoPredict_IC50_full.csv")), rs_dt[, .(sampleID, RiskGroup)], by = "sampleID")
drug_long <- rbindlist(lapply(seq_len(nrow(drug_top10)), function(i) {
  drug <- drug_top10$drug[i]
  if (!drug %in% colnames(ic50_plot)) return(NULL)
  data.table(sampleID = ic50_plot$sampleID, RiskGroup = ic50_plot$RiskGroup, drug = drug, value = ic50_plot[[drug]])
}))
drug_long <- drug_long[!is.na(value)]
drug_long[, drug := factor(drug, levels = drug_top10$drug)]
sig_10D <- drug_top10[, .(drug, raw_p = p_value)]
sig_10D[, sig := sig_from_p(raw_p)]
sig_10D[, drug := factor(drug, levels = drug_top10$drug)]
y_top10 <- drug_long[, .(ymax = max(value, na.rm = TRUE), ymin = min(value, na.rm = TRUE)), by = drug]
sig_10D <- merge(sig_10D, y_top10, by = "drug")
sig_10D[, ylabel := ymin + (ymax - ymin) * 0.92]
fig10D <- ggplot(drug_long, aes(x = RiskGroup, y = value, fill = RiskGroup)) +
  geom_boxplot(width = 0.55, outlier.size = 0.6, alpha = 0.85, colour = "black") +
  geom_text(data = sig_10D, aes(x = 1.5, y = ylabel, label = sig), inherit.aes = FALSE, size = 6, fontface = "bold") +
  facet_wrap(~ drug, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c(Low = COL_LOW, High = COL_HIGH), name = "Risk group") +
  labs(x = NULL, y = "Predicted ln(IC50)") +
  THEME_BASE +
  theme(legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.6),
        strip.text = element_text(size = 10.5, face = "bold"),
        panel.spacing = unit(0.6, "lines"))
save_figure(fig10D, "Fig10D_oncoPredict_top10_boxplot_v3", w_in = 12, h_in = 6)
cat(sprintf("=== end (%s) ===\n", Sys.time()))
sink(); close(log_con)

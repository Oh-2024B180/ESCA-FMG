###############################################################################
## Stage 11c: Figure 8 panels for the external evaluation on the training scale (stage 11b)
##
## Redraws the three panels of manuscript Figure 8 with the stage 11b score, using the code that drew them before:
##   (a) Fig7A KM + risk table   : stage 11 code (scripts/11_external_validation_GSE53625.R), 8 x 7 in
##   (b) Fig7B timeROC           : scripts/helpers/11_fix_Fig7B_legend.R, 7.2 x 6.6 in
##   (c) Fig7C score/status/heat : scripts/helpers/11_fix_Fig7CDE_cosmetic.R, 11 x 5.46 in
## The text that was later corrected in the assembled image by redrawing pixels is drawn here directly: en dashes in
## the HR and AUC intervals, the AUC intervals as AUC +/- 1.96 SD, italic gene names and the
## arrow in the heatmap axis title. Also redraws Fig7E (C-index, TCGA vs GSE53625; not used in the article) with
## the stage 11b value, using the code of scripts/helpers/11_fix_Fig7E_cindex.R. Fig7D is not used anywhere and is left.
## Inputs  tables/Tab11b_GSE53625_RS_per_sample_training_scale.csv, tables/Tab11b_GSE53625_metrics_training_scale.csv,
##         tables/Tab10_LASSO_selected_genes.csv, tables/Tab11_GSE53625_concordance_index.csv (TCGA row of Fig7E),
##         data/GSE53625_series_matrix.txt.gz, cache/GPL18109_pipe.rda (heatmap z-scores, as in stage 11)
## Outputs figures/Fig7A_*, Fig7B_*, Fig7C_*, Fig7E_* (.pdf/.png),
##         tables/Tab11b_GSE53625_concordance_index_training_scale.csv, logs/11c_external_figures.log
###############################################################################
# run from the repository root (see README.md)
FIGDIR <- "figures"
log_con <- file("logs/11c_external_figures.log", open = "wt")
sink(log_con, split = TRUE)
suppressPackageStartupMessages({
  library(data.table); library(survival); library(timeROC); library(ggplot2); library(survminer)
  library(cowplot); library(ragg); library(GEOquery); library(Biobase)
})
cat(sprintf("=== Stage 11c (%s): R %s, ggplot2 %s, survminer %s, cowplot %s ===\n", Sys.time(),
            paste(R.Version()[c("major", "minor")], collapse = "."), packageVersion("ggplot2"),
            packageVersion("survminer"), packageVersion("cowplot")))
RISK_COL <- c(Low = "#377EB8", High = "#E41A1C")
EN <- "\u2013"

rs <- fread("tables/Tab11b_GSE53625_RS_per_sample_training_scale.csv")
m <- fread("tables/Tab11b_GSE53625_metrics_training_scale.csv")[2]
stopifnot(grepl("^Stage 11b", m$score))
rs[, RiskGroup := factor(RiskGroup, levels = c("Low", "High"))]
cutoff_used <- rs$cutoff_used[1]
cat(sprintf("GSE53625: %d tumours, High/Low %d/%d, cutoff %.4f\n", nrow(rs), sum(rs$RiskGroup == "High"),
            sum(rs$RiskGroup == "Low"), cutoff_used))

# ---------------------------------------------------------------- (a) KM, stage 11 code
fit_km <- survfit(Surv(OS.time, OS) ~ RiskGroup, data = rs)
s <- summary(coxph(Surv(OS.time, OS) ~ RiskGroup, data = rs))
HR <- s$conf.int[1, 1]; HR_low <- s$conf.int[1, 3]; HR_high <- s$conf.int[1, 4]
lr_p <- pchisq(survdiff(Surv(OS.time, OS) ~ RiskGroup, data = rs)$chisq, 1, lower.tail = FALSE)
stopifnot(abs(HR - m$HR_high_vs_low) < 1e-12, abs(lr_p - m$logrank_p) < 1e-12)
ss <- ggsurvplot(fit_km, data = rs, pval = FALSE, risk.table = TRUE, conf.int = FALSE,
                 palette = c(unname(RISK_COL["Low"]), unname(RISK_COL["High"])),
                 legend.labs = c("Low risk", "High risk"), legend.title = "Risk group",
                 xlab = "Time (days)", ylab = "Overall survival", title = NULL,
                 ggtheme = theme_bw(base_size = 11), risk.table.fontsize = 3.5)
hr_lbl_7A <- sprintf("HR = %.2f (95%% CI %.2f%s%.2f)\nlog-rank p = %.3g", HR, HR_low, EN, HR_high, lr_p)
ss$plot <- ss$plot +
  annotate("label", x = Inf, y = Inf, label = hr_lbl_7A, hjust = 1.05, vjust = 1.2, size = 3.6, fontface = "plain",
           label.size = 0.3, fill = "white", label.r = unit(0, "pt"))
combined_7A <- cowplot::plot_grid(ss$plot + theme(legend.position = "top"),
                                  ss$table + theme(plot.title = element_blank()),
                                  ncol = 1, align = "v", axis = "lr", rel_heights = c(3, 1))
ggsave(file.path(FIGDIR, "Fig7A_GSE53625_KM_by_risk_group.pdf"), plot = combined_7A, width = 8, height = 7, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig7A_GSE53625_KM_by_risk_group.png"), plot = combined_7A, width = 8, height = 7, dpi = 600,
       device = ragg::agg_png)
cat("  [a] ", gsub("\n", " | ", hr_lbl_7A), "\n")

# ---------------------------------------------------------------- (b) timeROC, 11_fix_Fig7B_legend.R
roc_times_used <- c(365, 1095, 1825); roc_labels <- c("1-yr", "3-yr", "5-yr")
roc_res <- timeROC(T = rs$OS.time, delta = rs$OS, marker = rs$RiskScore, cause = 1, times = roc_times_used, iid = TRUE)
auc <- c(m$AUC_1y, m$AUC_3y, m$AUC_5y)
lo <- c(m$AUC_1y_low, m$AUC_3y_low, m$AUC_5y_low); hi <- c(m$AUC_1y_high, m$AUC_3y_high, m$AUC_5y_high)
stopifnot(max(abs(as.numeric(roc_res$AUC) - auc)) < 1e-12)
roc_df <- do.call(rbind, lapply(seq_along(roc_times_used), function(k) {
  data.frame(FPR = roc_res$FP[, k], TPR = roc_res$TP[, k],
             time_label = sprintf("%s (AUC=%.3f, 95%% CI %.3f%s%.3f)", roc_labels[k], auc[k], lo[k], EN, hi[k]))
}))
palette_used <- c("#377EB8", "#4DAF4A", "#E41A1C"); names(palette_used) <- unique(roc_df$time_label)
p7B <- ggplot(roc_df, aes(x = FPR, y = TPR, color = time_label)) +
  geom_line(linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = palette_used, name = NULL) +
  coord_equal() +
  labs(x = "False positive rate", y = "True positive rate") +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.62, 0.16), legend.justification = c(0.5, 0.5),
        legend.background = element_rect(fill = "white", color = "grey70"),
        legend.key.height = unit(0.95, "lines"), legend.text = element_text(size = 9),
        legend.margin = margin(4, 5, 4, 5))
ggsave(file.path(FIGDIR, "Fig7B_GSE53625_timeROC.pdf"), plot = p7B, width = 7.2, height = 6.6, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig7B_GSE53625_timeROC.png"), plot = p7B, width = 7.2, height = 6.6, dpi = 600,
       device = ragg::agg_png)
cat("  [b] ", paste(unique(roc_df$time_label), collapse = " | "), "\n")

# ---------------------------------------------------------------- (c) 3-panel, 11_fix_Fig7CDE_cosmetic.R
lasso <- fread("tables/Tab10_LASSO_selected_genes.csv"); mapped_genes <- lasso$gene
eset <- getGEO(filename = "data/GSE53625_series_matrix.txt.gz", AnnotGPL = FALSE, destdir = "cache", GSEMatrix = TRUE,
               getGPL = FALSE)
eset <- if (inherits(eset, "list")) eset[[1]] else eset
expr_raw <- exprs(eset)
e <- new.env(); load("cache/GPL18109_pipe.rda", envir = e); ids <- get("GPL18109_pipe", e)
ids <- ids[!is.na(ids$symbol) & ids$symbol != "", ]
expr_lasso <- t(sapply(mapped_genes, function(g) {
  p <- intersect(as.character(ids$probe_id[ids$symbol == g]), rownames(expr_raw))
  if (length(p) == 1) expr_raw[p, rs$sampleID] else colMeans(expr_raw[p, rs$sampleID, drop = FALSE], na.rm = TRUE)
}))
expr_lasso_z <- t(scale(t(expr_lasso), center = TRUE, scale = TRUE))
ord <- order(rs$RiskScore)
rs_sorted <- rs$RiskScore[ord]; samples_sort <- rs$sampleID[ord]; rg_sorted <- as.character(rs$RiskGroup)[ord]
os_event_s <- rs$OS[ord]; os_time_s <- rs$OS.time[ord]; N <- length(ord)
p1 <- ggplot(data.frame(rank = 1:N, RS = rs_sorted, RiskGroup = factor(rg_sorted, levels = c("Low", "High"))),
             aes(x = rank, y = RS, color = RiskGroup)) +
  geom_point(size = 0.8) +
  geom_hline(yintercept = cutoff_used, linetype = "dashed", color = "grey40") +
  annotate("text", x = N * 0.02, y = cutoff_used, label = sprintf("cutoff = %.4f (TCGA training)", cutoff_used),
           hjust = 0, vjust = 1.4, size = 3, color = "grey25", fontface = "italic") +
  scale_color_manual(values = c(Low = unname(RISK_COL["Low"]), High = unname(RISK_COL["High"])), name = "Risk group") +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.12))) +
  labs(x = NULL, y = "Risk Score") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = margin(8, 8, 4, 8))
p2 <- ggplot(data.frame(rank = 1:N, OS_time = os_time_s,
                        Event = factor(ifelse(os_event_s == 1, "Death", "Censored"), levels = c("Censored", "Death"))),
             aes(x = rank, y = OS_time, color = Event)) +
  geom_point(size = 0.9, alpha = 0.8) +
  scale_color_manual(values = c(Censored = "grey30", Death = "#E41A1C"), name = "OS event") +
  labs(x = NULL, y = "Time (days)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank())
heat_dt <- as.data.table(expr_lasso_z[, ord, drop = FALSE], keep.rownames = "gene")
heat_long <- melt(heat_dt, id.vars = "gene", variable.name = "sample", value.name = "z")
heat_long[, sample := factor(sample, levels = samples_sort)]
heat_long[, gene := factor(gene, levels = rev(mapped_genes))]
p3 <- ggplot(heat_long, aes(x = sample, y = gene, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#0571B0", mid = "white", high = "#CA0020", midpoint = 0, name = "z-score",
                       limits = c(-3, 3), oob = scales::squish) +
  labs(x = sprintf("GSE53625 samples sorted by Risk Score (low \u2192 high, n=%d)", N), y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8, face = "italic"), panel.grid = element_blank(),
        legend.position = "right")
combined_7C <- cowplot::plot_grid(p1, p2, p3, ncol = 1, align = "v", axis = "lr",
                                  rel_heights = c(1.25, 1.1, max(1.5, 0.18 * length(mapped_genes) + 0.6)))
H_7C <- 3.2 + 0.18 * length(mapped_genes) + 1
ggsave(file.path(FIGDIR, "Fig7C_GSE53625_RS_distribution_3panel.pdf"), plot = combined_7C, width = 11, height = H_7C,
       device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(FIGDIR, "Fig7C_GSE53625_RS_distribution_3panel.png"), plot = combined_7C, width = 11, height = H_7C,
       dpi = 600, device = ragg::agg_png, limitsize = FALSE)
cat(sprintf("  [c] score range %.3f to %.3f, cutoff %.4f, %.2f in high\n", min(rs_sorted), max(rs_sorted), cutoff_used, H_7C))

# ---------------------------------------------------------------- Fig7E (guide only): code of helpers/11_fix_Fig7E_cindex.R
old_c <- fread("tables/Tab11_GSE53625_concordance_index.csv")
cidx_dt <- rbind(old_c[cohort == "TCGA-ESCA training"],
                 data.table(cohort = "GSE53625 validation", n = nrow(rs), c_index = m$C, c_low_95 = m$C_low95,
                            c_high_95 = m$C_high95, c_se = (m$C_high95 - m$C) / 1.96))
fwrite(cidx_dt, "tables/Tab11b_GSE53625_concordance_index_training_scale.csv")
p7E <- ggplot(cidx_dt, aes(x = cohort, y = c_index, fill = cohort)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = c_low_95, ymax = c_high_95), width = 0.2, color = "black", linewidth = 0.4) +
  geom_text(aes(y = c_high_95, label = sprintf("%.3f", c_index)), vjust = -0.8, size = 4, fontface = "bold") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0.04))) +
  scale_fill_manual(values = c("TCGA-ESCA training" = "#377EB8", "GSE53625 validation" = "#E41A1C"), guide = "none") +
  labs(x = NULL, y = "Concordance index") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(face = "bold"), plot.margin = margin(8, 8, 8, 8))
ggsave(file.path(FIGDIR, "Fig7E_TCGA_vs_GSE53625_cindex_comparison.pdf"), plot = p7E, width = 7, height = 5)
ggsave(file.path(FIGDIR, "Fig7E_TCGA_vs_GSE53625_cindex_comparison.png"), plot = p7E, width = 7, height = 5, dpi = 600,
       device = ragg::agg_png)
print(cidx_dt)
cat(sprintf("=== end (%s) ===\n", Sys.time()))
sink(); close(log_con)

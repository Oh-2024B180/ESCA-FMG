###############################################################################
## the panels of Figures 3 and 8-11 whose text was below 7 pt at the printed size
## (audit 11-15: 3b, 8a, 8b, 9b, 10b, 11a), redrawn with the helpers of redraw_fig4_fig5_fig7.R.
##
##   Rscript redraw_fig3_fig8_to_11.R [3] [8] [9] [10] [11]        (R 4.4.1; default: all)
##
## Figures 8 and 9 are drawn whole: 8 takes the layout of the new Figure 7, and 9b needs more room than its old cell,
## so 9 is laid out again and 9a/9c (6-6.5 pt) are redrawn with it. For 3b, 10b and 11a only the panel is drawn, at
## the pixel size of its place in figures/assembled/fig_2, fig_9, fig_10.jpg; paste_panels.py puts them in.
## Sources (read-only), the code that drew the published panels:
##   8   scripts/11c_external_figures.R: Tab11b (stage 11b, training scale), GSE53625 series matrix + cache/GPL18109_pipe.rda
##   9   figures/remake/v3_fig8_10_remake.R: TCGA clinical matrix, Tab12_RS_clinical_test/immune_spearman/immune_wilcoxon,
##       Tab8A-C, expression
##   3b  scripts/05_TCGA_CNV_analysis.R: Tab5_TCGA_CNV_frequency.csv, gene order of Tab4
##   11a scripts/14a_tide_rerun.R: data/esca_tide_output_hgnc_centered.txt, Tab14A_TIDE_metrics_test.csv
##   10b nothing is recomputed: the data area of each published calibration panel (fig10b_data_*.png, cut from
##       figures/assembled/fig_9.jpg inside the frame; extents in fig10b_extents.csv) gets new 7-pt axes
## Stars from raw p, as in the published panels.
###############################################################################
which_figs <- "none"
source("scripts/figures/redraw_fig4_fig5_fig7.R")
suppressPackageStartupMessages(library(png))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) args <- c("3", "8", "9", "10", "11")
PXCM <- 4000 / 16.5                                  # px per cm in the assembled figures
save_panel <- function(p, stem, w_px, h_px) {       # a panel at the pixel size of its place in the assembled figure
  agg_png(file.path(OUT, paste0(stem, ".png")), width = w_px, height = h_px, units = "px", res = RES)
  print(p); invisible(dev.off())
  cat(sprintf("  %s: %d x %d px = %.2f x %.2f cm\n", stem, w_px, h_px, w_px / PXCM, h_px / PXCM))
}
CP <- tb("Tab8D_checkpoint_wilcoxon.csv")
IT_SYMBOLS <- CP[gene_display == gene_symbol, gene_display]   # checkpoint labels that are HGNC symbols (italic)

# risk score / survival status / heat map, as Figure 7e
risk3 <- function(rs, z, genes, hline, hline_lab, xlab, tag) {
  ord <- order(rs$RiskScore); N <- nrow(rs)
  r1 <- data.table(rank = seq_len(N), RS = rs$RiskScore[ord], grp = factor(rs$RiskGroup[ord], levels = c("Low", "High")),
                   time = rs$OS.time[ord], ev = factor(ifelse(rs$OS[ord] == 1, "Death", "Censored"), levels = c("Censored", "Death")))
  xs <- scale_x_continuous(expand = expansion(add = 1))
  bare_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
                  legend.justification = c(0, 0.5))
  e1 <- ggplot(r1, aes(rank, RS, colour = grp)) + geom_point(size = 0.55) +
    geom_hline(yintercept = hline, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    annotate("text", x = 3, y = hline, hjust = 0, vjust = -0.5, size = 7 / .pt, family = "Arial", colour = "grey30",
             label = hline_lab) +
    scale_colour_manual(values = COL_RISK, name = "Risk group") + xs + scale_y_continuous(labels = lab_num) +
    labs(y = "Risk score", tag = tag) + theme_fig + bare_x
  e2 <- ggplot(r1, aes(rank, time, colour = ev)) + geom_point(size = 0.55, alpha = 0.85) +
    scale_colour_manual(values = c(Censored = "grey35", Death = "#E41A1C"), name = "OS event") + xs +
    labs(y = "Time (days)") + theme_fig + bare_x
  hd <- data.table(rank = rep(seq_len(N), times = length(genes)), gene = factor(rep(genes, each = N), levels = rev(genes)),
                   z = as.vector(t(z[genes, ord, drop = FALSE])))
  e3 <- ggplot(hd, aes(rank, gene, fill = z)) + geom_raster() +
    scale_fill_gradient2(low = "#0571B0", mid = "white", high = "#CA0020", midpoint = 0, limits = c(-3, 3),
                         oob = scales::squish, name = "z-score") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0), labels = function(l) italic_labels(l, TRUE)) +
    labs(x = xlab, y = NULL) + theme_fig +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major = element_blank(),
          legend.justification = c(0, 0.5), legend.key.height = unit(3.5, "mm"))
  wrap_elements(full = e1 / e2 / e3 + plot_layout(heights = c(1, 1, 1.25)))
}
roc_plot <- function(T, delta, marker, times, labs, tag, inside = c(0.99, 0.02), just = c(1, 0)) {
  roc <- timeROC(T = T, delta = delta, marker = marker, cause = 1, times = times, iid = TRUE)
  rd <- rbindlist(lapply(seq_along(times), function(k) data.table(FPR = roc$FP[, k], TPR = roc$TP[, k], lab = labs[k])))
  p <- ggplot(rd, aes(FPR, TPR, colour = lab)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    geom_line(linewidth = 0.5) +
    scale_colour_manual(values = setNames(c("#377EB8", "#4DAF4A", "#E41A1C"), labs), breaks = labs, name = NULL) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "False positive rate", y = "True positive rate", tag = tag) + theme_fig +
    theme(legend.position = "inside", legend.position.inside = inside, legend.justification = just,
          legend.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.25),
          legend.key.size = unit(2.6, "mm"), legend.margin = margin(1.5, 2, 1.5, 2))
  list(plot = p, auc = as.numeric(roc$AUC))
}

# =====================================================================================================  Figure 8
if ("8" %in% args) {
  cat("Figure 8 (external cohort GSE53625, stage 11b scale)\n")
  suppressPackageStartupMessages({ library(GEOquery); library(Biobase) })
  rs8 <- tb("Tab11b_GSE53625_RS_per_sample_training_scale.csv")
  m8 <- tb("Tab11b_GSE53625_metrics_training_scale.csv")[2]
  stopifnot(grepl("^Stage 11b", m8$score))
  rs8[, RiskGroup := factor(RiskGroup, levels = c("Low", "High"))]
  s8 <- summary(coxph(Surv(OS.time, OS) ~ RiskGroup, data = rs8))
  lr8 <- pchisq(survdiff(Surv(OS.time, OS) ~ RiskGroup, data = rs8)$chisq, 1, lower.tail = FALSE)
  stopifnot(abs(s8$conf.int[1, 1] - m8$HR_high_vs_low) < 1e-12, abs(lr8 - m8$logrank_p) < 1e-12)
  xmax <- max(2000, rs8$OS.time)
  note8 <- list(
    annotate("text", x = xmax, y = 1, hjust = 1, vjust = 1, size = 7 / .pt, family = "Arial",
             label = sprintf("HR = %.2f (95%% CI %.2f%s%.2f)", m8$HR_high_vs_low, m8$HR_low95, NDASH, m8$HR_high95)),
    annotate("text", x = xmax, y = 0.9, hjust = 1, vjust = 1, size = 7 / .pt, family = "Arial", parse = TRUE,
             label = sprintf('"Log-rank" ~ italic(p) == "%.3f"', m8$logrank_p)))
  km8 <- km_panel(rs8$OS.time, rs8$OS, rs8$RiskGroup, COL_RISK, "Risk group", c(Low = "Low risk", High = "High risk"),
                  note8, "(a)", breaks = seq(0, 2000, 500))
  cat(sprintf("  (a) at risk: %s\n", km8$risk[, paste(n, collapse = "/"), by = group][, paste(group, V1, collapse = "; ")]))
  pA8 <- wrap_elements(full = km8$curve / km8$table + plot_layout(heights = c(3.4, 0.75)))
  auc8 <- c(m8$AUC_1y, m8$AUC_3y, m8$AUC_5y)
  lab8 <- sprintf("%d-year AUC %.3f (95%% CI %.3f%s%.3f)", c(1, 3, 5), auc8,
                  c(m8$AUC_1y_low, m8$AUC_3y_low, m8$AUC_5y_low), NDASH, c(m8$AUC_1y_high, m8$AUC_3y_high, m8$AUC_5y_high))
  r8 <- roc_plot(rs8$OS.time, rs8$OS, rs8$RiskScore, c(365, 1095, 1825), lab8, "(b)", inside = c(0.02, 0.98), just = c(0, 1))
  stopifnot(max(abs(r8$auc - auc8)) < 1e-12)
  # heat map z-scores exactly as stage 11c
  eset <- getGEO(filename = file.path(E, "data/GSE53625_series_matrix.txt.gz"), AnnotGPL = FALSE,
                 destdir = file.path(E, "cache"), GSEMatrix = TRUE, getGPL = FALSE)
  eset <- if (inherits(eset, "list")) eset[[1]] else eset
  xr8 <- exprs(eset)
  e8 <- new.env(); load(file.path(E, "cache/GPL18109_pipe.rda"), envir = e8); ids <- get("GPL18109_pipe", e8)
  ids <- ids[!is.na(ids$symbol) & ids$symbol != "", ]
  genes8 <- tb("Tab10_LASSO_selected_genes.csv")$gene
  ex8 <- t(sapply(genes8, function(g) {
    p <- intersect(as.character(ids$probe_id[ids$symbol == g]), rownames(xr8))
    if (length(p) == 1) xr8[p, rs8$sampleID] else colMeans(xr8[p, rs8$sampleID, drop = FALSE], na.rm = TRUE)
  }))
  z8 <- t(scale(t(ex8)))
  pC8 <- risk3(rs8, z8, genes8, rs8$cutoff_used[1],
               paste0("Cutoff = ", num(rs8$cutoff_used[1], 4), " (TCGA training median)"),
               sprintf("GSE53625 samples sorted by risk score (low \u2192 high, n = %d)", nrow(rs8)), "(c)")
  fig8 <- (pA8 | r8$plot) / pC8 + plot_layout(heights = c(5.6, 5.8))
  save_both(fig8, "Figure8_draft", 13.4)
}

# =====================================================================================================  Figure 9
if ("9" %in% args) {
  cat("Figure 9 (risk score vs clinical and immune features)\n")
  rs <- tb("Tab10_risk_scores_per_sample.csv"); rs[, RiskGroup := factor(RiskGroup, levels = c("Low", "High"))]
  cli <- fread(file.path(E, "data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv"))[sampleID %in% rs$sampleID]
  cli <- merge(cli, rs[, .(sampleID, RiskGroup)], by = "sampleID")
  collapse_stage <- function(s) { s[s == "" | grepl("^\\[", s)] <- NA; out <- rep(NA_character_, length(s))
    out[grepl("Stage I[ABC]?$", s)] <- "I"; out[grepl("Stage II[ABC]?$", s)] <- "II"
    out[grepl("Stage III[ABC]?$", s)] <- "III"; out[grepl("Stage IV[ABC]?$", s)] <- "IV"; out }
  collapse_residual <- function(s) { s[s == "" | grepl("^\\[", s)] <- NA; out <- rep(NA_character_, length(s))
    out[s == "R0"] <- "R0"; out[s %in% c("R1", "R2")] <- "R1+R2"; out[s == "RX"] <- "RX"; out }
  cli[, `:=`(age = suppressWarnings(as.numeric(age_at_initial_pathologic_diagnosis)),
             stage_c = collapse_stage(pathologic_stage), residual_c = collapse_residual(residual_tumor),
             hist_c = ifelse(histological_type == "Esophagus Squamous Cell Carcinoma", "ESCC",
                             ifelse(grepl("^Esophagus Adenocarcinoma", histological_type), "EAC", NA)),
             gender_c = ifelse(gender %in% c("FEMALE", "MALE"), gender, NA),
             alcohol_c = ifelse(alcohol_history_documented %in% c("NO", "YES"), alcohol_history_documented, NA),
             reflux_c = ifelse(reflux_history %in% c("NO", "YES"), reflux_history, NA),
             barretts_c = ifelse(barretts_esophagus %in% c("No", "Yes-UK", "Yes-USA"), barretts_esophagus, NA),
             grade_c = ifelse(neoplasm_histologic_grade %in% c("G1", "G2", "G3", "GX"), neoplasm_histologic_grade, NA),
             tobacco_c = ifelse(tobacco_smoking_history %in% c("1", "2", "3", "4", "Discrepancy"), tobacco_smoking_history, NA))]
  small <- theme(plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5, margin = margin(0, 0, 2, 0)),
                 panel.grid.major.x = element_blank(), plot.margin = margin(2, 3, 2, 2))
  cat_panel <- function(var, title, lev, pal) {
    d <- cli[!is.na(get(var)), .(RiskGroup, v = factor(get(var), levels = lev))]
    pct <- d[, .N, by = .(RiskGroup, v)][, pct := N / sum(N) * 100, by = RiskGroup]
    # no percentages in the bars (user, 11-15): the text cites Figure 9a only for "no significant difference", and at
    # this width the numbers hid the bar composition and did not fit the thin segments
    ggplot(pct, aes(RiskGroup, pct, fill = v)) + geom_col(width = 0.7) +
      scale_fill_manual(values = setNames(pal, lev), name = NULL, drop = TRUE) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.04)), labels = function(x) sprintf("%d%%", x)) +
      labs(title = title, x = NULL, y = NULL) + theme_fig + small +
      theme(legend.position = "bottom", legend.key.size = unit(2.4, "mm"), legend.margin = margin(-3, 0, 0, 0),
            legend.text = element_text(size = 7, margin = margin(l = 0.6, unit = "mm")), legend.key.spacing.x = unit(1, "mm"))
  }
  p_age <- ggplot(cli[!is.na(age)], aes(RiskGroup, age, fill = RiskGroup)) +
    geom_boxplot(width = 0.55, outlier.size = 0.4, alpha = 0.8, linewidth = 0.3) +
    scale_fill_manual(values = COL_RISK, guide = "none") + labs(title = "age", x = NULL, y = "Age (years)", tag = "(a)") +
    theme_fig + small
  pa9 <- list(p_age,
    cat_panel("gender_c", "gender", c("FEMALE", "MALE"), c("#E41A1C", "#377EB8")),
    cat_panel("alcohol_c", "alcohol", c("NO", "YES"), c("#E41A1C", "#377EB8")),
    cat_panel("tobacco_c", "tobacco", c("1", "2", "3", "4", "Discrepancy"), c("#FB8072", "#FDB462", "#FFFFB3", "#80B1D3", "#B3DE69")),
    cat_panel("reflux_c", "reflux", c("NO", "YES"), c("#E41A1C", "#377EB8")),
    cat_panel("barretts_c", "barretts", c("No", "Yes-UK", "Yes-USA"), c("#66C2A5", "#FC8D62", "#8DA0CB")),
    cat_panel("hist_c", "histology", c("ESCC", "EAC"), c("#1B9E77", "#D95F02")),
    cat_panel("grade_c", "grade", c("G1", "G2", "G3", "GX"), c("#FFFFB3", "#FDB462", "#FB8072", "#BEBADA")),
    cat_panel("stage_c", "pathologic stage", c("I", "II", "III", "IV"), c("#9ECAE1", "#6BAED6", "#3182BD", "#08519C")),
    cat_panel("residual_c", "residual tumor", c("R0", "R1+R2", "RX"), c("#A6CEE3", "#FB9A99", "#FDBF6F")))
  pa_grid <- wrap_elements(full = wrap_plots(pa9, nrow = 2))

  CAT_COL <- c(ESTIMATE = "#1B9E77", CIBERSORT = "#D95F02", ssGSEA_Bindea = "#7570B3", Immune_checkpoint = "#E7298A", HLA = "#66A61E")
  CAT_LAB <- c(ESTIMATE = "ESTIMATE", CIBERSORT = "CIBERSORTx", ssGSEA_Bindea = "ssGSEA (Bindea)",
               Immune_checkpoint = "Immune checkpoint", HLA = "HLA")
  feat_it <- function(f, cat) cat == "HLA" | (cat == "Immune_checkpoint" & f %in% IT_SYMBOLS)
  sp <- tb("Tab12_RS_immune_spearman.csv")[!is.na(rho)][order(p_value)][1:20]
  z <- atanh(sp$rho); se <- 1 / sqrt(sp$n - 3)
  sp[, `:=`(lo = tanh(z - qnorm(0.975) * se), hi = tanh(z + qnorm(0.975) * se), sig = sig_from_p(p_value),
            label = factor(feature, levels = rev(feature)))]
  it9 <- sp[feat_it(feature, feature_category), feature]
  rng <- max(abs(c(sp$lo, sp$hi)))
  sp[, star_x := fifelse(rho >= 0, hi, lo) + sign(rho) * rng * 0.03]
  pb9 <- ggplot(sp, aes(rho, label, colour = feature_category)) +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.45, linewidth = 0.35, orientation = "y", show.legend = FALSE) +
    geom_point(size = 1.2) +
    geom_text(aes(x = star_x, label = sig, hjust = fifelse(rho >= 0, 0, 1)), colour = "black", size = 7.5 / .pt,
              vjust = 0.72, family = "Arial") +
    scale_colour_manual(values = CAT_COL, labels = CAT_LAB, name = NULL) +
    scale_x_continuous(expand = expansion(mult = 0.1), labels = lab_num) +
    scale_y_discrete(labels = function(l) italic_labels(l, l %in% it9)) +
    guides(colour = guide_legend(ncol = 3, byrow = TRUE)) +
    labs(x = expression("Spearman" ~ rho ~ "(risk score vs feature) with 95% CI"), y = NULL, tag = "(b)") + theme_fig +
    theme(panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey92"), legend.position = "bottom",
          legend.margin = margin(-2, 0, 0, 0), legend.key.size = unit(2.6, "mm"))

  wx <- tb("Tab12_RS_immune_wilcoxon.csv")[!is.na(p_value)][order(p_value)][1:8]
  est <- tb("Tab8A_ESTIMATE_scores.csv"); cib <- tb("Tab8B_CIBERSORT_results.csv"); ssg <- tb("Tab8C_ssGSEA_Bindea.csv")
  cmap <- setNames(CP$gene_symbol, CP$gene_display)
  feat_vals <- function(f, cat) {
    if (cat == "ESTIMATE") return(data.table(sampleID = est$sampleID, value = est[[f]]))
    if (cat == "CIBERSORT") return(data.table(sampleID = cib$sampleID, value = cib[[f]]))
    if (cat == "ssGSEA_Bindea") { m <- which(ssg$cell_type == sub(" \\(Bindea\\)$", "", f)); sc <- setdiff(colnames(ssg), "cell_type")
      return(data.table(sampleID = sc, value = as.numeric(ssg[m, ..sc]))) }
    sym <- if (cat == "Immune_checkpoint") cmap[[f]] else f
    data.table(sampleID = rs$sampleID, value = expr_row(sym, rs$sampleID))
  }
  l9 <- rbindlist(lapply(seq_len(nrow(wx)), function(i) feat_vals(wx$feature[i], wx$feature_category[i])[, feature := wx$feature[i]]))
  l9 <- merge(l9, rs[, .(sampleID, RiskGroup)], by = "sampleID")
  l9[, feature := factor(feature, levels = wx$feature)]
  chk <- l9[, .(p = wilcox.test(value ~ RiskGroup, exact = FALSE)$p.value), by = feature]
  stopifnot(max(abs(chk$p - wx$p_value[match(chk$feature, wx$feature)])) < 1e-10)   # same values as Tab12
  s9 <- merge(wx[, .(feature = factor(feature, levels = wx$feature), sig = sig_from_p(p_value))],
              l9[, .(y = min(value) + (max(value) - min(value)) * 0.95), by = feature], by = "feature")
  it9c <- wx[feat_it(feature, feature_category), feature]
  pc9 <- ggplot(l9, aes(RiskGroup, value, fill = RiskGroup)) +
    geom_boxplot(width = 0.55, outlier.size = 0.3, alpha = 0.85, colour = "black", linewidth = 0.3) +
    geom_text(data = s9, aes(x = 1.5, y = y, label = sig), inherit.aes = FALSE, size = 7.5 / .pt, vjust = 0.6, family = "Arial") +
    facet_wrap(~feature, scales = "free_y", ncol = 4,     # names wrapped; none of the eight is a gene symbol (asserted)
               labeller = labeller(feature = label_wrap_gen(13))) +
    scale_fill_manual(values = COL_RISK, name = "Risk group") + scale_y_continuous(labels = lab_num) +
    labs(x = NULL, y = "Feature value", tag = "(c)") + theme_fig +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0), panel.grid.major.x = element_blank(),
          strip.text = element_text(size = 7), panel.spacing = unit(1.2, "mm"))
  cat(sprintf("  (b) top 20 by p, %d italic; (c) top 8: %s\n", length(it9), paste(wx$feature, collapse = ", ")))
  stopifnot(length(it9c) == 0)
  fig9 <- pa_grid / ((pb9 | pc9) + plot_layout(widths = c(0.45, 1))) + plot_layout(heights = c(8.2, 7.0))
  save_both(fig9, "Figure9_draft", 16.6)
}

# =====================================================================================================  Figure 3b
if ("3" %in% args) {
  cat("Figure 3b (copy-number frequency)\n")
  fq <- tb("Tab5_TCGA_CNV_frequency.csv")
  gord <- tb("Tab4_TCGA_mutation_frequency.csv")$gene
  stopifnot(setequal(fq$gene, gord))
  pd3 <- rbind(fq[, .(gene, alt = "gain", pct = gain_pct)], fq[, .(gene, alt = "high_amp", pct = high_amp_pct)],
               fq[, .(gene, alt = "loss", pct = -loss_pct)], fq[, .(gene, alt = "deep_del", pct = -deep_del_pct)])
  pd3[, `:=`(gene = factor(gene, levels = rev(gord)), alt = factor(alt, levels = c("gain", "high_amp", "loss", "deep_del")))]
  lim <- ceiling(max(abs(pd3[, sum(pct), by = .(gene, sign(pct))]$V1)) / 10) * 10
  br <- seq(-floor(lim / 30) * 30, floor(lim / 30) * 30, 30)
  p3b <- ggplot(pd3, aes(pct, gene, fill = alt)) +
    geom_col(width = 0.7, colour = NA, orientation = "y") +
    geom_vline(xintercept = 0, linewidth = 0.35) +
    scale_fill_manual(values = c(high_amp = "#CA0020", gain = "#F4A582", loss = "#92C5DE", deep_del = "#0571B0"),
                      breaks = c("high_amp", "gain", "loss", "deep_del"), name = NULL,
                      labels = c("High-level amplification (+2)", "Gain (+1)", paste0("Heterozygous loss (", MINUS, "1)"),
                                 paste0("Deep deletion (", MINUS, "2)"))) +
    scale_x_continuous(limits = c(-lim, lim), breaks = br,
                       labels = ifelse(br == 0, "0%", ifelse(br > 0, sprintf("+%d%%", br), paste0(MINUS, abs(br), "%")))) +
    scale_y_discrete(labels = function(l) italic_labels(l, TRUE)) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(x = "Frequency (%)", y = NULL, tag = "(b)") + theme_fig +
    theme(panel.grid.major.y = element_blank(), legend.position = "top", legend.key.size = unit(2.6, "mm"),
          legend.margin = margin(0, 0, -2, 0), legend.justification = "left")

  # (c) TCGA tumours vs GTEx mucosa, figures/remake/v3_fig2C_violin.R (raw-p stars from the same Wilcoxon tests)
  gt <- fread(file.path(E, "data/expr_GTEx_esophagus.csv")); setnames(gt, 1, "gene_id")
  muc <- intersect(fread(file.path(E, "data/GTEx_esophagus_subtype_map.csv"))[subtype == "Esophagus - Mucosa", sample], colnames(gt))
  tum <- setdiff(colnames(expr), "gene_id")
  vl <- rbindlist(lapply(gord, function(g) rbind(
    data.table(gene = g, group = "Normal", value = as.numeric(gt[gene_id == g][1, ..muc])),
    data.table(gene = g, group = "Tumor", value = as.numeric(expr[gene_id == g][1, ..tum])))))
  vl[, `:=`(group = factor(group, levels = c("Normal", "Tumor")), gene = factor(gene, levels = gord))]
  pv <- vl[, .(p = wilcox.test(value ~ group, exact = FALSE)$p.value), by = gene]
  sv <- merge(pv[, .(gene, sig = sig_from_p(p))], vl[, .(y = min(value) + (max(value) - min(value)) * 0.97), by = gene], by = "gene")
  cat(sprintf("  (c) GTEx mucosa n = %d, tumours n = %d; stars on %d of 18\n", length(muc), length(tum), sum(sv$sig != "")))
  p3c <- ggplot(vl, aes(group, value, fill = group)) +
    geom_violin(alpha = 0.6, colour = "black", linewidth = 0.3, trim = TRUE) +
    geom_boxplot(width = 0.18, fill = "white", colour = "black", outlier.size = 0.2, linewidth = 0.3) +
    geom_text(data = sv, aes(x = 1.5, y = y, label = sig), inherit.aes = FALSE, size = 7.5 / .pt, vjust = 0.6, family = "Arial") +
    facet_wrap(~gene, scales = "free_y", ncol = 6,
               labeller = as_labeller(function(f) sprintf('italic("%s")', f), label_parsed)) +
    scale_fill_manual(values = c(Normal = "#4DAF4A", Tumor = "#E41A1C"), name = NULL) +
    scale_y_continuous(labels = lab_num) +
    labs(x = NULL, y = expression(log[2] * "(RSEM + 1) expression"), tag = "(c)") + theme_fig +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0), panel.grid.major.x = element_blank(),
          strip.text = element_text(size = 7), panel.spacing = unit(1.2, "mm"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())   # groups by colour (legend); names did not fit

  # (a) oncoprint, scripts/04_TCGA_mutation_oncoprint.R (n = 184 primary tumours with MC3 calls)
  mc3 <- fread(file.path(E, "data/ESCA_mc3.csv"))
  cli3 <- fread(file.path(E, "data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv"))
  fmg <- tb("Tab1_2_FMG_list_primary.csv")$gene
  keep <- intersect(cli3[sample_type == "Primary Tumor", sampleID], unique(mc3$sample))
  eff <- c(Missense_Mutation = "Missense", Nonsense_Mutation = "Nonsense", Frame_Shift_Del = "Frameshift",
           Frame_Shift_Ins = "Frameshift", In_Frame_Del = "Inframe_indel", In_Frame_Ins = "Inframe_indel",
           Splice_Site = "Splice", Translation_Start_Site = "Other_functional", Nonstop_Mutation = "Other_functional")
  mf <- mc3[sample %in% keep & gene %in% fmg & effect %in% names(eff)][, alt := eff[effect]]
  ml <- mf[, .(alt = paste(unique(alt), collapse = ";")), by = .(sample, gene)]
  mat <- matrix("", length(fmg), length(keep), dimnames = list(fmg, keep))
  mat[cbind(ml$gene, ml$sample)] <- ml$alt
  t4 <- tb("Tab4_TCGA_mutation_frequency.csv")
  stopifnot(length(keep) == 184, all(rowSums(mat[t4$gene, ] != "") == t4$n_mut))   # = Table 2 / Tab4
  ALT <- c(Missense = "#377EB8", Nonsense = "#E41A1C", Frameshift = "#984EA3", Inframe_indel = "#FF7F00",
           Splice = "#4DAF4A", Other_functional = "#A65628")
  ALT <- ALT[names(ALT) %in% unlist(strsplit(ml$alt, ";"))]       # legend: the types that occur (published one lacked Nonsense)
  cat(sprintf("  (a) %d tumours; alteration types: %s\n", length(keep), paste(names(ALT), collapse = ", ")))
  af <- c(list(background = function(x, y, w, h) grid.rect(x, y, w * 0.92, h * 0.92, gp = gpar(fill = "#CCCCCC", col = NA))),
          lapply(ALT, function(cl) { force(cl); function(x, y, w, h) grid.rect(x, y, w * 0.92, h * 0.92, gp = gpar(fill = cl, col = NA)) }))
  hist3 <- cli3[match(keep, sampleID), fifelse(histological_type == "Esophagus Squamous Cell Carcinoma", "ESCC",
                                              fifelse(histological_type == "Esophagus Adenocarcinoma, NOS", "EAC", "Other"))]
  g7 <- gpar(fontsize = 7, fontfamily = "Arial")
  ht3 <- oncoPrint(mat[t4$gene, ], alter_fun = af, alter_fun_is_vectorized = FALSE, col = ALT,
                   top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(height = unit(1.2, "cm"),
                                                      axis_param = list(gp = g7)), show_annotation_name = FALSE),
                   bottom_annotation = HeatmapAnnotation(Histology = hist3, col = list(Histology = c(ESCC = "#1B9E77", EAC = "#D95F02", Other = "#7570B3")[unique(hist3)]),
                                                         annotation_name_side = "left", annotation_name_gp = g7,
                                                         simple_anno_size = unit(2.5, "mm"),
                                                         annotation_legend_param = list(Histology = list(title_gp = g7, labels_gp = g7))),
                   right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(width = unit(1.6, "cm"), axis_param = list(gp = g7)),
                                                    show_annotation_name = FALSE),
                   row_order = t4$gene, remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                   pct_side = "right", row_names_side = "left", pct_gp = g7,
                   row_names_gp = gpar(fontsize = 7, fontface = "italic", fontfamily = "Arial"), column_title = NULL,
                   heatmap_legend_param = list(title = "Alteration", at = names(ALT), labels = sub("_", " ", names(ALT)),
                                               title_gp = g7, labels_gp = g7))
  H3 <- 16.2; HA <- 7.0                                               # cm: whole figure, row (a)
  bottom3 <- (p3b | p3c) + plot_layout(widths = c(0.36, 1))
  draw3 <- function() {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(HA, H3 - HA), "cm"))))
    pushViewport(viewport(layout.pos.row = 1))
    draw(ht3, newpage = FALSE, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right",
         padding = unit(c(1, 3, 6, 2), "mm"))
    grid.text("(a)", x = unit(1, "mm"), y = unit(1, "npc") - unit(1, "mm"), just = c("left", "top"),
              gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"))
    popViewport()
    pushViewport(viewport(layout.pos.row = 2)); print(bottom3, newpage = FALSE); popViewport()
  }
  cairo_pdf(file.path(OUT, "Figure3_draft.pdf"), width = W_IN, height = H3 / 2.54); draw3(); invisible(dev.off())
  agg_png(file.path(OUT, "Figure3_draft.png"), width = W_IN, height = H3 / 2.54, units = "in", res = RES); draw3(); invisible(dev.off())
  cat(sprintf("  Figure3_draft: 16.5 x %.1f cm\n", H3))
}

# =====================================================================================================  Figure 11a
if ("11" %in% args) {
  cat("Figure 11a (TIDE metrics, TIDEpy re-run of stage 14a)\n")
  rs <- tb("Tab10_risk_scores_per_sample.csv"); rs[, RiskGroup := factor(RiskGroup, levels = c("Low", "High"))]
  tide <- fread(file.path(E, "data/esca_tide_output_hgnc_centered.txt"), header = TRUE); setnames(tide, 1, "sampleID")
  tt <- tb("Tab14A_TIDE_metrics_test.csv")
  col_of <- c(MSI.Score = "MSI Score", TAM.M2 = "TAM M2")
  tl <- rbindlist(lapply(tt$metric, function(m) {
    col <- if (m %in% names(col_of)) col_of[[m]] else m
    tide[, .(sampleID, metric = m, value = suppressWarnings(as.numeric(get(col))))][!is.na(value)]
  }))
  tl <- merge(tl, rs[, .(sampleID, RiskGroup)], by = "sampleID")
  tl[, metric := factor(metric, levels = tt$metric)]
  s11 <- merge(tt[, .(metric = factor(metric, levels = tt$metric), sig = sig_from_p(p_value))],
               tl[, .(y = min(value) + (max(value) - min(value)) * 0.95), by = metric], by = "metric")
  cat(sprintf("  stars: %s\n", paste(s11[sig != "", paste(metric, sig)], collapse = ", ")))
  p11 <- ggplot(tl, aes(RiskGroup, value, fill = RiskGroup)) +
    geom_boxplot(width = 0.55, outlier.size = 0.3, alpha = 0.85, colour = "black", linewidth = 0.3) +
    geom_text(data = s11, aes(x = 1.5, y = y, label = sig), inherit.aes = FALSE, size = 7.5 / .pt, vjust = 0.6, family = "Arial") +
    # three columns (four were too narrow for 7-pt titles); the legend fills the empty twelfth cell
    facet_wrap(~metric, scales = "free_y", ncol = 3,
               labeller = as_labeller(function(f) ifelse(f == "CD274", 'italic("CD274")', sprintf('"%s"', f)), label_parsed)) +
    scale_fill_manual(values = COL_RISK, name = "Risk group") + scale_y_continuous(labels = lab_num) +
    labs(x = NULL, y = "Score", tag = "(a)") + theme_fig +
    theme(legend.position = "inside", legend.position.inside = c(0.86, 0.11), panel.grid.major.x = element_blank(),
          strip.text = element_text(size = 7, face = "bold"), panel.spacing = unit(1.2, "mm"))
  # (b) IPS: per-sample values of figures/remake/v3_fig10C_IPS.R (IOBR, cached), tests of Tab14B
  ips <- readRDS(file.path(E, "cache/v3_fig10C_ips_long.rds"))
  t14b <- tb("Tab14B_IPS_metrics_test.csv")
  il <- merge(melt(ips[, c("sampleID", t14b$metric), with = FALSE], id.vars = "sampleID", variable.name = "metric"),
              rs[, .(sampleID, RiskGroup)], by = "sampleID")
  il[, metric := factor(metric, levels = t14b$metric)]
  chk <- il[, .(p = wilcox.test(value ~ RiskGroup, exact = FALSE)$p.value), by = metric]
  stopifnot(max(abs(chk$p - t14b$p_value[match(chk$metric, t14b$metric)])) < 1e-10)   # cache = the tested values
  sb <- merge(t14b[, .(metric = factor(metric, levels = t14b$metric), sig = sig_from_p(p_value))],
              il[, .(y = min(value) + (max(value) - min(value)) * 0.95), by = metric], by = "metric")
  p11b <- ggplot(il, aes(RiskGroup, value, fill = RiskGroup)) +
    geom_violin(alpha = 0.6, colour = "black", linewidth = 0.3, trim = TRUE) +
    geom_boxplot(width = 0.18, fill = "white", colour = "black", outlier.size = 0.3, linewidth = 0.3) +
    geom_text(data = sb, aes(x = 1.5, y = y, label = sig), inherit.aes = FALSE, size = 7.5 / .pt, vjust = 0.6, family = "Arial") +
    facet_wrap(~metric, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = COL_RISK, name = "Risk group", guide = "none") + scale_y_continuous(labels = lab_num) +
    labs(x = NULL, y = "IPS score", tag = "(b)") + theme_fig +
    theme(panel.grid.major.x = element_blank(), panel.spacing = unit(1.2, "mm"))
  # (c) oncoPredict: ln(IC50) of Tab14C_oncoPredict_IC50_full.csv, ten smallest p of Tab14C_oncoPredict_top30_drugs.csv
  dt14 <- tb("Tab14C_oncoPredict_top30_drugs.csv")[order(p_value)][1:10]
  ic <- merge(tb("Tab14C_oncoPredict_IC50_full.csv"), rs[, .(sampleID, RiskGroup)], by = "sampleID")
  dl <- rbindlist(lapply(dt14$drug, function(d) data.table(RiskGroup = ic$RiskGroup, drug = d, value = ic[[d]])))[!is.na(value)]
  dl[, drug := factor(drug, levels = dt14$drug)]
  chk <- dl[, .(p = wilcox.test(value ~ RiskGroup, exact = FALSE)$p.value), by = drug]
  stopifnot(max(abs(chk$p - dt14$p_value[match(chk$drug, dt14$drug)])) < 1e-10)
  sc <- merge(dt14[, .(drug = factor(drug, levels = dt14$drug), sig = sig_from_p(p_value))],
              dl[, .(y = min(value) + (max(value) - min(value)) * 0.95), by = drug], by = "drug")
  p11c <- ggplot(dl, aes(RiskGroup, value, fill = RiskGroup)) +
    geom_boxplot(width = 0.55, outlier.size = 0.3, alpha = 0.85, colour = "black", linewidth = 0.3) +
    geom_text(data = sc, aes(x = 1.5, y = y, label = sig), inherit.aes = FALSE, size = 7.5 / .pt, vjust = 0.6, family = "Arial") +
    facet_wrap(~drug, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = COL_RISK, name = "Risk group", guide = "none") + scale_y_continuous(labels = lab_num) +
    labs(x = NULL, y = "Predicted ln(IC50)", tag = "(c)") + theme_fig +
    theme(panel.grid.major.x = element_blank(), strip.text = element_text(size = 7, face = "plain"),
          panel.spacing = unit(1.2, "mm"))
  cat(sprintf("  stars (b): %s; (c): %s\n", paste(sb[sig != "", paste(metric, sig)], collapse = ", "),
              paste(sc[, paste0(drug, " ", sig)], collapse = ", ")))
  fig11 <- ((p11 | p11b) + plot_layout(widths = c(1, 1.1))) / p11c + plot_layout(heights = c(7.0, 5.0))
  save_both(fig11, "Figure11_draft", 15.6)
}

# =====================================================================================================  Figure 10b
if ("10" %in% args) {
  cat("Figure 10b (calibration; published data areas, new axes)\n")
  ext <- fread(file.path(OUT, "fig10b_extents.csv"))
  cal <- lapply(1:3, function(k) {
    img <- readPNG(file.path(OUT, ext$file[k]))
    ratio <- (dim(img)[1] / (ext$y_top[k] - ext$y_bottom[k])) / (dim(img)[2] / (ext$x_right[k] - ext$x_left[k]))
    ggplot(data.frame(x = c(-0.04, 1.04), y = c(-0.04, 1.04)), aes(x, y)) +
      geom_blank() +                                   # trains the scales, otherwise no tick labels are drawn
      annotation_raster(img, xmin = ext$x_left[k], xmax = ext$x_right[k], ymin = ext$y_bottom[k], ymax = ext$y_top[k],
                        interpolate = TRUE) +
      scale_x_continuous(breaks = seq(0, 1, 0.2), labels = sprintf("%.1f", seq(0, 1, 0.2))) +
      scale_y_continuous(breaks = seq(0, 1, 0.2), labels = sprintf("%.1f", seq(0, 1, 0.2))) +
      coord_fixed(ratio = ratio, xlim = c(-0.04, 1.04), ylim = c(-0.04, 1.04), expand = FALSE) +
      labs(title = c("1-year", "2-year", "3-year")[k], x = "Nomogram-predicted survival",
           y = if (k == 1) "Observed survival" else NULL) +
      theme_fig + theme(panel.grid.major = element_blank(), plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5))
  })
  pb10 <- wrap_plots(cal, nrow = 1)
  save_panel(pb10, "Figure10b_panel", 4000, 1344)

  # (a) nomogram: the cph model and nomogram() call of scripts/helpers/13_fix_Fig9_nomogram_calibration.R (deterministic),
  # drawn at the printed size (axis labels 7 pt, variable names 7.5 pt; device pointsize 12)
  suppressPackageStartupMessages(library(rms))
  rs <- tb("Tab10_risk_scores_per_sample.csv")
  c10 <- fread(file.path(E, "data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv"))
  st <- function(s) { s[s == "" | grepl("^\\[", s)] <- NA; o <- rep(NA_character_, length(s))
    o[grepl("^Stage I[ABC]?$|^Stage I$", s)] <- "I"; o[grepl("^Stage II", s)] <- "II"; o[grepl("^Stage III", s)] <- "III"
    o[grepl("^Stage IV", s)] <- "IV"; o }
  rr <- function(r) { r[r == "" | grepl("^\\[", r)] <- NA; o <- rep(NA_character_, length(r))
    o[r == "R0"] <- "R0"; o[r %in% c("R1", "R2")] <- "R1+R2"; o }
  c10[, `:=`(stage_collapsed = st(pathologic_stage), residual_collapsed = rr(residual_tumor),
             histology_collapsed = ifelse(histological_type == "Esophagus Squamous Cell Carcinoma", "ESCC",
                                          ifelse(histological_type == "Esophagus Adenocarcinoma, NOS", "EAC", NA)))]
  ph <- merge(rs[, .(sampleID, OS, OS.time, RiskScore)],
              c10[, .(sampleID, stage_collapsed, residual_collapsed, histology_collapsed)], by = "sampleID")
  ph <- ph[complete.cases(ph[, .(RiskScore, stage_collapsed, residual_collapsed, histology_collapsed)])]
  ph[, `:=`(stage_collapsed = factor(stage_collapsed, levels = c("I", "II", "III", "IV")),
            residual_collapsed = factor(residual_collapsed, levels = c("R0", "R1+R2")),
            histology_collapsed = factor(histology_collapsed, levels = c("ESCC", "EAC")))]
  phd <- as.data.frame(ph)
  label(phd$RiskScore) <- "Risk score"; label(phd$histology_collapsed) <- "Histology"
  label(phd$stage_collapsed) <- "Pathologic stage"; label(phd$residual_collapsed) <- "Residual tumor"
  dd <- datadist(phd); options(datadist = "dd")
  fit10 <- cph(Surv(OS.time, OS) ~ RiskScore + histology_collapsed + stage_collapsed + residual_collapsed,
               data = phd, x = TRUE, y = TRUE, surv = TRUE, time.inc = 365)
  sv <- Survival(fit10)
  nom <- nomogram(fit10, fun = list(function(x) sv(365, x), function(x) sv(730, x), function(x) sv(1095, x)),
                  funlabel = c("1-year survival", "2-year survival", "3-year survival"), lp = TRUE,
                  fun.at = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 0.95))
  cat(sprintf("  (a) nomogram on %d patients (%d deaths); coefficients: %s\n", nrow(phd), sum(phd$OS),
              paste(sprintf("%.4f", coef(fit10)), collapse = " ")))
  HA10 <- 12.0; H10 <- HA10 + 1344 / PXCM                           # cm
  draw10 <- function() {
    plot.new()
    par(fig = c(0, 1, 1 - HA10 / H10, 1), mar = c(0.5, 1, 1.5, 1), family = "Arial", new = TRUE)
    plot(nom, xfrac = 0.25, cex.axis = 7 / 12, cex.var = 7.5 / 12, lmgp = 0.15)
    print(pb10, vp = viewport(x = 0.5, y = 0, width = 1, height = (H10 - HA10) / H10, just = c("centre", "bottom")))
    tg <- gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial")
    grid.text("(a)", x = unit(1, "mm"), y = unit(1, "npc") - unit(1, "mm"), just = c("left", "top"), gp = tg)
    grid.text("(b)", x = unit(1, "mm"), y = unit((H10 - HA10) / H10, "npc") + unit(1, "mm"), just = c("left", "bottom"), gp = tg)
  }
  cairo_pdf(file.path(OUT, "Figure10_draft.pdf"), width = W_IN, height = H10 / 2.54); draw10(); invisible(dev.off())
  agg_png(file.path(OUT, "Figure10_draft.png"), width = W_IN, height = H10 / 2.54, units = "in", res = RES); draw10(); invisible(dev.off())
  cat(sprintf("  Figure10_draft: 16.5 x %.1f cm\n", H10))
}

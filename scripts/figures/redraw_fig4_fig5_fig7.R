###############################################################################
## Figures 4, 5 and 7 redrawn at the size they are printed (16.5 cm wide,
## text 7-7.5 pt, panel tags 10 pt bold), as was done for Figures 2 and 6 (redraw_fig2_fig6.R).
##
##   Rscript redraw_fig4_fig5_fig7.R [4] [5] [7]      (R 4.4.1; default: all three)
##
## Sources (read-only). Nothing reported in the text is recomputed differently:
##   Figure 4  (a)-(c) ccp_rerun.rds = scripts/07 STEP 2-6 run again with the same inputs and seed
##             (rerun_ccp_fig4.R: k = 2 assignment identical to Tab7_cluster_assignment.csv; PAC and delta area
##             equal to the points of figures/Fig3B_PAC_score.pdf to 4 decimals);
##             (d)-(f) scripts/07 STEP 10-13 on the same z-scored 18-FMG matrix; (g) figures/remake/v3_fig3G_violin.R
##   Figure 5  (a) Tab8A (script 08 Fig4A), stars from raw p;
##             (b)-(e) Tab8B-E + Hodges-Lehmann CIs exactly as figures/remake/v3_fig4_BCDE_lollipop_forest.R
##   Figure 7  (a)-(b) lasso_rerun.rds = the cv.glmnet step of scripts/10 run again (rerun_lasso_fig7.R:
##             coefficients and risk scores equal to Tab10 to 1e-15); (c)-(e) Tab10 risk scores, scripts/10 code
## Earlier pixel corrections now in the code: gene symbols in italics, "x 10^-k" and en dashes,
## AUC CIs in Fig 7d, raw-p stars in Fig 5a.
###############################################################################
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ragg); library(ggrepel); library(grid)
  library(survival); library(ComplexHeatmap); library(circlize); library(timeROC); library(glmnet)
})
set.seed(1262118388)
# grid.grabExpr measures text on a pdf(NULL) device; Arial has the Helvetica metrics, so register it there
pdfFonts(Arial = pdfFonts()$Helvetica)
E <- getwd()   # repository root (see README.md)
OUT <- file.path(E, "figures", "final")
W_IN <- 16.5 / 2.54
RES <- 4000 / W_IN
# a script that sources this file for its helpers sets which_figs <- "none" first
if (!exists("which_figs")) {
  which_figs <- commandArgs(trailingOnly = TRUE)
  if (!length(which_figs)) which_figs <- c("4", "5", "7")
}
tb <- function(f) fread(file.path(E, "tables", f))

theme_fig <- theme_bw(base_size = 7, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
        panel.border = element_rect(colour = "black", linewidth = 0.4, fill = NA),
        axis.text = element_text(size = 7, colour = "black"),
        axis.title = element_text(size = 7.5, colour = "black"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 7.5),
        legend.key.size = unit(3.2, "mm"),
        strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 0.4),
        strip.text = element_text(size = 7.5, face = "bold", colour = "black"),
        plot.tag = element_text(size = 10, face = "bold", family = "Arial"),
        plot.margin = margin(3, 6, 3, 3))
tags <- function(n) plot_annotation(tag_levels = list(sprintf("(%s)", letters[seq_len(n)])))
save_both <- function(p, stem, h_cm) {
  ggsave(file.path(OUT, paste0(stem, ".pdf")), p, width = W_IN, height = h_cm / 2.54, device = cairo_pdf)
  agg_png(file.path(OUT, paste0(stem, ".png")), width = W_IN, height = h_cm / 2.54, units = "in", res = RES)
  print(p); invisible(dev.off())
  cat(sprintf("  %s: 16.5 x %.1f cm\n", stem, h_cm))
}
sig_from_p <- function(p) fifelse(is.na(p), "", fifelse(p < 0.001, "***", fifelse(p < 0.01, "**", fifelse(p < 0.05, "*", ""))))
COL_CL <- c(C1 = "#1B9E77", C2 = "#D95F02")
COL_RISK <- c(Low = "#377EB8", High = "#E41A1C")
MINUS <- "\u2212"; NDASH <- "\u2013"
num <- function(x, d) sub("^-", MINUS, formatC(x, format = "f", digits = d))
# "a x 10^-b" as plotmath, e.g. 7.65e-06 -> 7.65 %*% 10^-6
sci <- function(p) { e <- floor(log10(p)); sprintf("%.2f %%*%% 10^%d", p / 10^e, e) }
italic_labels <- function(lab, it) parse(text = ifelse(rep_len(it, length(lab)), sprintf('italic("%s")', lab), sprintf('"%s"', lab)))
lab_num <- scales::label_number(style_negative = "minus", big.mark = "")   # true minus sign on numeric axes

expr <- fread(file.path(E, "data/expr_ESCA_tumor.csv")); setnames(expr, 1, "gene_id")
expr_row <- function(g, samples) as.numeric(expr[match(g, gene_id), ..samples])
clu <- tb("Tab7_cluster_assignment.csv")

# KM curve + number-at-risk table as one panel; the tag goes on the curve (wrap_elements would put it over the axis title)
km_panel <- function(time, event, group, cols, legend_title, labs_grp, note, tag, breaks = seq(0, 4000, 1000)) {
  d <- data.frame(time, event, group = factor(group, levels = names(cols)))
  fit <- survfit(Surv(time, event) ~ group, data = d)
  st <- factor(rep(names(cols), fit$strata), levels = names(cols))
  cur <- rbindlist(lapply(names(cols), function(g) {
    i <- st == g; data.table(group = factor(g, levels = names(cols)), time = c(0, fit$time[i]), surv = c(1, fit$surv[i]))
  }))
  cen <- data.table(group = st, time = fit$time, surv = fit$surv)[fit$n.censor > 0]
  risk <- summary(fit, times = breaks, extend = TRUE)
  rt <- data.table(group = factor(sub("group=", "", risk$strata), levels = rev(names(cols))), time = risk$time, n = risk$n.risk)
  xl <- c(0, max(breaks, time))       # never cut follow-up beyond the last break (scale limits drop data)
  p_curve <- ggplot(cur, aes(time, surv, colour = group)) +
    geom_step(linewidth = 0.5) +
    geom_point(data = cen, shape = 3, size = 0.9, stroke = 0.35) +
    note +
    scale_colour_manual(values = cols, labels = labs_grp, name = legend_title) +
    scale_x_continuous(breaks = breaks, limits = xl, expand = expansion(mult = 0.03)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Time (days)", y = "Overall survival", tag = tag) + theme_fig +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0), legend.key = element_blank())
  p_tab <- ggplot(rt, aes(time, group, label = n, colour = group)) +
    geom_text(size = 7 / .pt, family = "Arial", show.legend = FALSE) +
    scale_colour_manual(values = cols) +
    scale_y_discrete(labels = labs_grp[levels(rt$group)]) +
    scale_x_continuous(breaks = breaks, limits = xl, expand = expansion(mult = 0.03)) +
    labs(x = NULL, y = NULL, title = "Number at risk") + theme_fig +
    theme(panel.grid.major = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
          axis.text.y = element_text(colour = cols[levels(rt$group)], face = "bold"),
          plot.title = element_text(size = 7, margin = margin(0, 0, 1, 0)), panel.border = element_blank())
  list(curve = p_curve, table = p_tab, fit = fit, risk = rt)
}

# =====================================================================================================  Figure 5
if ("5" %in% which_figs) {
  cat("Figure 5 (tumour immune microenvironment)\n")
  est <- tb("Tab8A_ESTIMATE_scores.csv")
  long <- melt(est, id.vars = c("sampleID", "cluster", "histology"),
               measure.vars = c("StromalScore", "ImmuneScore", "ESTIMATEScore"), variable.name = "score")
  st <- tb("Tab8A_ESTIMATE_stats.csv")
  st[, lab := fifelse(sig_from_p(p_value) == "", "ns", sig_from_p(p_value))]
  st <- merge(st, long[, .(y = max(value)), by = score], by = "score")
  st[, `:=`(y = y + 0.06 * abs(y), score = factor(score, levels = levels(long$score)))]
  cat(sprintf("  (a) raw p: %s\n", paste(sprintf("%s %.2g %s", st$score, st$p_value, st$lab), collapse = "; ")))
  pa <- ggplot(long, aes(cluster, value, fill = cluster)) +
    geom_violin(width = 0.8, alpha = 0.6, colour = "black", linewidth = 0.3, scale = "width") +
    geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3, colour = "black", linewidth = 0.3) +
    geom_text(data = st, aes(x = 1.5, y = y, label = lab), inherit.aes = FALSE, size = 7.5 / .pt,
              fontface = "bold", family = "Arial") +
    facet_wrap(~score, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = COL_CL, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.16)), labels = lab_num) +
    labs(x = NULL, y = "ESTIMATE score") + theme_fig + theme(panel.grid.major.x = element_blank())

  dir_scale <- scale_colour_manual(values = COL_CL, labels = c(C1 = "Higher in C1", C2 = "Higher in C2"),
                                   name = NULL, drop = FALSE)
  rank_panel <- function(dt, xlab, ci = FALSE, it = FALSE) {
    dt <- copy(dt)
    dt[, `:=`(sig = sig_from_p(p_value), dir = factor(fifelse(effect >= 0, "C1", "C2"), levels = c("C1", "C2")))]
    if (!"facet" %in% names(dt)) dt[, label := factor(label, levels = dt[order(abs(effect))]$label)]
    ends <- if (ci) c(dt$lo, dt$hi, dt$effect) else dt$effect
    rng <- max(abs(ends), na.rm = TRUE)
    dt[, star_x := (if (ci) fifelse(effect >= 0, pmax(hi, effect, na.rm = TRUE), pmin(lo, effect, na.rm = TRUE))
                    else effect) + sign(effect) * rng * 0.04]
    p <- ggplot(dt, aes(effect, label, colour = dir)) +
      geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3)
    # legend keys from the points only, so that the lollipop and forest panels share one legend
    p <- if (ci) p + geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.45, linewidth = 0.35, orientation = "y", show.legend = FALSE)
         else p + geom_segment(aes(x = 0, xend = effect, yend = label), linewidth = 0.45, show.legend = FALSE)
    p <- p + geom_point(size = 1.2) +
      geom_text(aes(x = star_x, label = sig, hjust = fifelse(effect >= 0, 0, 1)), colour = "black",
                size = 7.5 / .pt, vjust = 0.72, family = "Arial") +
      dir_scale + scale_x_continuous(expand = expansion(mult = 0.13), labels = lab_num) +
      labs(x = xlab, y = NULL) + theme_fig + theme(panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey92"))
    if (isTRUE(it)) p <- p + scale_y_discrete(labels = function(l) italic_labels(l, TRUE))
    else if (is.character(it)) p <- p + scale_y_discrete(labels = function(l) italic_labels(l, l %in% it))
    p
  }
  cib <- tb("Tab8B_CIBERSORT_celltype_stats.csv")
  cib[, effect := fifelse(median_C1 == 0 & median_C2 == 0, mean_C1 - mean_C2, median_C1 - median_C2)]
  pb <- rank_panel(cib[, .(label = cell_type, effect, p_value)], paste0("Median fraction difference (C1 ", MINUS, " C2)"))
  bin <- tb("Tab8C_ssGSEA_Bindea_celltype_stats.csv")
  pc <- rank_panel(bin[, .(label = cell_type, effect = median_C1 - median_C2, p_value)],
                   paste0("Median score difference (C1 ", MINUS, " C2)"))
  hl_ci <- function(sym) {                       # Hodges-Lehmann 95% CI of C1 - C2, as in v3_fig4_BCDE
    t(vapply(sym, function(g) {
      v <- expr_row(g, clu$sampleID)
      w <- wilcox.test(v[clu$cluster == "C1"], v[clu$cluster == "C2"], conf.int = TRUE, exact = FALSE)
      w$conf.int[1:2]
    }, numeric(2)))
  }
  cp <- tb("Tab8D_checkpoint_wilcoxon.csv")
  ci <- hl_ci(cp$gene_symbol)
  cpd <- cp[, .(label = gene_display, effect = log2FC_C1_minus_C2, lo = ci[, 1], hi = ci[, 2], p_value)]
  fc_lab <- bquote(log[2] ~ "fold change (C1" ~ .(MINUS) ~ "C2)")
  pd <- rank_panel(cpd, fc_lab, ci = TRUE, it = cp[gene_display == gene_symbol, gene_display])
  hla <- tb("Tab8E_HLA_wilcoxon.csv")
  ci <- hl_ci(hla$gene_symbol)
  hl <- hla[, .(label = gene_symbol, effect = log2FC_C1_minus_C2, lo = ci[, 1], hi = ci[, 2], p_value,
                facet = factor(fifelse(gene_symbol %in% paste0("HLA-", c("A", "B", "C", "E", "F", "G")),
                                       "Class I", "Class II"), levels = c("Class I", "Class II")))]
  hl <- hl[order(facet, abs(effect))]
  hl[, label := factor(label, levels = label)]
  pe <- rank_panel(hl, fc_lab, ci = TRUE, it = TRUE) +
    facet_grid(rows = vars(facet), scales = "free_y", space = "free_y", switch = "y") +
    theme(strip.placement = "outside", strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 0, size = 7.5, face = "bold"), panel.spacing.y = unit(1.5, "mm"))
  cat(sprintf("  stars (raw p) b %d/22, c %d/28, d %d/30, e %d/19\n", sum(cib$p_value < 0.05), sum(bin$p_value < 0.05),
              sum(cp$p_value < 0.05), sum(hla$p_value < 0.05)))
  # free(): (a) and (d) are not aligned with the long cell-type names of (b), so their panels use the width
  fig5 <- wrap_plots(A = free(pa), G = guide_area(), B = pb, C = pc, D = free(pd), F = pe,
                     design = "AA\nGG\nBC\nDF") +
    plot_layout(heights = c(2.6, 0.25, 29.2 / 4, 31.2 / 4), guides = "collect") + tags(5) &
    theme(legend.position = "top", legend.direction = "horizontal")
  save_both(fig5, "Figure5_draft", 21.6)
}

# =====================================================================================================  Figure 7
if ("7" %in% which_figs) {
  cat("Figure 7 (prognostic model)\n")
  cv <- readRDS(file.path(OUT, "lasso_rerun.rds"))
  fit <- cv$glmnet.fit
  gene_order <- rownames(fit$beta)                       # screening order = heat-map order of the published panel
  path <- data.table(gene = rep(gene_order, ncol(fit$beta)), x = rep(-log(fit$lambda), each = nrow(fit$beta)),
                     coef = as.vector(as.matrix(fit$beta)))
  # top axis as glmnet's plotCoef: value at each break by approx(method = "constant", f = 1); gives 2, 7, 7, 7 as published
  at <- function(v, lam, val) { y <- rep(NA_character_, length(v)); ok <- !is.na(v)
    y[ok] <- approx(-log(lam), val, xout = v[ok], method = "constant", f = 1, rule = 2)$y; y }
  gcol <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF"), gene_order)
  xr <- range(path$x)
  pA <- ggplot(path, aes(x, coef, colour = gene)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) +
    geom_line(linewidth = 0.5) +
    geom_text_repel(data = path[x == max(x)], aes(label = gene), fontface = "italic", size = 7 / .pt, family = "Arial",
                    direction = "y", hjust = 0, nudge_x = 0.12, xlim = c(xr[2] + 0.05, NA), segment.size = 0.2,
                    min.segment.length = 0, box.padding = 0.08, seed = 1262118388) +
    scale_colour_manual(values = gcol, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.30)),
                       sec.axis = dup_axis(name = "Number of nonzero coefficients", labels = function(b) at(b, fit$lambda, fit$df))) +
    scale_y_continuous(labels = lab_num) +
    labs(x = expression(-log(lambda)), y = "Coefficient", tag = "(a)") + theme_fig + theme(axis.title.x.top = element_text(size = 7))
  cvd <- data.table(x = -log(cv$lambda), m = cv$cvm, lo = cv$cvlo, hi = cv$cvup)
  # dotted lines at lambda.1se and lambda.min as in the published panel (unlabelled)
  vl <- data.table(x = -log(c(cv$lambda.1se, cv$lambda.min)))
  pB <- ggplot(cvd, aes(x, m)) +
    geom_vline(data = vl, aes(xintercept = x), linetype = "dotted", linewidth = 0.4) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.04, colour = "grey55", linewidth = 0.25) +
    geom_point(colour = "#E41A1C", size = 0.8) +
    scale_x_continuous(sec.axis = dup_axis(name = "Number of nonzero coefficients", labels = function(b) at(b, cv$lambda, cv$nzero))) +
    labs(x = expression(-log(lambda)), y = "Partial likelihood deviance", tag = "(b)") + theme_fig +
    theme(axis.title.x.top = element_text(size = 7))
  cat(sprintf("  (b) lambda.min -log %.3f, lambda.1se -log %.3f; min deviance %.3f\n", -log(cv$lambda.min),
              -log(cv$lambda.1se), min(cv$cvm)))

  rs <- tb("Tab10_risk_scores_per_sample.csv")
  kmr <- tb("Tab10_KM_logrank_results.csv")
  chk <- coxph(Surv(OS.time, OS) ~ factor(RiskGroup, levels = c("Low", "High")), data = rs)
  lr <- survdiff(Surv(OS.time, OS) ~ RiskGroup, data = rs)
  lr_p <- 1 - pchisq(lr$chisq, 1)
  stopifnot(abs(exp(coef(chk)) - kmr$HR) < 1e-8, abs(lr_p - kmr$logrank_p) < 1e-12)
  note7 <- list(
    annotate("text", x = 4000, y = 1, hjust = 1, vjust = 1, size = 7 / .pt, family = "Arial",
             label = sprintf("HR = %.2f (95%% CI %.2f%s%.2f)", kmr$HR, kmr$HR_low_95, NDASH, kmr$HR_high_95)),
    annotate("text", x = 4000, y = 0.9, hjust = 1, vjust = 1, size = 7 / .pt, family = "Arial", parse = TRUE,
             label = sprintf('"Log-rank" ~ italic(p) == %s', sci(kmr$logrank_p))))
  km7 <- km_panel(rs$OS.time, rs$OS, rs$RiskGroup, COL_RISK, "Risk group", c(Low = "Low risk", High = "High risk"), note7, "(c)")
  cat(sprintf("  (c) at risk: %s\n", km7$risk[, paste(n, collapse = "/"), by = group][, paste(group, V1, collapse = "; ")]))
  pC <- wrap_elements(full = km7$curve / km7$table + plot_layout(heights = c(3.4, 0.75)))

  times <- c(365, 730, 1095)
  roc <- timeROC(T = rs$OS.time, delta = rs$OS, marker = rs$RiskScore, cause = 1, times = times, iid = TRUE)
  auc_t <- tb("Tab10_timeROC_AUC.csv")
  se <- as.numeric(roc$inference$vect_sd_1)
  stopifnot(max(abs(roc$AUC - auc_t$AUC)) < 1e-10, max(abs(se - auc_t$SE)) < 1e-10)
  roc_lab <- sprintf("%d-year AUC %.3f (95%% CI %.3f%s%.3f)", 1:3, auc_t$AUC, auc_t$AUC_low_95, NDASH, auc_t$AUC_high_95)
  rd <- rbindlist(lapply(1:3, function(k) data.table(FPR = roc$FP[, k], TPR = roc$TP[, k], lab = roc_lab[k])))
  pD <- ggplot(rd, aes(FPR, TPR, colour = lab)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    geom_line(linewidth = 0.5) +
    scale_colour_manual(values = setNames(c("#377EB8", "#4DAF4A", "#E41A1C"), roc_lab), name = NULL) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "False positive rate", y = "True positive rate", tag = "(d)") + theme_fig +
    theme(legend.position = "inside", legend.position.inside = c(0.99, 0.02), legend.justification = c(1, 0),
          legend.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.25),
          legend.key.size = unit(2.6, "mm"), legend.margin = margin(1.5, 2, 1.5, 2))

  ord <- order(rs$RiskScore)
  N <- nrow(rs); med <- median(rs$RiskScore)
  r1 <- data.table(rank = seq_len(N), RS = rs$RiskScore[ord], grp = factor(rs$RiskGroup[ord], levels = c("Low", "High")),
                   time = rs$OS.time[ord], ev = factor(ifelse(rs$OS[ord] == 1, "Death", "Censored"), levels = c("Censored", "Death")))
  xs <- scale_x_continuous(expand = expansion(add = 1))
  bare_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
                  legend.justification = c(0, 0.5))
  e1 <- ggplot(r1, aes(rank, RS, colour = grp)) + geom_point(size = 0.55) +
    geom_hline(yintercept = med, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    annotate("text", x = 3, y = med, hjust = 0, vjust = -0.5, size = 7 / .pt, family = "Arial", colour = "grey30",
             label = paste0("Median = ", num(med, 4))) +
    scale_colour_manual(values = COL_RISK, name = "Risk group") + xs + scale_y_continuous(labels = lab_num) +
    labs(y = "Risk score", tag = "(e)") + theme_fig + bare_x
  e2 <- ggplot(r1, aes(rank, time, colour = ev)) + geom_point(size = 0.55, alpha = 0.85) +
    scale_colour_manual(values = c(Censored = "grey35", Death = "#E41A1C"), name = "OS event") + xs +
    labs(y = "Time (days)") + theme_fig + bare_x
  hz <- scale(sapply(gene_order, function(g) expr_row(g, rs$sampleID[ord])))   # samples x genes, z per gene
  hd <- data.table(rank = rep(seq_len(N), times = length(gene_order)), gene = factor(rep(gene_order, each = N), levels = rev(gene_order)),
                   z = as.vector(hz))
  e3 <- ggplot(hd, aes(rank, gene, fill = z)) + geom_raster() +
    scale_fill_gradient2(low = "#0571B0", mid = "white", high = "#CA0020", midpoint = 0, limits = c(-3, 3),
                         oob = scales::squish, name = "z-score") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0), labels = function(l) italic_labels(l, TRUE)) +
    labs(x = sprintf("Samples sorted by risk score (low \u2192 high, n = %d)", N), y = NULL) + theme_fig +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major = element_blank(),
          legend.justification = c(0, 0.5), legend.key.height = unit(3.5, "mm"))
  pE <- wrap_elements(full = e1 / e2 / e3 + plot_layout(heights = c(1, 1, 1.25)))
  fig7 <- (pA | pB) / (pC | pD) / pE + plot_layout(heights = c(4.6, 5.6, 5.8))     # tags set on the plots
  save_both(fig7, "Figure7_draft", 19.6)
}

# =====================================================================================================  Figure 4
if ("4" %in% which_figs) {
  cat("Figure 4 (consensus clustering)\n")
  cc <- readRDS(file.path(OUT, "ccp_rerun.rds"))
  dimnames(cc$cm2) <- list(cc$samples, cc$samples)     # CCP keeps the input column order
  stopifnot(identical(unname(cc$class2[clu$sampleID]), clu$cluster))
  cdf <- as.data.table(cc$cdf)
  pa <- ggplot(cdf, aes(x, y, colour = factor(k))) + geom_step(linewidth = 0.4) +
    scale_colour_brewer(palette = "Dark2", name = "k") +
    guides(colour = guide_legend(nrow = 2, title.position = "left", keywidth = unit(2.4, "mm"), keyheight = unit(2.4, "mm"))) +
    labs(x = "Consensus index", y = "CDF", tag = "(a)") + theme_fig +
    theme(legend.position = "inside", legend.position.inside = c(0.98, 0.03), legend.justification = c(1, 0),
          legend.background = element_blank(), legend.title = element_text(size = 7),
          legend.key.spacing.x = unit(0.6, "mm"), legend.text = element_text(size = 7, margin = margin(l = 0.8, unit = "mm")))
  mt <- data.table(k = rep(2:9, 2), v = c(cc$delta, cc$pac),
                   m = factor(rep(c("Delta area", "PAC"), each = 8), levels = c("Delta area", "PAC")))
  pb <- ggplot(mt, aes(k, v)) +
    geom_vline(xintercept = 2, linetype = "dashed", colour = COL_CL[["C2"]], linewidth = 0.35) +
    geom_line(colour = COL_CL[["C1"]], linewidth = 0.45) + geom_point(colour = COL_CL[["C1"]], size = 1.1) +
    facet_wrap(~m, scales = "free_y") + scale_x_continuous(breaks = 2:9) +
    labs(x = "k", y = NULL, tag = "(b)") + theme_fig
  cl_ord <- cc$samples[order(cc$class2[cc$samples])]        # scripts/07 STEP 10 order
  # (c) consensus matrix as in scripts/07 STEP 10 (white -> #08519C, cluster bar on top), drawn as a ggplot raster so
  # that its legend is laid out like the other panels; row 1 at the top as in the ComplexHeatmap original
  n <- length(cl_ord); n1 <- sum(cc$class2 == "C1")
  cmd <- data.table(x = rep(seq_len(n), each = n), y = rep(n:1, times = n), v = as.vector(cc$cm2[cl_ord, cl_ord]))
  pc <- ggplot(cmd, aes(x, y, fill = v)) + geom_raster() +
    annotate("rect", xmin = 0.5, xmax = n1 + 0.5, ymin = n + 3, ymax = n + 12, fill = COL_CL[["C1"]]) +
    annotate("rect", xmin = n1 + 0.5, xmax = n + 0.5, ymin = n + 3, ymax = n + 12, fill = COL_CL[["C2"]]) +
    scale_fill_gradient(low = "white", high = "#08519C", limits = c(0, 1), breaks = c(0, 0.5, 1), name = "Consensus") +
    coord_equal(expand = FALSE) + labs(x = NULL, y = NULL, tag = "(c)") + theme_fig +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(),
          panel.border = element_blank(), legend.key.width = unit(2.5, "mm"), legend.key.height = unit(4, "mm"))
  tag_theme <- theme(plot.tag = element_text(size = 10, face = "bold", family = "Arial"))

  # (d) PCA and (e) KM on the same z-scored 18-FMG matrix as scripts/07
  fmg <- tb("Tab1_2_FMG_list_primary.csv")$gene
  mz <- as.matrix(expr[match(fmg, gene_id), cc$samples, with = FALSE]); rownames(mz) <- fmg
  mz <- t(scale(t(mz)))
  pca <- prcomp(t(mz), center = FALSE, scale. = FALSE)
  ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  pdd <- data.table(PC1 = pca$x[, 1], PC2 = pca$x[, 2], cluster = cc$class2[cc$samples])
  pd <- ggplot(pdd, aes(PC1, PC2, colour = cluster)) + geom_point(size = 0.9, alpha = 0.85) +
    stat_ellipse(linewidth = 0.4) + scale_colour_manual(values = COL_CL, name = "Cluster") +
    scale_x_continuous(labels = lab_num) + scale_y_continuous(labels = lab_num) +
    labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]), tag = "(d)") + theme_fig +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0), legend.key = element_blank())
  lr4 <- survdiff(Surv(OS.time, OS) ~ cluster, data = clu)
  p4 <- 1 - pchisq(lr4$chisq, 1)
  note4 <- list(annotate("text", x = 4000, y = 1, hjust = 1, vjust = 1, size = 7 / .pt, family = "Arial", parse = TRUE,
                         label = sprintf('"Log-rank" ~ italic(p) == "%.2f"', p4)))
  km4 <- km_panel(clu$OS.time, clu$OS, clu$cluster, COL_CL, "Cluster", c(C1 = "C1", C2 = "C2"), note4, "(e)")
  cat(sprintf("  (e) log-rank p = %.4f; at risk: %s
", p4, km4$risk[, paste(n, collapse = "/"), by = group][, paste(group, V1, collapse = "; ")]))
  pe <- wrap_elements(full = km4$curve / km4$table + plot_layout(heights = c(3.4, 0.75)))

  # (f) FMG heat map, scripts/07 STEP 13 (rows clustered as there; gene names in italics)
  ord2 <- order(cc$class2[cc$samples])
  hist <- clu$histology[match(cc$samples[ord2], clu$sampleID)]
  HIST <- c(ESCC = "#7570B3", EAC = "#E7298A")
  # row clustering exactly as ComplexHeatmap did it for the published panel (same call as scripts/07 STEP 13);
  # the panel itself is a ggplot so that it is laid out like the others
  pdf(NULL)
  hd <- draw(Heatmap(mz[, ord2], cluster_columns = FALSE, cluster_rows = TRUE,
                     column_split = factor(cc$class2[cc$samples][ord2])))
  hc <- as.hclust(row_dend(hd)); invisible(dev.off())
  g_ord <- rownames(mz)[hc$order]
  stopifnot(identical(g_ord, c("FGFR2", "TP53", "RUNX1", "PIK3CA", "FAT3", "NOTCH3", "LRP1B", "ROS1", "FANCD2", "BRCA2",
                               "KMT2D", "ERCC6", "SLX4", "NOTCH1", "ASXL1", "ATM", "PKHD1", "NTRK1")))   # published order
  ng <- length(g_ord); ns <- length(ord2); n1 <- sum(cc$class2 == "C1"); GAP <- 2
  xs_f <- seq_len(ns) + ifelse(seq_len(ns) > n1, GAP, 0)          # column split between C1 and C2
  ypos <- setNames(ng:1, g_ord)                                    # first gene at the top
  hm <- data.table(x = rep(xs_f, each = ng), y = rep(ypos[rownames(mz)], times = ns), z = as.vector(mz[, ord2]))
  ann <- data.table(x = rep(xs_f, 2), y = rep(c(ng + 2.25, ng + 1.15), each = ns),
                    col = c(COL_CL[cc$class2[cc$samples][ord2]], c(ESCC = "#7570B3", EAC = "#E7298A")[hist]))
  # dendrogram segments (heights scaled into 9 columns' width left of the first column)
  pos <- numeric(0); nx <- ny <- numeric(ng - 1); seg <- list()
  leaf_y <- ypos[rownames(mz)]
  for (k in seq_len(ng - 1)) {
    m <- hc$merge[k, ]
    yy <- vapply(m, function(i) if (i < 0) leaf_y[[-i]] else ny[i], 1)
    hh <- vapply(m, function(i) if (i < 0) 0 else hc$height[i], 1)
    ny[k] <- mean(yy)
    seg[[k]] <- data.table(h0 = c(hh, hc$height[k]), h1 = rep(hc$height[k], 3), y0 = c(yy, yy[1]), y1 = c(yy, yy[2]))
  }
  seg <- rbindlist(seg)
  to_x <- function(h) -1.5 - 9 * h / max(hc$height)
  pf <- ggplot(hm, aes(x, y)) +
    geom_tile(aes(fill = z), width = 1, height = 1) +
    geom_tile(data = ann, aes(x, y, fill = NULL), fill = ann$col, height = 0.95, width = 1) +
    geom_segment(data = seg, aes(x = to_x(h0), xend = to_x(h1), y = y0, yend = y1), linewidth = 0.25, inherit.aes = FALSE) +
    geom_point(data = data.table(h = c("ESCC", "EAC")), aes(x = -Inf, y = -Inf, colour = h), shape = 15, size = 2.6,
               alpha = 0, inherit.aes = FALSE) +                      # legend only
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-2, 2),
                         oob = scales::squish, name = "z-score") +
    scale_colour_manual(values = c(ESCC = "#7570B3", EAC = "#E7298A"), breaks = c("EAC", "ESCC"), name = "Histology") +
    scale_x_continuous(expand = c(0, 0.5)) +
    scale_y_continuous(breaks = c(ypos, ng + 1.15, ng + 2.25), labels = italic_labels(c(names(ypos), "Histology", "Cluster"),
                       c(rep(TRUE, ng), FALSE, FALSE)), position = "right", expand = c(0, 0.3)) +
    labs(x = NULL, y = NULL, tag = "(f)") + theme_fig +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_blank(), legend.key.width = unit(2.5, "mm"), legend.key.height = unit(3.5, "mm"),
          legend.box.spacing = unit(1, "mm"))

  # (g) FMG expression by cluster, figures/remake/v3_fig3G_violin.R (raw-p stars)
  gord <- tb("Tab4_TCGA_mutation_frequency.csv")$gene
  gl <- rbindlist(lapply(gord, function(g) data.table(gene = g, value = expr_row(g, clu$sampleID), cluster = clu$cluster)))
  gl[, gene := factor(gene, levels = gord)]
  pg_dt <- gl[, .(p = wilcox.test(value ~ cluster, exact = FALSE)$p.value), by = gene]
  t7c <- tb("Tab7C_FMG_expr_diff_clusters_C1vsC2.csv")
  stopifnot(max(abs(pg_dt$p - t7c$p_value[match(pg_dt$gene, t7c$gene)])) < 1e-12)
  pg_dt[, `:=`(sig = sig_from_p(p), y = max(gl$value) + diff(range(gl$value)) * 0.03)]
  cat(sprintf("  (g) stars on %d of 18 genes\n", sum(pg_dt$sig != "")))
  pg <- ggplot(gl, aes(gene, value, fill = cluster)) +
    geom_boxplot(width = 0.7, colour = "black", outlier.size = 0.25, linewidth = 0.25, position = position_dodge(0.8)) +
    geom_text(data = pg_dt, aes(gene, y, label = sig), inherit.aes = FALSE, size = 7.5 / .pt, vjust = 0.6, family = "Arial") +
    scale_fill_manual(values = COL_CL, name = "Cluster") +
    scale_x_discrete(labels = function(l) italic_labels(l, TRUE)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = NULL, y = expression("Expression (log"[2] * "(RSEM + 1))"), tag = "(g)") + theme_fig +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.x = element_blank(),
          legend.position = "top", legend.margin = margin(0, 0, -4, 0))

  # rows nested, so (c) gets the width its legend needs without moving the panels of rows 2-4; tags set on the plots
  fig4 <- ((pa | pb | pc) + plot_layout(widths = c(1.3, 1.55, 0.85))) / ((pd | pe) + plot_layout(widths = c(1, 1.3))) /
    pf / free(pg, side = "lr") + plot_layout(heights = c(3.4, 4.0, 5.6, 2.9))     # free: (g) is not narrowed by the legends of (f)
  save_both(fig4, "Figure4_draft", 21.6)
}

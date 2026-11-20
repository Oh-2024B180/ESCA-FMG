###############################################################################
##  Stage 13 helper — Fig9A nomogram + Fig9B calibration cosmetic patch
##
##  Input (read-only):
##    tables/Tab10_risk_scores_per_sample.csv
##    data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv
##
##  Output (overwrite):
##    figures/Fig9A_nomogram.{pdf,png}
##    figures/Fig9B_calibration.{pdf,png}
##    figures/Fig9C_C_index_comparison.{pdf,png}

# run from the repository root (see README.md)

suppressPackageStartupMessages({
  library(data.table)
  library(rms)
  library(Hmisc)
  library(survival)
  library(ggplot2)
  library(ragg)
})

set.seed(1262118388)
t0 <- Sys.time()

INPUT_RS  <- "tables/Tab10_risk_scores_per_sample.csv"
INPUT_CLI <- "data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv"
FIGDIR    <- "figures"

CALIB_B      <- 1000
AT_RISK_MIN  <- 10

# ── 1. Load + same collapse logic as 13_nomogram_construction.R ──────────────────────────
cat("[1] Load + build ph_nom (complete-case for RS+histology+stage+residual)\n")
rs_dt  <- fread(INPUT_RS)
cli_dt <- fread(INPUT_CLI)

collapse_stage <- function(s) {
  s[s == "" | grepl("^\\[", s)] <- NA
  out <- rep(NA_character_, length(s))
  out[grepl("^Stage I[ABC]?$|^Stage I$", s)] <- "I"
  out[grepl("^Stage II", s)]  <- "II"
  out[grepl("^Stage III", s)] <- "III"
  out[grepl("^Stage IV", s)]  <- "IV"
  out
}
collapse_residual <- function(r) {
  r[r == "" | grepl("^\\[", r)] <- NA
  out <- rep(NA_character_, length(r))
  out[r == "R0"] <- "R0"
  out[r %in% c("R1","R2")] <- "R1+R2"
  out
}
collapse_histology <- function(h) {
  ifelse(h == "Esophagus Squamous Cell Carcinoma", "ESCC",
  ifelse(h == "Esophagus Adenocarcinoma, NOS",     "EAC", NA))
}

cli_dt[, stage_collapsed     := collapse_stage(pathologic_stage)]
cli_dt[, residual_collapsed  := collapse_residual(residual_tumor)]
cli_dt[, histology_collapsed := collapse_histology(histological_type)]

ph <- merge(rs_dt[, .(sampleID, OS, OS.time, RiskScore)],
            cli_dt[, .(sampleID, stage_collapsed, residual_collapsed,
                       histology_collapsed)],
            by = "sampleID")

ph_nom <- ph[complete.cases(ph[, .(RiskScore, stage_collapsed,
                                    residual_collapsed, histology_collapsed)])]
ph_nom[, stage_collapsed     := factor(stage_collapsed,     levels = c("I","II","III","IV"))]
ph_nom[, residual_collapsed  := factor(residual_collapsed,  levels = c("R0","R1+R2"))]
ph_nom[, histology_collapsed := factor(histology_collapsed, levels = c("ESCC","EAC"))]

ph_nom_df <- as.data.frame(ph_nom)
N_NOM <- nrow(ph_nom_df)
cat(sprintf("    N (nomogram CC): %d\n", N_NOM))

# ── 2. Hmisc::label() → friendly variable names on nomogram ────────────────
cat("[2] Apply Hmisc::label() friendly names\n")
label(ph_nom_df$RiskScore)            <- "Risk score"
label(ph_nom_df$histology_collapsed)  <- "Histology"
label(ph_nom_df$stage_collapsed)      <- "Pathologic stage"
label(ph_nom_df$residual_collapsed)   <- "Residual tumor"

dd <- datadist(ph_nom_df); options(datadist = "dd")

nom_fit <- cph(Surv(OS.time, OS) ~ RiskScore + histology_collapsed +
                 stage_collapsed + residual_collapsed,
               data = ph_nom_df, x = TRUE, y = TRUE, surv = TRUE,
               time.inc = 365)

# ── 3. Fig9A nomogram ─────────────────────────────────────────────────────
cat("[3] Render Fig9A nomogram (with friendly labels)\n")
surv_obj <- Survival(nom_fit)
surv_1yr <- function(x) surv_obj(365,  x)
surv_2yr <- function(x) surv_obj(730,  x)
surv_3yr <- function(x) surv_obj(1095, x)

nom <- nomogram(nom_fit,
                fun = list(surv_1yr, surv_2yr, surv_3yr),
                funlabel = c("1-year survival", "2-year survival",
                             "3-year survival"),
                lp = TRUE,
                fun.at = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 0.95))

W <- 11; H <- 8
draw_nom <- function() {
  par(mar = c(3, 4, 1, 1))
  plot(nom, xfrac = 0.30, cex.axis = 0.85, cex.var = 0.95, lmgp = 0.2)
}
cairo_pdf(file.path(FIGDIR, "Fig9A_nomogram.pdf"), width = W, height = H)
draw_nom(); dev.off()
ragg::agg_png(file.path(FIGDIR, "Fig9A_nomogram.png"), width = W, height = H,
              units = "in", res = 600)
draw_nom(); dev.off()
cat("    [OK] Fig9A_nomogram.{pdf,png}\n")

# ── 4. Fig9B calibration ─────────────────────────────
cat("[4] Render Fig9B calibration (1-yr + 3-yr with B=1000 bootstrap)\n")
at_risk_5yr <- sum(ph_nom_df$OS.time >= 1825)
cat(sprintf("    at-risk @ 1-yr: %d, 2-yr: %d, 3-yr: %d, 5-yr: %d\n",
            sum(ph_nom_df$OS.time >= 365), sum(ph_nom_df$OS.time >= 730),
            sum(ph_nom_df$OS.time >= 1095), at_risk_5yr))
cat(sprintf("    5-yr fallback rule: skip 5-yr if at-risk < %d\n", AT_RISK_MIN))
if (at_risk_5yr < AT_RISK_MIN) {
  times_calib  <- c(365, 730, 1095)
  labels_calib <- c("1-year","2-year","3-year")
} else {
  times_calib  <- c(365, 730, 1095, 1825)
  labels_calib <- c("1-year","2-year","3-year","5-year")
}

draw_calib <- function() {
  par(mfrow = c(1, length(times_calib)),
      mar  = c(5, 4.5, 4, 1.2),
      oma  = c(0.5, 0, 2.5, 0))
  for (i in seq_along(times_calib)) {
    t_d <- times_calib[i]
    fit_t <- cph(Surv(OS.time, OS) ~ RiskScore + histology_collapsed +
                   stage_collapsed + residual_collapsed,
                 data = ph_nom_df, x = TRUE, y = TRUE, surv = TRUE,
                 time.inc = t_d)
    cal <- calibrate(fit_t, cmethod = "KM", method = "boot", u = t_d, m = 50,
                     B = CALIB_B)
    # subtitles=FALSE: suppresses rms default "n=… d=… p=… 50 subjects per group",
    # "X - resampling optimism added, B=…", "Gray: ideal", "Based on observed-predicted"
    plot(cal,
         xlim = c(0, 1), ylim = c(0, 1),
         xlab = "Nomogram predicted survival",
         ylab = "Observed survival",
         # Per-panel main label retained as a sub-panel identifier (1-yr/3-yr).
         main = labels_calib[i],
         riskdist  = FALSE,
         subtitles = FALSE)
    abline(0, 1, lty = 3, col = "grey50")
  }
}

W2 <- 6 * length(times_calib); H2 <- 6
cairo_pdf(file.path(FIGDIR, "Fig9B_calibration.pdf"), width = W2, height = H2)
draw_calib(); dev.off()
ragg::agg_png(file.path(FIGDIR, "Fig9B_calibration.png"), width = W2, height = H2,
              units = "in", res = 600)
draw_calib(); dev.off()
cat("    [OK] Fig9B_calibration.{pdf,png}\n")

# ── 5. Fig9C C-index bar ─────────────────────
cat("[5] Render Fig9C C-index comparison (anchor labels at c_high_95)\n")
INPUT_CIDX <- "tables/Tab13_C_index_comparison.csv"
cidx_dt <- fread(INPUT_CIDX)
cidx_dt[, model := factor(model, levels = model)]

p9C <- ggplot(cidx_dt, aes(x = model, y = c_index, fill = model)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = c_low_95, ymax = c_high_95),
                width = 0.2, color = "black", linewidth = 0.4) +
  geom_text(aes(y = c_high_95, label = sprintf("%.3f", c_index)),
            vjust = -0.8, size = 4, fontface = "bold") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  scale_y_continuous(limits = c(0, 1.05),
                     breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0, 0.04))) +
  scale_fill_manual(values = c("#E41A1C","#377EB8","#4DAF4A"), guide = "none") +
  labs(x = NULL, y = "Concordance index") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(face = "bold", size = 8),
        plot.margin = margin(8, 8, 8, 8))

ggsave(file.path(FIGDIR, "Fig9C_C_index_comparison.pdf"),
       plot = p9C, width = 9, height = 5, device = cairo_pdf)
ggsave(file.path(FIGDIR, "Fig9C_C_index_comparison.png"),
       plot = p9C, width = 9, height = 5, dpi = 600, device = ragg::agg_png)
cat("    [OK] Fig9C_C_index_comparison.{pdf,png}\n")

elapsed <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("\n=== Fig9A/9B/9C re-render done in %.1f sec ===\n", elapsed))
cat("Tab13_*.csv NOT touched.\n")

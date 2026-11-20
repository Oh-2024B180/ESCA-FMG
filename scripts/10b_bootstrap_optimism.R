###############################################################################
## Stage 10b: Bootstrap optimism and selection stability of the seven-gene risk score
##
## Manuscript: Methods "Prognostic risk model construction" (last two sentences),
##             Results "Construction of the prognostic risk model" (last paragraph), Table 4 footnote a.
##
## Procedure (Harrell 1996; Meinshausen & Buhlmann 2010):
##   B = 1,000 resamples drawn with replacement from the 181 TCGA-ESCA patients. In each resample the
##   whole selection pipeline of Stage 10 is repeated: univariate Cox screen (Wald p < 0.01) over the
##   2,304 protein-coding DEG, then LASSO-Cox (glmnet, alpha = 1, 10-fold CV, lambda.min).
##   optimism = mean( C of the resample model in its own resample - C of the same model in the original cohort )
##   optimism-corrected C = apparent C of the final model - optimism
##   C is always taken for the fixed score with "higher score = shorter survival"
##   (concordance(..., reverse = TRUE)); a model is never refitted in the data it is evaluated on.
##
## Inputs  : tables/Tab10_risk_scores_per_sample.csv, tables/Tab10_univariate_Cox_results.csv (gene list),
##           data/expr_ESCA_tumor.csv
## Outputs : cache/boot_optimism_B1000.rds            (one list element per resample)
##           tables/Tab10b_bootstrap_optimism.csv     (summary used in the manuscript)
##
## History : the script of the first run (optimism 0.147, corrected C 0.561, re-selection 18-45%, median
##           23 genes, IQR 19-27) was not kept. This file was written on 2026-10-06 from the Methods text and
##           gives 0.149 / 0.560 / 18-45% / 23 (19-27); the small difference is the random-number stream
##           (Monte Carlo SE of the optimism 0.001). Since the revision the manuscript reports these values.
## Runtime : about 36 min with 12 workers (cv.glmnet takes ~30 s per resample).
##           Results depend on the seed AND on the number of workers (one random stream per worker).
## RERUN   : FALSE reads the cached resamples and only rebuilds the summary table.
###############################################################################
suppressPackageStartupMessages({library(survival); library(glmnet); library(data.table); library(parallel)})
ROOT <- ""   # repository root = working directory (see README.md)
B <- 1000; N_WORKERS <- 12; SEED <- 1262118388; RERUN <- FALSE
CACHE_RDS <- paste0(ROOT, "cache/boot_optimism_B1000.rds")
OUT_CSV <- paste0(ROOT, "tables/Tab10b_bootstrap_optimism.csv")
SIG <- c("FOLR3", "ZBTB7C", "RPRML", "HEPACAM2", "CLCA1", "CACNG8", "LIPC")

rs <- fread(paste0(ROOT, "tables/Tab10_risk_scores_per_sample.csv"))
uni <- fread(paste0(ROOT, "tables/Tab10_univariate_Cox_results.csv"))
ex <- fread(paste0(ROOT, "data/expr_ESCA_tumor.csv")); setnames(ex, 1, "gene")
ex <- ex[gene %in% uni$gene]
X <- t(as.matrix(ex[, rs$sampleID, with = FALSE])); colnames(X) <- ex$gene
time <- rs$OS.time; event <- rs$OS
cat(sprintf("matrix: %d patients x %d genes, %d events\n", nrow(X), ncol(X), sum(event)))

# univariate Cox screen: Wald p value of each gene (same test as summary(coxph()) in Stage 10)
screen <- function(Xb, tb, eb) {
  y <- Surv(tb, eb); ctl <- coxph.control(); p <- rep(NA_real_, ncol(Xb))
  for (j in seq_len(ncol(Xb))) {
    x <- Xb[, j, drop = FALSE]
    if (sd(x) == 0) next
    f <- tryCatch(coxph.fit(x, y, strata = NULL, offset = NULL, init = NULL, control = ctl, method = "efron", rownames = NULL),
                  error = function(e) NULL, warning = function(w) NULL)
    if (is.null(f) || is.na(f$coefficients) || f$var[1] <= 0) next
    p[j] <- 2 * pnorm(-abs(f$coefficients / sqrt(f$var[1])))
  }
  p
}
# concordance of a fixed score, higher score = shorter survival
cidx <- function(score, tt, ee) concordance(Surv(tt, ee) ~ score, reverse = TRUE)$concordance

# the original cohort must give the published model
p0 <- screen(X, time, event); g0 <- colnames(X)[which(p0 < 0.01)]
stopifnot(setequal(g0, SIG))
set.seed(SEED)
cv0 <- cv.glmnet(X[, g0, drop = FALSE], Surv(time, event), family = "cox", alpha = 1, nfolds = 10, standardize = TRUE, type.measure = "deviance")
b0 <- as.numeric(coef(cv0, s = "lambda.min")); names(b0) <- g0
stopifnot(all(b0 != 0))
app0 <- cidx(as.numeric(X[, g0, drop = FALSE] %*% b0), time, event)
cat(sprintf("final model reproduced (%s); apparent C = %.4f\n", paste(sprintf("%s %+.4f", g0, b0), collapse = "; "), app0))

one <- function(b) {
  idx <- sample.int(nrow(X), replace = TRUE)
  Xb <- X[idx, , drop = FALSE]; tb <- time[idx]; eb <- event[idx]
  p <- screen(Xb, tb, eb); pass <- which(p < 0.01)
  if (length(pass) < 2) return(list(ok = FALSE, npass = length(pass)))
  cv <- tryCatch(cv.glmnet(Xb[, pass, drop = FALSE], Surv(tb, eb), family = "cox", alpha = 1, nfolds = 10, standardize = TRUE, type.measure = "deviance"),
                 error = function(e) NULL)
  if (is.null(cv)) return(list(ok = FALSE, npass = length(pass)))
  be <- as.numeric(coef(cv, s = "lambda.min")); sel <- which(be != 0)
  if (!length(sel)) return(list(ok = FALSE, npass = length(pass)))
  g <- colnames(X)[pass[sel]]; be <- be[sel]
  list(ok = TRUE, npass = length(pass), nsel = length(g), genes = g,
       app = cidx(as.numeric(Xb[, g, drop = FALSE] %*% be), tb, eb),
       test = cidx(as.numeric(X[, g, drop = FALSE] %*% be), time, event))
}

if (RERUN || !file.exists(CACHE_RDS)) {
  cl <- makeCluster(N_WORKERS)
  invisible(clusterEvalQ(cl, {suppressPackageStartupMessages({library(survival); library(glmnet)})}))
  clusterExport(cl, c("X", "time", "event", "screen", "cidx", "one"))
  clusterSetRNGStream(cl, SEED)
  t0 <- Sys.time(); res <- parLapply(cl, seq_len(B), one); stopCluster(cl)
  cat(sprintf("B = %d done in %.1f min\n", B, as.numeric(Sys.time() - t0, units = "mins")))
  saveRDS(res, CACHE_RDS)
} else {
  res <- readRDS(CACHE_RDS); cat("read cached resamples:", CACHE_RDS, "\n")
}

ok <- Filter(function(r) isTRUE(r$ok), res)
app <- sapply(ok, `[[`, "app"); test <- sapply(ok, `[[`, "test"); ns <- sapply(ok, `[[`, "nsel")
opt <- mean(app - test)
fr <- sapply(SIG, function(g) 100 * mean(sapply(ok, function(r) g %in% r$genes)))
out <- rbind(
  data.table(item = "resamples drawn", value = length(res)),
  data.table(item = "resamples with a model", value = length(ok)),
  data.table(item = "apparent C of the final model", value = app0),
  data.table(item = "mean C of resample models in their own resample", value = mean(app)),
  data.table(item = "mean C of resample models in the original cohort", value = mean(test)),
  data.table(item = "optimism", value = opt),
  data.table(item = "Monte Carlo SE of the optimism", value = sd(app - test) / sqrt(length(ok))),
  data.table(item = "optimism-corrected C", value = app0 - opt),
  data.table(item = "genes per resampled signature, median", value = median(ns)),
  data.table(item = "genes per resampled signature, first quartile", value = unname(quantile(ns, .25))),
  data.table(item = "genes per resampled signature, third quartile", value = unname(quantile(ns, .75))),
  data.table(item = paste0("re-selection frequency (%), ", SIG), value = fr))
fwrite(out, OUT_CSV)
print(out[, .(item, value = round(value, 4))])
cat("written:", OUT_CSV, "\n")

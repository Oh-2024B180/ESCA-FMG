# Stage 14b (2026-10-19): does the oncoPredict result depend on feeding GDSC2_Res (ln IC50) to calcPhenotype()
# without exp()? The package vignette exponentiates first ("the IC50s are actual ln values ... calcPhenotype will do a
# power transformation ... better to not have both transformations"). Stage 14C did not.
#   run A = stage 14C settings (ln IC50 in)            -> must reproduce tables/Tab14C_oncoPredict_IC50_full.csv
#   run B = vignette (exp(ln IC50) in, i.e. IC50)      -> predictions back on the IC50 scale; compared on ln scale
# Writes A.rds and B.rds into OUT and run B as tables/Tab14Cb_oncoPredict_IC50_vignette.csv (read by stage 14c).
suppressPackageStartupMessages({ library(data.table); library(oncoPredict) })
E <- ""   # repository root = working directory (see README.md)
OUT <- "cache/oncopredict_runs/"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("R %s, oncoPredict %s\n", paste(R.Version()[c("major","minor")], collapse = "."), packageVersion("oncoPredict")))

gdsc_expr <- readRDS(paste0(E, "data/GDSC2_Expr (RMA Normalized and Log Transformed).rds"))
gdsc_res <- readRDS(paste0(E, "data/GDSC2_Res.rds"))
rs <- fread(paste0(E, "tables/Tab10_risk_scores_per_sample.csv"))
expr_dt <- fread(paste0(E, "data/expr_ESCA_tumor.csv"))
genes <- expr_dt[[1]]
expr_mat <- as.matrix(expr_dt[, -1, with = FALSE]); mode(expr_mat) <- "numeric"; rownames(expr_mat) <- genes
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]
idx <- match(rs$sampleID, colnames(expr_mat)); idx <- idx[!is.na(idx)]
expr_test <- expr_mat[, idx, drop = FALSE]
cat(sprintf("test %d x %d, training expr %d x %d, IC50 %d x %d (range %.2f to %.2f)\n", nrow(expr_test), ncol(expr_test),
            nrow(gdsc_expr), ncol(gdsc_expr), nrow(gdsc_res), ncol(gdsc_res), min(gdsc_res, na.rm = TRUE), max(gdsc_res, na.rm = TRUE)))

run <- function(ptype, tag) {
  d <- paste0(OUT, tag); dir.create(d, showWarnings = FALSE)
  old <- setwd(d); on.exit(setwd(old))
  t0 <- Sys.time()
  calcPhenotype(trainingExprData = gdsc_expr, trainingPtype = ptype, testExprData = expr_test,
                batchCorrect = "eb", powerTransformPhenotype = TRUE, removeLowVaryingGenes = 0.4,
                minNumSamples = 10, printOutput = FALSE, removeLowVaringGenesFrom = "rawData")
  cat(sprintf("run %s done in %.1f min\n", tag, as.numeric(Sys.time() - t0, units = "mins")))
  p <- fread(file.path(d, "calcPhenotype_Output", "DrugPredictions.csv")); setnames(p, 1, "sampleID"); p
}
A <- run(gdsc_res, "A_ln")
orig <- fread(paste0(E, "tables/Tab14C_oncoPredict_IC50_full.csv"))
drugs <- setdiff(colnames(orig), "sampleID")
stopifnot(all(drugs %in% colnames(A)))
a <- as.matrix(A[match(orig$sampleID, A$sampleID), ..drugs]); o <- as.matrix(orig[, ..drugs])
cat(sprintf("A vs published Tab14C: max |difference| = %.3g\n", max(abs(a - o), na.rm = TRUE)))
saveRDS(A, paste0(OUT, "A.rds"))

B <- run(exp(gdsc_res), "B_exp")
saveRDS(B, paste0(OUT, "B.rds"))
fwrite(B[match(rs$sampleID, B$sampleID)], paste0(E, "tables/Tab14Cb_oncoPredict_IC50_vignette.csv"))   # read by 14c
cat("=== done ===\n")

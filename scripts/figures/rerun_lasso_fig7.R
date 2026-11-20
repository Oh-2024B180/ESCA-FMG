###############################################################################
## Figure 7 (a)-(b): the cv.glmnet object behind the published panels was not kept (the files in
## cache/ are from another run), so the LASSO-Cox step of scripts/10_LASSO_Cox_model.R is run again on the same
## 7 screened genes (tables/Tab10_univariate_Cox_p001_genes.csv), samples and seed. The result is kept ONLY if
## lambda.min selects the same genes with the same coefficients as tables/Tab10_LASSO_selected_genes.csv and the
## risk score equals tables/Tab10_risk_scores_per_sample.csv.
##   Rscript rerun_lasso_fig7.R     (R 4.4.1)
## Output: figures/lasso_rerun.rds (cv.glmnet object) + lasso_rerun.log
###############################################################################
suppressPackageStartupMessages({ library(data.table); library(survival); library(glmnet) })
E <- getwd()   # repository root (see README.md)
OUT <- file.path(E, "figures", "final")
log <- file(file.path(OUT, "lasso_rerun.log"), open = "wt"); sink(log, split = TRUE)
cat(sprintf("LASSO re-run %s | R %s | glmnet %s | survival %s\n", Sys.time(), getRversion(),
            packageVersion("glmnet"), packageVersion("survival")))
expr <- fread(file.path(E, "data/expr_ESCA_tumor.csv")); setnames(expr, 1, "gene_id")
surv <- fread(file.path(E, "data/survival_ESCA_survival.csv")); setnames(surv, "sample", "sampleID")
cli  <- fread(file.path(E, "data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv"))
s3 <- Reduce(intersect, list(setdiff(colnames(expr), "gene_id"), cli[sample_type == "Primary Tumor", sampleID], surv$sampleID))
keep <- surv[sampleID %in% s3][order(match(sampleID, s3))][!is.na(OS) & !is.na(OS.time) & OS.time > 0]
screened <- fread(file.path(E, "tables/Tab10_univariate_Cox_p001_genes.csv"))$gene
x <- t(as.matrix(expr[match(screened, expr$gene_id), keep$sampleID, with = FALSE])); colnames(x) <- screened
cat(sprintf("%d samples, %d events, genes: %s\n", nrow(x), sum(keep$OS), paste(screened, collapse = ", ")))
set.seed(1262118388)
cv_fit <- cv.glmnet(x, Surv(keep$OS.time, keep$OS), family = "cox", alpha = 1, nfolds = 10, standardize = TRUE,
                    type.measure = "deviance")
b <- as.numeric(coef(cv_fit, s = "lambda.min")); names(b) <- screened
saved <- fread(file.path(E, "tables/Tab10_LASSO_selected_genes.csv"))
d_coef <- max(abs(b[saved$gene] - saved$coefficient))
rs <- as.numeric(x[, names(b)[b != 0]] %*% b[b != 0])
rs_saved <- fread(file.path(E, "tables/Tab10_risk_scores_per_sample.csv"))
d_rs <- max(abs(rs - rs_saved$RiskScore[match(keep$sampleID, rs_saved$sampleID)]))
cat(sprintf("lambda.min %.6g (-log %.4f), lambda.1se %.6g (-log %.4f)\n", cv_fit$lambda.min, -log(cv_fit$lambda.min),
            cv_fit$lambda.1se, -log(cv_fit$lambda.1se)))
cat(sprintf("selected %d genes: %s\nmax |coef - saved| = %.3g, max |risk score - saved| = %.3g\n", sum(b != 0),
            paste(names(b)[b != 0], collapse = ", "), d_coef, d_rs))
stopifnot(setequal(names(b)[b != 0], saved$gene), d_coef < 1e-6, d_rs < 1e-6)
saveRDS(cv_fit, file.path(OUT, "lasso_rerun.rds"))
cat("saved lasso_rerun.rds\n"); sink(); close(log)

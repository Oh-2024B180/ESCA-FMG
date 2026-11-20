###############################################################################
## Figure 4 (a)-(c): the CCP numbers were never saved (only the panel images), so the consensus
## clustering of scripts/07_consensus_clustering.R (STEP 2-6) is run again with the same inputs, parameters and
## seed, and its result is kept ONLY if the k = 2 assignment equals tables/Tab7_cluster_assignment.csv exactly.
##   Rscript rerun_ccp_fig4.R      (R 4.4.1, ConsensusClusterPlus 1.68.0; the analysis used 1.70.0)
## Output: figures/ccp_rerun.rds (CDF per k, PAC, delta area, k = 2 consensus matrix and classes) + ccp_rerun.log
###############################################################################
suppressPackageStartupMessages({ library(data.table); library(ConsensusClusterPlus) })
E <- getwd()   # repository root (see README.md)
OUT <- file.path(E, "figures", "final")
SEED <- 1262118388
log <- file(file.path(OUT, "ccp_rerun.log"), open = "wt"); sink(log, split = TRUE)
cat(sprintf("CCP re-run %s | R %s | ConsensusClusterPlus %s | cluster %s\n", Sys.time(), getRversion(),
            packageVersion("ConsensusClusterPlus"), packageVersion("cluster")))
set.seed(SEED)
expr <- fread(file.path(E, "data/expr_ESCA_tumor.csv"))
fmg  <- fread(file.path(E, "tables/Tab1_2_FMG_list_primary.csv"))
surv <- fread(file.path(E, "data/survival_ESCA_survival.csv"))
cli  <- fread(file.path(E, "data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv"))
samples_keep <- Reduce(intersect, list(setdiff(colnames(expr), "sample"),
                                       cli[sample_type == "Primary Tumor", sampleID], surv$sample))
expr_fmg <- expr[match(fmg$gene, expr$sample), c("sample", samples_keep), with = FALSE]
mat_raw <- as.matrix(expr_fmg[, -1]); rownames(mat_raw) <- expr_fmg$sample; storage.mode(mat_raw) <- "numeric"
mat_z <- t(scale(t(mat_raw), center = TRUE, scale = TRUE))
cat(sprintf("matrix %d genes x %d samples\n", nrow(mat_z), ncol(mat_z)))
tmp <- file.path(tempdir(), "ccp_native")
ccp <- ConsensusClusterPlus(d = mat_z, maxK = 9, reps = 1000, pItem = 0.8, pFeature = 1.0, clusterAlg = "pam",
                            distance = "spearman", innerLinkage = "ward.D2", finalLinkage = "ward.D2", seed = SEED,
                            plot = "pdf", writeTable = FALSE, title = tmp, verbose = FALSE)
compute_cdf <- function(M) { v <- M[lower.tri(M, diag = FALSE)]; data.frame(x = sort(v), y = ecdf(v)(sort(v))) }
compute_pac <- function(M) { v <- M[lower.tri(M, diag = FALSE)]; e <- ecdf(v); e(0.9) - e(0.1) }
ks <- 2:9
cdf_df <- do.call(rbind, lapply(ks, function(k) { d <- compute_cdf(ccp[[k]]$consensusMatrix); d$k <- k; d }))
pac_vec <- setNames(vapply(ks, function(k) compute_pac(ccp[[k]]$consensusMatrix), 1), ks)
auc_vec <- vapply(ks, function(k) { d <- compute_cdf(ccp[[k]]$consensusMatrix); sum(diff(c(0, d$x, 1)) * c(0, d$y)) }, 1)
delta_auc <- setNames(c(auc_vec[1], diff(auc_vec) / auc_vec[-length(auc_vec)]), ks)
cat("PAC by k:\n"); print(round(pac_vec, 4)); cat("Delta area by k:\n"); print(round(delta_auc, 4))
cl <- setNames(paste0("C", ccp[[2]]$consensusClass), colnames(mat_z))
saved <- fread(file.path(E, "tables/Tab7_cluster_assignment.csv"))
same <- identical(unname(cl[saved$sampleID]), saved$cluster)
cat(sprintf("k = 2 sizes: %s | identical to Tab7_cluster_assignment.csv: %s\n",
            paste(names(table(cl)), table(cl), sep = "=", collapse = ", "), same))
stopifnot(same)
saveRDS(list(cdf = cdf_df, pac = pac_vec, delta = delta_auc, cm2 = ccp[[2]]$consensusMatrix, class2 = cl,
             samples = colnames(mat_z), versions = c(R = as.character(getRversion()),
             CCP = as.character(packageVersion("ConsensusClusterPlus")), cluster = as.character(packageVersion("cluster")))),
        file.path(OUT, "ccp_rerun.rds"))
cat("saved ccp_rerun.rds\n"); sink(); close(log)

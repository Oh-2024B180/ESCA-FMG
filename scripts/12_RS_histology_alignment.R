###############################################################################
## Data-driven Discussion material:
##  A. Is the risk score itself aligned with cluster (C1/C2) and histology (ESCC/EAC)?
##  B. Which of the 7 signature genes drive that alignment (per-gene C1-vs-C2 and
##     ESCC-vs-EAC direction), and does each gene's contribution to RS push toward
##     one histology?
##  C. How well does the RS discriminate histology (AUC)?  (contrast with survival)
##  D. CODA vs TCGA 18-FMG prevalence correlation (population divergence).
###############################################################################
# run from the repository root (see README.md)
suppressPackageStartupMessages({library(data.table)})

expr <- fread("data/expr_ESCA_tumor.csv"); setnames(expr, 1, "gene")
clu  <- fread("tables/Tab7_cluster_assignment.csv")   # sampleID, cluster, histology, OS, OS.time
rs   <- fread("tables/Tab10_risk_scores_per_sample.csv")

GENES <- c("FOLR3","ZBTB7C","RPRML","HEPACAM2","CLCA1","CACNG8","LIPC")
BETA  <- c(FOLR3=0.2216, ZBTB7C=-0.1222, RPRML=0.1203, HEPACAM2=0.0757,
           CLCA1=0.0314, CACNG8=0.0314, LIPC=0.0301)

## histology label in cluster table
clu[, hist := toupper(as.character(histology))]
clu[grepl("SQUAMOUS|ESCC", hist), H := "ESCC"]
clu[grepl("ADENO|EAC",     hist), H := "EAC"]
clu[, CL := ifelse(grepl("1", as.character(cluster)), "C1", "C2")]   # normalize C1/C2
cat("histology in cluster table:\n"); print(table(clu$H, useNA="ifany"))
cat("cluster sizes:\n"); print(table(clu$cluster))
cat("cluster x histology:\n"); print(table(clu$cluster, clu$H))

## merge RS + cluster + histology
m <- merge(rs[, .(sampleID, RiskScore, RiskGroup, OS, OS.time)],
           clu[, .(sampleID, CL, H)], by="sampleID")
cat("\nmerged samples:", nrow(m), "\n")

## ---- A. RS by cluster and by histology ----
cat("\n=== A. Risk score distribution ===\n")
cat(sprintf("RS by cluster: C1 median %.3f (n=%d) vs C2 median %.3f (n=%d), Wilcoxon p=%.3g\n",
    median(m[CL=="C1"]$RiskScore), nrow(m[CL=="C1"]),
    median(m[CL=="C2"]$RiskScore), nrow(m[CL=="C2"]),
    wilcox.test(RiskScore~CL, m)$p.value))
mh <- m[!is.na(H)]
cat(sprintf("RS by histology: ESCC median %.3f (n=%d) vs EAC median %.3f (n=%d), Wilcoxon p=%.3g\n",
    median(mh[H=="ESCC"]$RiskScore), nrow(mh[H=="ESCC"]),
    median(mh[H=="EAC"]$RiskScore),  nrow(mh[H=="EAC"]),
    wilcox.test(RiskScore~H, mh)$p.value))

## ---- C. AUC of RS for discriminating histology (ESCC vs EAC) ----
auc <- function(score, pos){ r<-rank(score); (sum(r[pos])-sum(pos)*(sum(pos)+1)/2)/(sum(pos)*sum(!pos)) }
a_hist <- auc(mh$RiskScore, mh$H=="ESCC")
a_hist <- max(a_hist, 1-a_hist)
cat(sprintf("\n=== C. RS discriminating ESCC vs EAC: AUC = %.3f ===\n", a_hist))
cat("(compare: RS discriminating survival gave training C-index 0.709)\n")

## ---- B. per-gene C1vsC2 and ESCC vs EAC + contribution sign ----
cat("\n=== B. Per-gene alignment (log2 expr) ===\n")
E <- expr[gene %in% GENES]
E <- as.data.frame(E); rownames(E) <- E$gene; E$gene <- NULL
samp <- intersect(colnames(E), m$sampleID)
E <- E[GENES, samp, drop=FALSE]
mm <- m[match(samp, m$sampleID)]
res <- data.table(gene=GENES, beta=BETA[GENES])
for (g in GENES){
  x <- as.numeric(E[g, ]); cl <- mm$CL; H <- mm$H
  d_clu <- median(x[cl=="C1"]) - median(x[cl=="C2"])      # C1 - C2
  p_clu <- wilcox.test(x~cl)$p.value
  ok <- !is.na(H)
  d_his <- median(x[ok & H=="ESCC"]) - median(x[ok & H=="EAC"])   # ESCC - EAC
  p_his <- wilcox.test(x[ok]~H[ok])$p.value
  res[gene==g, `:=`(C1_minus_C2=round(d_clu,3), p_clu=signif(p_clu,2),
                    median_ESCC=round(median(x[ok & H=="ESCC"]),2), median_EAC=round(median(x[ok & H=="EAC"]),2),  # Table 3
                    ESCC_minus_EAC=round(d_his,3), p_his=signif(p_his,2),
                    higher_in=ifelse(d_his>0,"ESCC","EAC"),
                    # sign of this gene's push on RS toward ESCC:
                    RS_push=ifelse(BETA[g]*d_his>0,"raises RS in ESCC","raises RS in EAC"))]
}
print(res)

## summarize: of the 7 genes, how many raise RS in ESCC vs EAC, weighted
push_escc <- sum(res$beta * res$ESCC_minus_EAC)   # net directional push
cat(sprintf("\nNet directional term Sigma(beta * (ESCC-EAC median diff)) = %.4f\n", push_escc))
cat("(positive => the signature tends to give ESCC higher risk scores overall)\n")

## ---- D. CODA vs TCGA 18-FMG prevalence correlation ----
cat("\n=== D. CODA vs TCGA FMG prevalence ===\n")
coda <- c(TP53=32.7,ATM=9.2,FANCD2=9.2,LRP1B=9.2,NOTCH1=9.2,RUNX1=9.2,NOTCH3=8.2,
          PKHD1=7.1,BRCA2=6.1,NTRK1=6.1,PIK3CA=6.1,ASXL1=5.1,ERCC6=5.1,FAT3=5.1,
          FGFR2=5.1,KMT2D=5.1,ROS1=5.1,SLX4=5.1)
tcga <- c(TP53=85.9,ATM=4.9,FANCD2=1.1,LRP1B=15.2,NOTCH1=8.2,RUNX1=0.5,NOTCH3=6.0,
          PKHD1=2.2,BRCA2=3.8,NTRK1=2.2,PIK3CA=9.2,ASXL1=4.9,ERCC6=1.6,FAT3=10.9,
          FGFR2=0.0,KMT2D=11.4,ROS1=3.3,SLX4=1.1)
cat(sprintf("Spearman rho (CODA vs TCGA prevalence, 18 FMG) = %.3f, p=%.3g\n",
    cor(coda,tcga,method="spearman"), cor.test(coda,tcga,method="spearman")$p.value))
cat(sprintf("Pearson r = %.3f\n", cor(coda,tcga)))

###############################################################################
## Stage 4c (2026-10-19, released 2026-11-16): mutation status vs expression of the same FMG
## (Results 3.1, last paragraph; Table S1)
##
## For each FMG with >= 3 mutated tumors, expression (log2(norm_count + 1)) is compared between tumors with and without
## a functional mutation in that gene (Wilcoxon rank-sum, R default). Functional variants as in stage 04.
##   _180 columns: the 180 tumors with both MC3 calls and RNA-seq data — the values reported in the manuscript.
##   _181 columns: all 181 expression samples, with the one RNA-only EAC tumor (no MC3 calls) counted as wild-type,
##                 kept for comparison.
## Inputs  data/expr_ESCA_tumor.csv, data/ESCA_mc3.csv, tables/Tab1_2_FMG_list_primary.csv
## Output  tables/Tab04c_FMG_mutation_vs_expression.csv
###############################################################################
expr <- read.csv("data/expr_ESCA_tumor.csv", check.names = FALSE)
rownames(expr) <- expr$sample; expr$sample <- NULL
mc3 <- read.csv("data/ESCA_mc3.csv", stringsAsFactors = FALSE)
fmg <- read.csv("tables/Tab1_2_FMG_list_primary.csv", stringsAsFactors = FALSE)$gene
if (is.null(fmg)) stop("FMG list not found")
FUNC <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
          "In_Frame_Del", "In_Frame_Ins", "Translation_Start_Site", "Nonstop_Mutation")
mut <- unique(mc3[mc3$gene %in% fmg & mc3$effect %in% FUNC, c("sample", "gene")])
run <- function(samples, label) {
  out <- do.call(rbind, lapply(fmg, function(g) {
    m <- samples %in% mut$sample[mut$gene == g]
    if (sum(m) < 3) return(NULL)
    x <- as.numeric(expr[g, samples])
    data.frame(gene = g, n_mut = sum(m), n_wt = sum(!m), med_mut = round(median(x[m]), 2),
               med_wt = round(median(x[!m]), 2), p = signif(wilcox.test(x[m], x[!m])$p.value, 3))
  }))
  out <- out[order(out$p), ]
  cat("\n==", label, ":", length(samples), "samples,", nrow(out), "FMG with >= 3 mutated\n")
  print(out, row.names = FALSE)
  invisible(out)
}
all181 <- colnames(expr)
with_mc3 <- intersect(all181, unique(mc3$sample))
a <- run(all181, "all 181 expression samples (RNA-only tumor as wild-type)")
b <- run(with_mc3, "180 samples with MC3 data (manuscript)")
write.csv(merge(a, b, by = "gene", suffixes = c("_181", "_180")),
          "tables/Tab04c_FMG_mutation_vs_expression.csv", row.names = FALSE)
cat("[OK] tables/Tab04c_FMG_mutation_vs_expression.csv\n")

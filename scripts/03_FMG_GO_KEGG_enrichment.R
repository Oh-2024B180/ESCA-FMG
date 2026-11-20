###############################################################################
## Stage 3: FMG GO/KEGG Functional Enrichment
##
## - 3A: clusterProfiler R package (>=4.0) [PMID: 34557778]
## - 3B: P-value correction = Benjamini-Hochberg, adj.p < 0.05
## - 3C: Background = org.Hs.eg.db default universe
## SYMBOL -> ENTREZID conversion before enrichment
## - 3D: Visualization = dotplot top 10 per ontology (BP/MF/CC + KEGG)
##
## Input:
## tables/Tab1_2_FMG_list_primary.csv (18 Primary FMG)
##
## Outputs (tables/):
## Tab3_FMG_GO_BP.csv (full enrichGO BP result + significant flag)
## Tab3_FMG_GO_MF.csv (full enrichGO MF result + significant flag)
## Tab3_FMG_GO_CC.csv (full enrichGO CC result + significant flag)
## Tab3_FMG_KEGG.csv (full enrichKEGG result + significant flag)
##
## Outputs (figures/, PDF + PNG @ 1200 dpi):
## Fig1B_FMG_GO_BP_dotplot.{pdf,png} (body Fig 1B, panel: GO BP)
## Fig1B_FMG_GO_CC_dotplot.{pdf,png} (body Fig 1B, panel: GO CC)
## Fig1B_FMG_GO_MF_dotplot.{pdf,png} (body Fig 1B, panel: GO MF supplementary)
## Fig1C_FMG_KEGG_dotplot.{pdf,png} (body Fig 1C)
##
## Reference PMID:
## - 34557778 : Wu T et al. 2021 Innovation -- clusterProfiler 4.0
##
## Random seed: 1262118388 (Liu et al. 2024)
##
## Notes:
## - 18 FMG are a small enrichment input; if no term passes BH adj.p < 0.05 (e.g. KEGG),
## the dot plot is skipped and this is reported.
## - enrichKEGG calls the KEGG REST API on the first run (internet access required).
## - The full (unfiltered) results are written to CSV; column `significant_BH_05` marks significance.
###############################################################################

# 0. Working directory + output directory
# run from the repository root (see README.md)

INPUT_FMG <- "tables/Tab1_2_FMG_list_primary.csv"
TBLDIR <- "tables"
FIGDIR <- "figures"
LOGDIR <- "logs"

for (d in c(TBLDIR, FIGDIR, LOGDIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

ADJ_P_CUTOFF <- 0.05 # Stage 3B
TOP_N <- 10 # Stage 3D

# 0.5. Reproducibility
set.seed(1262118388)

# 0.6. Logging start
log_file <- file.path(LOGDIR, "03_diagnostics.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat(sprintf("=== Stage 3 Pipeline start (%s) ===\n", Sys.time()))
cat(sprintf("Input : %s\n", INPUT_FMG))
cat(sprintf("BH adj.p cut : %.3f\n", ADJ_P_CUTOFF))
cat(sprintf("Top N (plot) : %d\n\n", TOP_N))

# 1. load
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(ragg)
})

cat(sprintf("clusterProfiler version : %s\n", packageVersion("clusterProfiler")))
cat(sprintf("org.Hs.eg.db version : %s\n", packageVersion("org.Hs.eg.db")))
cat(sprintf("enrichplot version : %s\n", packageVersion("enrichplot")))
cat(sprintf("ggplot2 version : %s\n\n", packageVersion("ggplot2")))

# 2. Load the FMG list + SYMBOL -> ENTREZID
cat(">>> STEP 2: Load FMG list and convert SYMBOL -> ENTREZID\n")
fmg <- read.csv(INPUT_FMG, stringsAsFactors = FALSE)
stopifnot("gene" %in% colnames(fmg))
genes_symbol <- unique(fmg$gene)
cat(sprintf(" Input FMG count: %d\n", length(genes_symbol)))

gene_map <- bitr(genes_symbol,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

dropped <- setdiff(genes_symbol, gene_map$SYMBOL)
if (length(dropped) > 0) {
  cat(sprintf(" WARNING: %d gene(s) failed SYMBOL->ENTREZID mapping: %s\n",
              length(dropped), paste(dropped, collapse = ", ")))
} else {
  cat(" All FMG mapped to ENTREZID successfully.\n")
}
genes_entrez <- gene_map$ENTREZID
cat(sprintf(" ENTREZID count : %d\n\n", length(genes_entrez)))

# 3. Enrichment (GO BP/MF/CC + KEGG)
# result pvalueCutoff=1, qvalueCutoff=1 use.
# Significance column(significant_BH_05) display.
run_enrichGO <- function(ont) {
  enrichGO(
    gene = genes_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
}

cat(">>> STEP 3a: enrichGO BP\n")
ego_bp <- run_enrichGO("BP")
cat(sprintf(" total terms returned: %d\n", nrow(as.data.frame(ego_bp))))

cat(">>> STEP 3b: enrichGO MF\n")
ego_mf <- run_enrichGO("MF")
cat(sprintf(" total terms returned: %d\n", nrow(as.data.frame(ego_mf))))

cat(">>> STEP 3c: enrichGO CC\n")
ego_cc <- run_enrichGO("CC")
cat(sprintf(" total terms returned: %d\n", nrow(as.data.frame(ego_cc))))

cat(">>> STEP 3d: enrichKEGG (organism = hsa)\n")
ekegg <- tryCatch({
  res <- enrichKEGG(
    gene = genes_entrez,
    organism = "hsa",
    keyType = "kegg",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  if (!is.null(res)) res <- setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  res
}, error = function(e) {
  cat(sprintf(" ERROR: enrichKEGG failed: %s\n", conditionMessage(e)))
  NULL
})
if (!is.null(ekegg)) {
  cat(sprintf(" total terms returned: %d\n\n", nrow(as.data.frame(ekegg))))
} else {
  cat(" KEGG result is NULL.\n\n")
}

# 4. Result tables
write_enrich_table <- function(enrich_obj, out_file) {
  if (is.null(enrich_obj)) {
    cat(sprintf(" [SKIP] %s : enrichment object is NULL\n", basename(out_file)))
    return(invisible(NULL))
  }
  df <- as.data.frame(enrich_obj)
  if (nrow(df) == 0) {
    cat(sprintf(" [WARN] %s : 0 enrichment rows\n", basename(out_file)))
    write.csv(df, out_file, row.names = FALSE)
    return(invisible(NULL))
  }
  df$significant_BH_05 <- df$p.adjust < ADJ_P_CUTOFF
  write.csv(df, out_file, row.names = FALSE)
  n_sig <- sum(df$significant_BH_05, na.rm = TRUE)
  cat(sprintf(" [OK] %s : %d rows (%d significant @ BH adj.p<%.2f)\n",
              basename(out_file), nrow(df), n_sig, ADJ_P_CUTOFF))
  invisible(df)
}

cat(">>> STEP 4: Write enrichment tables\n")
df_bp <- write_enrich_table(ego_bp, file.path(TBLDIR, "Tab3_FMG_GO_BP.csv"))
df_mf <- write_enrich_table(ego_mf, file.path(TBLDIR, "Tab3_FMG_GO_MF.csv"))
df_cc <- write_enrich_table(ego_cc, file.path(TBLDIR, "Tab3_FMG_GO_CC.csv"))
df_kegg <- write_enrich_table(ekegg, file.path(TBLDIR, "Tab3_FMG_KEGG.csv"))
cat("\n")

# 5. Dot plots (top N significant terms)
make_dotplot <- function(enrich_obj, title_str, fig_stem) {
  if (is.null(enrich_obj)) {
    cat(sprintf(" [SKIP] %s : object NULL\n", fig_stem))
    return(invisible(NULL))
  }
  full_df <- as.data.frame(enrich_obj)
  sig_df <- full_df[!is.na(full_df$p.adjust) & full_df$p.adjust < ADJ_P_CUTOFF, , drop = FALSE]
  if (nrow(sig_df) == 0) {
    cat(sprintf(" [SKIP] %s : 0 significant terms (BH adj.p<%.2f) -> dotplot not generated\n",
                fig_stem, ADJ_P_CUTOFF))
    return(invisible(NULL))
  }

  obj_sig <- enrich_obj
  obj_sig@result <- sig_df

  show_n <- min(TOP_N, nrow(sig_df))
  p <- dotplot(obj_sig, showCategory = show_n) +
       theme_bw(base_size = 10)

  pdf_file <- file.path(FIGDIR, paste0(fig_stem, ".pdf"))
  png_file <- file.path(FIGDIR, paste0(fig_stem, ".png"))

  ggsave(pdf_file, plot = p, width = 7, height = 6, device = cairo_pdf)
  ggsave(png_file, plot = p, width = 7, height = 6, dpi = 1200, device = ragg::agg_png)

  cat(sprintf(" [OK] %s.{pdf,png} : top %d of %d sig terms\n",
              fig_stem, show_n, nrow(sig_df)))
  invisible(p)
}

cat(">>> STEP 5: Generate dotplots\n")
make_dotplot(ego_bp, "FMG GO Biological Process", "Fig1B_FMG_GO_BP_dotplot")
make_dotplot(ego_cc, "FMG GO Cellular Component", "Fig1B_FMG_GO_CC_dotplot")
make_dotplot(ego_mf, "FMG GO Molecular Function", "Fig1B_FMG_GO_MF_dotplot")
make_dotplot(ekegg, "FMG KEGG Pathway", "Fig1C_FMG_KEGG_dotplot")
cat("\n")

# 6. Summary Report
top_terms <- function(df, label, n = 5) {
  if (is.null(df) || nrow(df) == 0) {
    cat(sprintf(" %s: (no rows)\n", label))
    return(invisible(NULL))
  }
  sig <- df[!is.na(df$p.adjust) & df$p.adjust < ADJ_P_CUTOFF, , drop = FALSE]
  cat(sprintf("\n --- %s ---\n", label))
  cat(sprintf(" total terms : %d\n", nrow(df)))
  cat(sprintf(" significant (BH<%.2f): %d\n", ADJ_P_CUTOFF, nrow(sig)))
  if (nrow(sig) > 0) {
    show_n <- min(n, nrow(sig))
    cat(sprintf(" Top %d (by p.adjust):\n", show_n))
    show_df <- sig[order(sig$p.adjust), c("ID","Description","Count","GeneRatio","p.adjust"), drop = FALSE]
    print(head(show_df, show_n), row.names = FALSE)
  }
}

cat("\n")
cat("========================================================================\n")
cat(" Stage 3 FMG GO/KEGG Summary Report \n")
cat("========================================================================\n")
cat(sprintf(" Input FMG : %d\n", length(genes_symbol)))
cat(sprintf(" ENTREZID mapped : %d\n", length(genes_entrez)))
cat(sprintf(" SYMBOL -> ENTREZID drop : %d\n", length(dropped)))
cat(sprintf(" BH adj.p cutoff : %.3f\n", ADJ_P_CUTOFF))
top_terms(df_bp, "GO Biological Process")
top_terms(df_mf, "GO Molecular Function")
top_terms(df_cc, "GO Cellular Component")
top_terms(df_kegg, "KEGG Pathway")
cat("\n========================================================================\n")

# 7. sessionInfo
sess_file <- file.path(LOGDIR, "sessionInfo/sessionInfo_stage03.txt")
sink(sess_file)
cat(sprintf("Session info captured: %s\n\n", Sys.time()))
print(sessionInfo())
sink()
cat(sprintf("\n>>> sessionInfo saved: %s\n", sess_file))

cat(sprintf("\n=== Stage 3 Pipeline done (%s) ===\n", Sys.time()))

# end: log file close
sink(type = "message")
sink()
close(log_con)

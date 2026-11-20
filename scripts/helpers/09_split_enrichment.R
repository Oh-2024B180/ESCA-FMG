###############################################################################
##  Stage 9B — Direction-resolved enrichment (16 enrichments)
##
##  Splits primary + sensitivity sig DEGs by logFC sign and runs enrichGO
##  BP/MF/CC + enrichKEGG separately for each of 4 subsets:
##    primary_C1up, primary_C2up, sensitivity_C1up, sensitivity_C2up
##
##  Inputs:
##    tables/Tab9_DEG_significant_C1_vs_C2.csv      (Primary sig)
##    tables/Tab9_DEG_sensitivity_significant.csv   (Sensitivity sig)
##
##  Outputs:
##    Tab9_DEG_<ontology>_<track>_<dir>.csv  (16 tables)
##    Fig5E_DEG_<ontology>_<track>_<dir>_dotplot.{pdf,png}  (16 figures × 2)
##  Log: 09B_split_enrichment_diagnostics.txt
##
##  Random seed: 1262118388
###############################################################################

# run from the repository root (see README.md)
suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(ggplot2)
  library(ragg)
})
set.seed(1262118388)

INPUT_PRIMARY     <- "tables/Tab9_DEG_significant_C1_vs_C2.csv"
INPUT_SENSITIVITY <- "tables/Tab9_DEG_sensitivity_significant.csv"
TBLDIR <- "tables"; FIGDIR <- "figures"; LOGDIR <- "logs"
DIAG_LOG <- file.path(LOGDIR, "09B_split_enrichment_diagnostics.txt")
ADJ_P_CUT  <- 0.05
TOP_N_PLOT <- 10
KEGG_RETRY_MAX <- 1

t_start <- Sys.time()
diag_lines <- c(sprintf("=== Stage 9B split enrichment === %s",
                         format(t_start, "%Y-%m-%d %H:%M:%S")), "")

cat("[1] Loading sig DEG tables\n")
read_sig <- function(path, label) {
  dt <- fread(path)
  dt[, source_track := label]
  dt[, c("gene","log2FC_C1minusC2","p_adj_BH","source_track"), with = FALSE]
}
primary     <- read_sig(INPUT_PRIMARY,     "primary")
sensitivity <- read_sig(INPUT_SENSITIVITY, "sensitivity")
cat(sprintf("    Primary sig: %d  | Sensitivity sig: %d\n",
            nrow(primary), nrow(sensitivity)))

cat("\n[2] Splitting by direction\n")
subsets <- list(
  primary_C1up     = primary[log2FC_C1minusC2 >  0, gene],
  primary_C2up     = primary[log2FC_C1minusC2 <  0, gene],
  sensitivity_C1up = sensitivity[log2FC_C1minusC2 >  0, gene],
  sensitivity_C2up = sensitivity[log2FC_C1minusC2 <  0, gene]
)
for (nm in names(subsets)) cat(sprintf("    %-22s : %d\n", nm, length(subsets[[nm]])))
diag_lines <- c(diag_lines,
                sprintf("Subset gene counts: %s",
                        paste(sprintf("%s=%d", names(subsets), lengths(subsets)),
                              collapse = " | ")), "")

cat("\n[3] SYMBOL -> ENTREZID per subset\n")
entrez_lists <- list()
for (nm in names(subsets)) {
  res <- suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = subsets[[nm]],
                          columns = c("SYMBOL","ENTREZID"), keytype = "SYMBOL"))
  res <- res[!is.na(res$ENTREZID), ]
  res <- res[!duplicated(res$SYMBOL), ]
  entrez_lists[[nm]] <- res$ENTREZID
  cat(sprintf("    %-22s : %d -> %d ENTREZID (drop %d)\n",
              nm, length(subsets[[nm]]), nrow(res), length(subsets[[nm]]) - nrow(res)))
}

cat("\n[4] Running 16 enrichments (4 subsets × 4 ontologies)\n")
ontologies <- c("BP","MF","CC","KEGG")
result_summary  <- data.table()
top_results     <- list()
kegg_retries_used <- 0

run_enrichment <- function(entrez_ids, ontology) {
  if (length(entrez_ids) < 5) return(NULL)
  if (ontology == "KEGG") {
    res <- NULL
    for (try_i in seq_len(KEGG_RETRY_MAX + 1)) {
      res <- tryCatch({
        r <- enrichKEGG(gene = entrez_ids, organism = "hsa", keyType = "kegg",
                        pAdjustMethod = "BH", pvalueCutoff = ADJ_P_CUT)
        if (!is.null(r)) r <- setReadable(r, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
        r
      }, error = function(e) {
        cat(sprintf("      [warn] KEGG attempt %d: %s\n", try_i, conditionMessage(e)))
        Sys.sleep(5); NULL
      })
      if (!is.null(res)) break
      if (try_i > 1) kegg_retries_used <<- kegg_retries_used + 1
    }
    return(res)
  }
  enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
           ont = ontology, pAdjustMethod = "BH",
           pvalueCutoff = ADJ_P_CUT, qvalueCutoff = 1, readable = TRUE)
}

write_enrich_table <- function(obj, out_file) {
  if (is.null(obj)) { fwrite(data.table(message = "NULL or insufficient input"), out_file); return(0L) }
  df <- as.data.frame(obj)
  if (nrow(df) == 0) { fwrite(data.table(message = "0 sig"), out_file); return(0L) }
  fwrite(df, out_file); nrow(df)
}

make_dotplot <- function(obj, title_str, fig_stem, n_sig) {
  if (is.null(obj) || n_sig == 0) {
    cat(sprintf("        [SKIP] %s : 0 sig\n", basename(fig_stem))); return(invisible(NULL))
  }
  show_n <- min(TOP_N_PLOT, n_sig)
  pp <- tryCatch({
    dotplot(obj, showCategory = show_n) +
      theme_bw(base_size = 10)
  }, error = function(e) NULL)
  if (is.null(pp)) return(invisible(NULL))
  ggsave(paste0(fig_stem, ".pdf"), pp, width = 7, height = 6, device = cairo_pdf)
  ggsave(paste0(fig_stem, ".png"), pp, width = 7, height = 6, dpi = 600, device = ragg::agg_png)
  cat(sprintf("        [OK] %s.{pdf,png} (top %d of %d)\n",
              basename(fig_stem), show_n, n_sig))
}

i_counter <- 0
total_n <- length(subsets) * length(ontologies)
for (subset_name in names(subsets)) {
  for (ontology in ontologies) {
    i_counter <- i_counter + 1
    cat(sprintf("\n  [%2d/%d] subset=%s ont=%s (n_entrez=%d)\n",
                i_counter, total_n, subset_name, ontology, length(entrez_lists[[subset_name]])))
    obj <- run_enrichment(entrez_lists[[subset_name]], ontology)
    tab_file <- file.path(TBLDIR,
      sprintf("Tab9_DEG_%s_%s.csv",
              if (ontology == "KEGG") "KEGG" else paste0("GO_", ontology),
              subset_name))
    n_sig <- write_enrich_table(obj, tab_file)
    cat(sprintf("        [OK] %s : %d sig (BH<%.2f)\n", basename(tab_file), n_sig, ADJ_P_CUT))

    fig_stem <- file.path(FIGDIR,
      sprintf("Fig5E_DEG_%s_%s_dotplot",
              if (ontology == "KEGG") "KEGG" else paste0("GO_", ontology),
              subset_name))
    title_str <- sprintf("DEG %s — %s",
                          if (ontology == "KEGG") "KEGG"
                          else sprintf("GO %s",
                                       switch(ontology, BP = "Biological Process",
                                              MF = "Molecular Function",
                                              CC = "Cellular Component")),
                          subset_name)
    make_dotplot(obj, title_str, fig_stem, n_sig)

    result_summary <- rbind(result_summary,
                             data.table(subset = subset_name, ontology = ontology, n_sig = n_sig))
    if (!is.null(obj) && n_sig > 0) {
      df <- as.data.frame(obj)
      top_results[[paste(subset_name, ontology, sep="_")]] <-
        head(df[order(df$p.adjust),
                c("ID","Description","Count","GeneRatio","p.adjust"), drop = FALSE], 10)
    }
  }
}

t_end <- Sys.time()
elapsed_min <- as.numeric(t_end - t_start, units = "mins")
cat("\n========================================================================\n")
cat("        Stage 9B Split Enrichment Summary\n")
cat("========================================================================\n")
cat(sprintf("  Elapsed: %.1f min, KEGG retries: %d\n\n", elapsed_min, kegg_retries_used))
ssum <- dcast(result_summary, subset ~ ontology, value.var = "n_sig")
print(ssum)
cat("\n  Sensitivity_C1up GO BP top 10:\n")
if (!is.null(top_results[["sensitivity_C1up_BP"]])) print(top_results[["sensitivity_C1up_BP"]])
cat("\n  Primary_C2up GO BP top 5:\n")
if (!is.null(top_results[["primary_C2up_BP"]])) print(head(top_results[["primary_C2up_BP"]], 5))
cat("========================================================================\n")

diag_lines <- c(diag_lines,
  "Significant term counts (BH<0.05):",
  capture.output(print(ssum)), "",
  sprintf("Elapsed: %.1f min, KEGG retries: %d", elapsed_min, kegg_retries_used),
  "", sprintf("=== END %s ===", format(t_end, "%Y-%m-%d %H:%M:%S")))
writeLines(diag_lines, DIAG_LOG)
cat(sprintf(">>> diagnostics: %s\n", DIAG_LOG))
cat(sprintf("\n=== Stage 9B done in %.1f min ===\n", elapsed_min))

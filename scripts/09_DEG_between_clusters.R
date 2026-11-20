###############################################################################
## Stage 9: DEG Between Clusters (TWO-TRACK + biotype filter)
##
## Inputs:
## data/expr_ESCA_tumor.csv
## tables/Tab7_cluster_assignment.csv
##
## Outputs (tables/, flat):
## Tab9_DEG_full_C1_vs_C2.csv
## Tab9_DEG_significant_C1_vs_C2.csv
## Tab9_DEG_sensitivity_full.csv
## Tab9_DEG_sensitivity_significant.csv
## Tab9_DEG_overlap_summary.csv
## Tab9_biotype_filter_summary.csv (biotype distribution)
## Tab9_DEG_GO_BP / MF / CC.csv Tab9_DEG_KEGG.csv
##
## Outputs (figures/):
## Fig5A_DEG_volcano.{pdf,png}
## Fig5B_DEG_heatmap_top25up_top25down.{pdf,png}
## Fig5C_DEG_GO_{BP,MF,CC}_dotplot.{pdf,png}
## Fig5D_DEG_KEGG_dotplot.{pdf,png}
###############################################################################

# run from the repository root (see README.md)

INPUT_EXPR <- "data/expr_ESCA_tumor.csv"
INPUT_CLUSTER <- "tables/Tab7_cluster_assignment.csv"

TBLDIR <- "tables"
FIGDIR <- "figures"
LOGDIR <- "logs"
CACHE <- "cache"
for (d in c(TBLDIR, FIGDIR, LOGDIR, CACHE)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

TAB_FULL_P <- file.path(TBLDIR, "Tab9_DEG_full_C1_vs_C2.csv")
TAB_SIG_P <- file.path(TBLDIR, "Tab9_DEG_significant_C1_vs_C2.csv")
TAB_FULL_S <- file.path(TBLDIR, "Tab9_DEG_sensitivity_full.csv")
TAB_SIG_S <- file.path(TBLDIR, "Tab9_DEG_sensitivity_significant.csv")
TAB_OVERLAP <- file.path(TBLDIR, "Tab9_DEG_overlap_summary.csv")
TAB_BIOTYPE <- file.path(TBLDIR, "Tab9_biotype_filter_summary.csv")
TAB_GO_BP <- file.path(TBLDIR, "Tab9_DEG_GO_BP.csv")
TAB_GO_MF <- file.path(TBLDIR, "Tab9_DEG_GO_MF.csv")
TAB_GO_CC <- file.path(TBLDIR, "Tab9_DEG_GO_CC.csv")
TAB_KEGG <- file.path(TBLDIR, "Tab9_DEG_KEGG.csv")

FIG_VOLCANO <- file.path(FIGDIR, "Fig5A_DEG_volcano")
FIG_HEATMAP <- file.path(FIGDIR, "Fig5B_DEG_heatmap_top25up_top25down")
FIG_GO_BP <- file.path(FIGDIR, "Fig5C_DEG_GO_BP_dotplot")
FIG_GO_MF <- file.path(FIGDIR, "Fig5C_DEG_GO_MF_dotplot")
FIG_GO_CC <- file.path(FIGDIR, "Fig5C_DEG_GO_CC_dotplot")
FIG_KEGG <- file.path(FIGDIR, "Fig5D_DEG_KEGG_dotplot")

LOGFC_CUT <- 1.0
ADJ_P_CUT <- 0.05
ZERO_FRAC_CUT <- 0.70
TOP_N_GO <- 10

CLUSTER_COL <- c(C1 = "#1B9E77", C2 = "#D95F02") # Stage 7 unify
HIST_COL <- c(ESCC = "#7570B3", EAC = "#E7298A", Other = "#66A61E")

# biomaRt hosts: the October 2024 Ensembl archive first. It is the release given in the Methods (2.7) and the
# only host that answered in the original run (16,921 protein-coding genes, 2,304 DEG); the current release could
# classify gene biotypes differently. The other hosts stay as fallbacks.
BIOMART_HOSTS <- c(
  "https://oct2024.archive.ensembl.org",
  "https://www.ensembl.org",
  "https://useast.ensembl.org",
  "https://uswest.ensembl.org",
  "https://asia.ensembl.org",
  "https://nov2024.archive.ensembl.org"
)
BIOTYPE_CACHE <- file.path(CACHE, "cache_biomart_HGNC_biotype_stage09.rds")

set.seed(1262118388)

# Logging
log_file <- file.path(LOGDIR, "09_diagnostics.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
cat(sprintf("=== Stage 9 Pipeline start (%s) ===\n\n", Sys.time()))

# 1. Packages
cat(">>> STEP 1: Loading packages\n")
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(biomaRt)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(ragg)
  library(grid)
})
cat(sprintf(" limma : %s\n", packageVersion("limma")))
cat(sprintf(" biomaRt : %s\n", packageVersion("biomaRt")))
cat(sprintf(" clusterProfiler : %s\n\n", packageVersion("clusterProfiler")))

# 2. Inputs
cat(">>> STEP 2: Load expression + cluster\n")
expr_dt <- fread(INPUT_EXPR)
cluster_dt <- fread(INPUT_CLUSTER)
cat(sprintf(" expr : %d genes x %d cols\n", nrow(expr_dt), ncol(expr_dt)))
cat(sprintf(" cluster: %d samples (C1=%d, C2=%d)\n",
            nrow(cluster_dt),
            sum(cluster_dt$cluster == "C1"), sum(cluster_dt$cluster == "C2")))

key_g <- colnames(expr_dt)[1]
setnames(expr_dt, key_g, "gene_id")
samples_keep <- intersect(cluster_dt$sampleID, colnames(expr_dt))
stopifnot(length(samples_keep) == nrow(cluster_dt))

expr_mat <- as.matrix(expr_dt[, ..samples_keep])
rownames(expr_mat) <- expr_dt$gene_id

pheno <- cluster_dt[match(samples_keep, sampleID)]
stopifnot(all(pheno$sampleID == colnames(expr_mat)))
pheno[, cluster := factor(cluster, levels = c("C1","C2"))]
pheno[, histology := factor(histology, levels = c("ESCC","EAC","Other"))]
pheno[, histology := droplevels(histology)]
cat(sprintf(" Histology (post-droplevels): %s\n",
            paste(sprintf("%s=%d", levels(pheno$histology),
                          as.integer(table(pheno$histology))), collapse = ", ")))

# 3. STEP 3a: Zero-expression filter
cat("\n>>> STEP 3a: Zero-expression filter (>70% zero drop)\n")
n_genes_in <- nrow(expr_mat)
zero_frac <- rowSums(expr_mat == 0) / ncol(expr_mat)
keep_zero <- zero_frac <= ZERO_FRAC_CUT
expr_z <- expr_mat[keep_zero, ]
cat(sprintf(" Genes before zero filter : %d\n", n_genes_in))
cat(sprintf(" Genes dropped (zero>70%%) : %d (%.1f%%)\n",
            n_genes_in - nrow(expr_z),
            100 * (n_genes_in - nrow(expr_z)) / n_genes_in))
cat(sprintf(" Genes kept : %d\n", nrow(expr_z)))

# STEP 3b: biomaRt biotype query
cat("\n>>> STEP 3b: biomaRt biotype query\n")

if (file.exists(BIOTYPE_CACHE)) {
  cat(sprintf(" [cache hit] %s\n", BIOTYPE_CACHE))
  biotype_dt <- readRDS(BIOTYPE_CACHE)
} else {
  cat(" [cache miss] querying Ensembl biomaRt for hgnc_symbol gene_biotype\n")
  cat(sprintf(" query size: %d unique SYMBOLs\n", length(unique(rownames(expr_z)))))

  biotype_raw <- NULL
  for (h in BIOMART_HOSTS) {
    cat(sprintf(" trying host: %s ...\n", h))
    biotype_raw <- tryCatch({
      mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                      dataset = "hsapiens_gene_ensembl", host = h)
      getBM(
        attributes = c("hgnc_symbol","gene_biotype"),
        filters = "hgnc_symbol",
        values = unique(rownames(expr_z)),
        mart = mart
      )
    }, error = function(e) {
      cat(sprintf(" [fail] %s %s\n", h, conditionMessage(e)))
      NULL
    })
    if (!is.null(biotype_raw) && nrow(biotype_raw) > 0) {
      cat(sprintf(" [SUCCESS] %s returned %d rows\n", h, nrow(biotype_raw)))
      break
    }
    Sys.sleep(2)
  }
  if (is.null(biotype_raw) || nrow(biotype_raw) == 0) {
    stop("biomaRt failed on all mirrors cannot apply biotype filter.")
  }

  biotype_dt <- as.data.table(biotype_raw)
  setnames(biotype_dt, c("hgnc_symbol","gene_biotype"), c("gene","biotype"))
  saveRDS(biotype_dt, BIOTYPE_CACHE)
  cat(sprintf(" [cached] %s\n", BIOTYPE_CACHE))
}
cat(sprintf(" biomaRt mappings: %d (gene, biotype) rows\n", nrow(biotype_dt)))

# Some SYMBOLs may have multiple biotypes (different Ensembl IDs).
# Rule: keep gene if ANY entry is protein_coding.
gene_pc_dt <- biotype_dt[, .(any_pc = any(biotype == "protein_coding")), by = gene]
pc_symbols <- gene_pc_dt[any_pc == TRUE, gene]
cat(sprintf(" unique SYMBOLs in biomaRt : %d\n", nrow(gene_pc_dt)))
cat(sprintf(" protein_coding SYMBOLs : %d\n", length(pc_symbols)))

# Biotype distribution table for documentation
biotype_dist <- biotype_dt[, .N, by = biotype][order(-N)]
cat(" Biotype distribution (top 15):\n")
print(head(biotype_dist, 15))
fwrite(biotype_dist, TAB_BIOTYPE)

# STEP 3c: Apply biotype filter
cat("\n>>> STEP 3c: Apply biotype == 'protein_coding' filter\n")
keep_pc <- rownames(expr_z) %in% pc_symbols
expr_filt <- expr_z[keep_pc, ]
n_after_zero <- nrow(expr_z)
n_after_pc <- nrow(expr_filt)
n_dropped_pc <- n_after_zero - n_after_pc
cat(sprintf(" Genes after zero filter : %d\n", n_after_zero))
cat(sprintf(" Genes dropped (non-PC) : %d (%.1f%%)\n",
            n_dropped_pc, 100 * n_dropped_pc / n_after_zero))
cat(sprintf(" Genes kept (PC only) : %d\n", n_after_pc))
cat(sprintf(" Total drop from input : %d %d (%.1f%% kept)\n\n",
            n_genes_in, n_after_pc, 100 * n_after_pc / n_genes_in))

# 4. PRIMARY DEG
cat(">>> STEP 4: PRIMARY DEG (~ 0 + cluster, contrast C1 - C2, eBayes trend=TRUE)\n")
design_p <- model.matrix(~ 0 + cluster, data = pheno)
colnames(design_p) <- gsub("^cluster", "", colnames(design_p))
contrast_p <- makeContrasts(C1_vs_C2 = C1 - C2, levels = design_p)
fit_p <- lmFit(expr_filt, design_p)
fit2_p <- contrasts.fit(fit_p, contrast_p)
fit2_p <- eBayes(fit2_p, trend = TRUE)
deg_p <- topTable(fit2_p, coef = "C1_vs_C2", number = Inf, sort.by = "none",
                   adjust.method = "BH")
deg_p_dt <- data.table(gene = rownames(deg_p), deg_p)
setnames(deg_p_dt,
         old = c("logFC","P.Value","adj.P.Val"),
         new = c("log2FC_C1minusC2","p_value","p_adj_BH"))
deg_p_dt[, sig := abs(log2FC_C1minusC2) >= LOGFC_CUT & p_adj_BH < ADJ_P_CUT]
deg_p_dt[, direction := fifelse(!sig, "ns",
                        fifelse(log2FC_C1minusC2 > 0, "C1>C2",
                                                                 "C2>C1"))]
setorder(deg_p_dt, p_adj_BH)
cat(sprintf(" Total tested : %d\n", nrow(deg_p_dt)))
n_primary_sig <- sum(deg_p_dt$sig)
cat(sprintf(" Significant DEG : %d (|log2FC|>=%.1f AND BH<%.2f)\n",
            n_primary_sig, LOGFC_CUT, ADJ_P_CUT))
cat(sprintf(" C1>C2 (up) : %d\n", sum(deg_p_dt$direction == "C1>C2")))
cat(sprintf(" C2>C1 (down) : %d\n", sum(deg_p_dt$direction == "C2>C1")))
fwrite(deg_p_dt, TAB_FULL_P)
fwrite(deg_p_dt[sig == TRUE], TAB_SIG_P)
cat(sprintf(" [OK] %s\n [OK] %s\n", TAB_FULL_P, TAB_SIG_P))

# 5. SENSITIVITY DEG
cat("\n>>> STEP 5: SENSITIVITY DEG (~ 0 + cluster + histology, same contrast)\n")
design_s <- model.matrix(~ 0 + cluster + histology, data = pheno)
colnames(design_s) <- gsub("^cluster", "", colnames(design_s))
colnames(design_s) <- gsub("^histology", "hist_", colnames(design_s))
contrast_s <- makeContrasts(C1_vs_C2 = C1 - C2, levels = design_s)
fit_s <- lmFit(expr_filt, design_s)
fit2_s <- contrasts.fit(fit_s, contrast_s)
fit2_s <- eBayes(fit2_s, trend = TRUE)
deg_s <- topTable(fit2_s, coef = "C1_vs_C2", number = Inf, sort.by = "none",
                   adjust.method = "BH")
deg_s_dt <- data.table(gene = rownames(deg_s), deg_s)
setnames(deg_s_dt,
         old = c("logFC","P.Value","adj.P.Val"),
         new = c("log2FC_C1minusC2","p_value","p_adj_BH"))
deg_s_dt[, sig := abs(log2FC_C1minusC2) >= LOGFC_CUT & p_adj_BH < ADJ_P_CUT]
deg_s_dt[, direction := fifelse(!sig, "ns",
                        fifelse(log2FC_C1minusC2 > 0, "C1>C2",
                                                                 "C2>C1"))]
setorder(deg_s_dt, p_adj_BH)
n_sens_sig <- sum(deg_s_dt$sig)
cat(sprintf(" Total tested : %d\n", nrow(deg_s_dt)))
cat(sprintf(" Significant DEG : %d\n", n_sens_sig))
cat(sprintf(" C1>C2 (up) : %d\n", sum(deg_s_dt$direction == "C1>C2")))
cat(sprintf(" C2>C1 (down) : %d\n", sum(deg_s_dt$direction == "C2>C1")))
fwrite(deg_s_dt, TAB_FULL_S)
fwrite(deg_s_dt[sig == TRUE], TAB_SIG_S)
cat(sprintf(" [OK] %s\n [OK] %s\n", TAB_FULL_S, TAB_SIG_S))

# 6. Overlap summary
cat("\n>>> STEP 6: Primary vs Sensitivity overlap\n")
p_min <- deg_p_dt[, .(gene, primary_logFC = log2FC_C1minusC2,
                       primary_p_adj = p_adj_BH, primary_sig = sig)]
s_min <- deg_s_dt[, .(gene, sensitivity_logFC = log2FC_C1minusC2,
                       sensitivity_p_adj = p_adj_BH, sensitivity_sig = sig)]
overlap_dt <- merge(p_min, s_min, by = "gene", all = TRUE)
overlap_dt[, overlap_status := fifelse(primary_sig & sensitivity_sig, "both",
                              fifelse(primary_sig & !sensitivity_sig, "primary_only",
                              fifelse(!primary_sig & sensitivity_sig, "sensitivity_only",
                                                                     "neither")))]
setorder(overlap_dt, primary_p_adj)
fwrite(overlap_dt, TAB_OVERLAP)
n_both <- sum(overlap_dt$overlap_status == "both", na.rm = TRUE)
n_primary_only <- sum(overlap_dt$overlap_status == "primary_only", na.rm = TRUE)
n_sens_only <- sum(overlap_dt$overlap_status == "sensitivity_only", na.rm = TRUE)
n_union <- n_both + n_primary_only + n_sens_only
jaccard <- if (n_union > 0) n_both / n_union else NA_real_
cat(sprintf(" both / primary_only / sens_only / Jaccard: %d / %d / %d / %.3f\n\n",
            n_both, n_primary_only, n_sens_only, jaccard))

# 7. Volcano (Fig5A)
cat(">>> STEP 7: Fig5A volcano (primary)\n")
vp <- copy(deg_p_dt)
vp[, neg_log10_padj := -log10(pmax(p_adj_BH, .Machine$double.xmin))]
vp[, color_grp := fifelse(direction == "C1>C2", "C1>C2 (up)",
                  fifelse(direction == "C2>C1", "C2>C1 (down)", "ns"))]
vp[, color_grp := factor(color_grp, levels = c("C1>C2 (up)","C2>C1 (down)","ns"))]
top_label_p <- vp[sig == TRUE]
setorder(top_label_p, p_adj_BH)
top_label_p <- head(top_label_p, 10)

p_volcano <- ggplot(vp, aes(x = log2FC_C1minusC2, y = neg_log10_padj, color = color_grp)) +
  geom_point(alpha = 0.5, size = 0.7) +
  geom_vline(xintercept = c(-LOGFC_CUT, LOGFC_CUT), linetype = "dashed", linewidth = 0.3, color = "grey50") +
  geom_hline(yintercept = -log10(ADJ_P_CUT), linetype = "dashed", linewidth = 0.3, color = "grey50") +
  geom_text(data = top_label_p, aes(label = gene), size = 2.8, fontface = "bold",
            vjust = -0.5, show.legend = FALSE) +
  scale_color_manual(values = c("C1>C2 (up)" = unname(CLUSTER_COL["C1"]),
                                 "C2>C1 (down)" = unname(CLUSTER_COL["C2"]),
                                 "ns" = "grey70"),
                     name = NULL) +
  labs(x = expression(log[2]*" fold change (C1 - C2)"),
       y = expression(-log[10]*"(BH-adjusted p-value)")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")
ggsave(paste0(FIG_VOLCANO, ".pdf"), plot = p_volcano, width = 9, height = 7, device = cairo_pdf)
ggsave(paste0(FIG_VOLCANO, ".png"), plot = p_volcano, width = 9, height = 7, dpi = 600, device = ragg::agg_png)
cat(sprintf(" [OK] %s.{pdf,png}\n", FIG_VOLCANO))

# 8. Heatmap (Fig5B)
cat("\n>>> STEP 8: Fig5B heatmap (top 25 up + 25 down by |log2FC|)\n")
sig_dt <- deg_p_dt[sig == TRUE]
top_up <- sig_dt[direction == "C1>C2"][order(-log2FC_C1minusC2)][1:min(25, .N)]
top_down <- sig_dt[direction == "C2>C1"][order( log2FC_C1minusC2)][1:min(25, .N)]
heatmap_genes <- c(top_up$gene, top_down$gene)
cat(sprintf(" Top up (C1>C2) : %d (req 25)\n", nrow(top_up)))
cat(sprintf(" Top down (C2>C1) : %d (req 25)\n", nrow(top_down)))

heat_mat <- expr_filt[heatmap_genes, , drop = FALSE]
col_ord <- order(pheno$cluster, pheno$histology)
heat_mat <- heat_mat[, col_ord]
heat_z <- t(scale(t(heat_mat), center = TRUE, scale = TRUE))

col_anno <- HeatmapAnnotation(
  Cluster = pheno$cluster[col_ord],
  Histology = pheno$histology[col_ord],
  col = list(Cluster = CLUSTER_COL,
             Histology = HIST_COL[levels(pheno$histology)]),
  annotation_name_side = "left",
  annotation_legend_param = list(Cluster = list(title = "Cluster"),
                                  Histology = list(title = "Histology"))
)
row_split <- factor(ifelse(rownames(heat_z) %in% top_up$gene, "C1>C2 (up)", "C2>C1 (down)"),
                    levels = c("C1>C2 (up)","C2>C1 (down)"))
heat_col <- circlize::colorRamp2(c(-2, 0, 2), c("#0571B0", "white", "#CA0020"))
ht <- Heatmap(heat_z,
              name = "z-score", col = heat_col,
              top_annotation = col_anno,
              cluster_columns = FALSE, cluster_rows = TRUE,
              row_split = row_split,
              show_column_names = FALSE,
              row_names_side = "left", row_names_gp = gpar(fontsize = 7),
              row_title_gp = gpar(fontsize = 9, fontface = "bold"),
              column_title = NULL,
              heatmap_legend_param = list(title = "Row z-score"))
W <- 12; H <- 11
cairo_pdf(paste0(FIG_HEATMAP, ".pdf"), width = W, height = H); draw(ht, merge_legend = TRUE); dev.off()
ragg::agg_png(paste0(FIG_HEATMAP, ".png"), width = W, height = H, units = "in", res = 600)
draw(ht, merge_legend = TRUE); dev.off()
cat(sprintf(" [OK] %s.{pdf,png}\n", FIG_HEATMAP))

# 9. clusterProfiler GO/KEGG (combined, primary sig)
cat("\n>>> STEP 9: clusterProfiler GO BP/MF/CC + KEGG (primary sig DEG combined)\n")
sig_genes <- deg_p_dt[sig == TRUE, gene]
gene_map <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat(sprintf(" Input gene list : %d (primary sig)\n", length(sig_genes)))
cat(sprintf(" ENTREZID mapped : %d / %d\n", nrow(gene_map), length(sig_genes)))

run_enrichGO <- function(ont) {
  enrichGO(gene = gene_map$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
           ont = ont, pAdjustMethod = "BH",
           pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
}

cat(" enrichGO BP ...\n"); ego_bp <- run_enrichGO("BP")
cat(" enrichGO MF ...\n"); ego_mf <- run_enrichGO("MF")
cat(" enrichGO CC ...\n"); ego_cc <- run_enrichGO("CC")

cat(" enrichKEGG (with up to 3 retries) ...\n")
ekegg <- NULL
for (try_i in 1:3) {
  ekegg <- tryCatch({
    res <- enrichKEGG(gene = gene_map$ENTREZID, organism = "hsa", keyType = "kegg",
                      pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
    if (!is.null(res)) res <- setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    res
  }, error = function(e) {
    cat(sprintf(" [warn] KEGG attempt %d failed: %s\n", try_i, conditionMessage(e)))
    Sys.sleep(5); NULL
  })
  if (!is.null(ekegg)) { cat(sprintf(" KEGG attempt %d OK.\n", try_i)); break }
}

write_enrich <- function(obj, out_file, label) {
  if (is.null(obj)) { cat(sprintf(" [SKIP] %s NULL\n", basename(out_file))); return(invisible(NULL)) }
  df <- as.data.frame(obj)
  if (nrow(df) == 0) { fwrite(df, out_file); cat(sprintf(" [WARN] %s : 0 rows\n", basename(out_file))); return(invisible(NULL)) }
  df$significant_BH_05 <- df$p.adjust < ADJ_P_CUT
  fwrite(df, out_file)
  cat(sprintf(" [OK] %s : %d rows (%d sig)\n", basename(out_file), nrow(df), sum(df$significant_BH_05)))
  invisible(df)
}
df_bp <- write_enrich(ego_bp, TAB_GO_BP, "BP")
df_mf <- write_enrich(ego_mf, TAB_GO_MF, "MF")
df_cc <- write_enrich(ego_cc, TAB_GO_CC, "CC")
df_kegg <- write_enrich(ekegg, TAB_KEGG, "KEGG")

make_dotplot <- function(obj, title_str, fig_stem) {
  if (is.null(obj)) return(invisible(NULL))
  full_df <- as.data.frame(obj)
  sig_df <- full_df[!is.na(full_df$p.adjust) & full_df$p.adjust < ADJ_P_CUT, , drop = FALSE]
  if (nrow(sig_df) == 0) return(invisible(NULL))
  obj_sig <- obj; obj_sig@result <- sig_df
  show_n <- min(TOP_N_GO, nrow(sig_df))
  pp <- dotplot(obj_sig, showCategory = show_n) +
        theme_bw(base_size = 10)
  ggsave(paste0(fig_stem, ".pdf"), plot = pp, width = 7, height = 6, device = cairo_pdf)
  ggsave(paste0(fig_stem, ".png"), plot = pp, width = 7, height = 6, dpi = 600, device = ragg::agg_png)
  cat(sprintf(" [OK] %s.{pdf,png} (top %d of %d)\n", basename(fig_stem), show_n, nrow(sig_df)))
}
make_dotplot(ego_bp, "DEG GO Biological Process", FIG_GO_BP)
make_dotplot(ego_mf, "DEG GO Molecular Function", FIG_GO_MF)
make_dotplot(ego_cc, "DEG GO Cellular Component", FIG_GO_CC)
make_dotplot(ekegg, "DEG KEGG Pathway", FIG_KEGG)

# 10. Summary
cat("\n========================================================================\n")
cat(" Stage 9 (biotype filter) Summary Report \n")
cat("========================================================================\n")
cat(sprintf(" N samples : %d (C1=%d, C2=%d)\n",
            ncol(expr_filt), sum(pheno$cluster == "C1"), sum(pheno$cluster == "C2")))
cat(sprintf(" Pre-filter cascade:\n"))
cat(sprintf(" Input : %d genes\n", n_genes_in))
cat(sprintf(" after zero filter : %d\n", n_after_zero))
cat(sprintf(" after biotype filter: %d (protein_coding only)\n", n_after_pc))
cat(sprintf("\n PRIMARY (~cluster) : %d sig DEG\n", n_primary_sig))
cat(sprintf(" SENSITIVITY (+histology) : %d sig DEG\n", n_sens_sig))
cat(sprintf(" Overlap both : %d\n", n_both))
cat(sprintf(" Jaccard : %.3f\n", jaccard))
cat(sprintf("\n Top 5 PRIMARY sig DEG (smallest BH adj.p):\n"))
print(head(deg_p_dt[sig == TRUE, .(gene, log2FC_C1minusC2, p_adj_BH, direction)], 5))
cat("\n========================================================================\n")

# sessionInfo
sess_file <- file.path(LOGDIR, "sessionInfo/sessionInfo_stage09.txt")
sink(sess_file); cat(sprintf("Session info captured: %s\n\n", Sys.time())); print(sessionInfo()); sink()
cat(sprintf(">>> sessionInfo saved: %s\n", sess_file))

cat(sprintf("\n=== Stage 9 Pipeline done (%s) ===\n", Sys.time()))
sink(type = "message"); sink(); close(log_con)

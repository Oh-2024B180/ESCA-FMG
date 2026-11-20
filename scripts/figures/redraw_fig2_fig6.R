###############################################################################
## Figures 2 and 6 redrawn at the size they are printed
## (16.5 cm wide) from the saved result tables — no statistics are recomputed.
##
##   Rscript redraw_fig2_fig6.R            (R 4.4.1; plotting only)
##
## Inputs (read-only): tables/Tab3_FMG_GO_{BP,CC,MF}.csv, Tab3_FMG_KEGG.csv (Figure 2),
##   tables/Tab9_DEG_GO_{BP,CC,MF}.csv, Tab9_DEG_KEGG.csv, Tab9_DEG_full_C1_vs_C2.csv (Figure 6),
##   data/expr_ESCA_tumor.csv + tables/Tab7_cluster_assignment.csv (heatmap, Figure S1 draft)
## Rules carried over from the published panels:
##   GO/KEGG: top 10 terms by BH-adjusted p, shown by GeneRatio, size = Count, colour = adjusted p (red -> blue,
##            log scale) — figures/remake/v3_fig1BC_5CD_dotplot.R. NEW: only terms with adjusted p < 0.05
##            (Methods 2.4); in Figure 2 this removes six non-significant CC terms (p.adjust 0.055-0.097).
##   Volcano: |log2FC| >= 1 and BH p < 0.05, top 10 labels per direction by adjusted p, top 3 framed,
##            y capped at 20 — scripts/helpers/09_fix_Fig5A_volcano.R. Gene labels in italics.
##   Heatmap: top 25 up / 25 down significant DEG by log2FC, row z-scores, columns by cluster then histology,
##            rows clustered within each direction — scripts/09_DEG_between_clusters.R step 8.
## Text: 7-7.5 pt at the printed size, panel tags 10 pt bold; adjusted p written as powers of ten.
###############################################################################
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(stringr); library(patchwork); library(ragg)
  library(ggrepel); library(ComplexHeatmap); library(circlize); library(grid)
})
set.seed(1262118388)
E <- getwd()   # repository root (see README.md)
OUT <- file.path(E, "figures", "final")
dir.create(OUT, showWarnings = FALSE)
W_IN <- 16.5 / 2.54                       # printed width
RES <- 4000 / W_IN                         # 4000 px across, as the assembled figures

theme_fig <- theme_bw(base_size = 7, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
        panel.border = element_rect(colour = "black", linewidth = 0.4, fill = NA),
        axis.text = element_text(size = 7, colour = "black"),
        axis.title = element_text(size = 7.5, colour = "black"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 7.5),
        legend.key.size = unit(3.2, "mm"),
        plot.tag = element_text(size = 10, face = "bold", family = "Arial"),
        plot.margin = margin(3, 6, 3, 3))

save_both <- function(p, stem, h_in) {
  ggsave(file.path(OUT, paste0(stem, ".pdf")), p, width = W_IN, height = h_in, device = cairo_pdf)
  agg_png(file.path(OUT, paste0(stem, ".png")), width = W_IN, height = h_in, units = "in", res = RES)
  print(p); invisible(dev.off())
  cat(sprintf("  %s: %.1f x %.1f cm\n", stem, W_IN * 2.54, h_in * 2.54))
}

# ------------------------------------------------------------------ GO / KEGG dot plots
parse_gr <- function(x) vapply(strsplit(as.character(x), "/"), function(p) as.numeric(p[1]) / as.numeric(p[2]), 1)
read_top <- function(f, top_n = 10) {
  dt <- fread(file.path(E, "tables", f))
  setorder(dt, p.adjust)
  n_all <- nrow(dt[seq_len(min(top_n, .N))])
  dt <- dt[p.adjust < 0.05][seq_len(min(top_n, sum(dt$p.adjust < 0.05)))]
  cat(sprintf("  %-24s %2d terms (top %d by adjusted p, of which %d with p < 0.05)\n", f, nrow(dt), n_all, nrow(dt)))
  dt[, gr := parse_gr(GeneRatio)]
  dt
}
pow10 <- function(x) parse(text = sprintf("10^%d", as.integer(round(log10(x)))))
dot_panel <- function(dt, clim, slim, wrap = 36, horiz = FALSE) {
  dt <- copy(dt)
  dt[, lab := str_wrap(Description, wrap)]
  setorder(dt, gr)
  dt[, lab := factor(lab, levels = unique(lab))]
  lo <- ceiling(log10(clim[1])); hi <- floor(log10(clim[2]))
  ex <- lo:hi
  if (length(ex) > 5) { ex <- pretty(c(lo, hi), n = 4); ex <- ex[ex >= lo & ex <= hi] }
  br <- 10^ex
  ggplot(dt, aes(gr, lab, size = Count, colour = p.adjust)) +
    geom_point() +
    scale_colour_gradient(low = "#E41A1C", high = "#377EB8", trans = "log10", limits = clim, breaks = br,
                          labels = pow10, name = expression("Adjusted" ~ italic(p)),
                          guide = if (horiz) guide_colourbar(direction = "horizontal", barwidth = unit(32, "mm"),
                                                             barheight = unit(2.5, "mm"), title.vjust = 0.8)
                                  else guide_colourbar(reverse = TRUE, barheight = unit(18, "mm"), barwidth = unit(3, "mm"))) +
    scale_size_continuous(limits = slim, range = c(1.2, 4.2), name = "Count") +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    labs(x = "GeneRatio", y = NULL) + theme_fig
}
go_set <- function(files, horiz = FALSE) {
  dts <- lapply(files, read_top)
  allp <- unlist(lapply(dts, `[[`, "p.adjust")); allc <- unlist(lapply(dts, `[[`, "Count"))
  clim <- c(10^floor(log10(min(allp))), 0.05)
  slim <- range(allc)
  lapply(dts, dot_panel, clim = clim, slim = slim, horiz = horiz)
}

cat("Figure 2 (FMG enrichment)\n")
f2 <- go_set(c("Tab3_FMG_GO_BP.csv", "Tab3_FMG_GO_CC.csv", "Tab3_FMG_GO_MF.csv", "Tab3_FMG_KEGG.csv"))
fig2 <- (f2[[1]] | f2[[2]]) / (f2[[3]] | f2[[4]]) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)"))) &
  theme(legend.position = "right")
save_both(fig2, "Figure2_draft", h_in = 13.0 / 2.54)

# ------------------------------------------------------------------ volcano
cat("Figure 6 (DEG between clusters)\n")
deg <- fread(file.path(E, "tables/Tab9_DEG_full_C1_vs_C2.csv"))
CUT_LFC <- 1.0; CUT_P <- 0.05
deg[, dir_class := fifelse(p_adj_BH < CUT_P & log2FC_C1minusC2 >= CUT_LFC, "up_C1",
                    fifelse(p_adj_BH < CUT_P & log2FC_C1minusC2 <= -CUT_LFC, "up_C2", "ns"))]
deg[, neglog10p := -log10(p_adj_BH)]
top_c1 <- deg[dir_class == "up_C1"][order(p_adj_BH)][1:10]
top_c2 <- deg[dir_class == "up_C2"][order(p_adj_BH)][1:10]
top_all <- rbindlist(list(top_c1, top_c2))
top_all[, top_rank := seq_len(.N), by = dir_class]
top_all[, framed := top_rank <= 3]
top_all[, fill_key := fifelse(framed, dir_class, "plain")]
top_all[, face := fifelse(framed, "bold.italic", "italic")]
cat(sprintf("  DEG up in C1 %d, up in C2 %d; framed: %s\n", sum(deg$dir_class == "up_C1"), sum(deg$dir_class == "up_C2"),
            paste(top_all[framed == TRUE, gene], collapse = ", ")))
y_cap <- min(max(deg$neglog10p, na.rm = TRUE), 20)
COL <- c(up_C1 = "#1B9E77", up_C2 = "#D95F02", ns = "grey75")
volcano <- ggplot(deg, aes(log2FC_C1minusC2, neglog10p, colour = dir_class)) +
  geom_vline(xintercept = c(-CUT_LFC, CUT_LFC), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = -log10(CUT_P), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_point(size = 0.35, alpha = 0.55) +
  # one label layer so framed and plain labels repel each other (two layers let boxes cover plain names);
  # framed = top 3 per direction, filled box and bold, as in the published panel
  geom_label_repel(data = top_all, aes(label = gene, fill = fill_key, fontface = face), size = 7 / .pt,
                   family = "Arial", label.size = NA, label.padding = unit(0.12, "lines"), max.overlaps = Inf,
                   force = 10, force_pull = 0.5, box.padding = 0.35, point.padding = 0.2, segment.alpha = 0.5,
                   segment.size = 0.2, min.segment.length = 0, seed = 1262118388, show.legend = FALSE) +
  scale_colour_manual(values = COL, labels = c(up_C1 = "C1 > C2 (up)", up_C2 = "C2 > C1 (down)", ns = "ns"),
                      breaks = c("up_C1", "up_C2", "ns"), name = NULL) +
  scale_fill_manual(values = c(up_C1 = "#A8D8C9", up_C2 = "#F8C29B", plain = NA), na.value = NA, guide = "none") +
  coord_cartesian(ylim = c(0, y_cap + 3)) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = expression(log[2] ~ "fold change (C1" ~ "\u2212" ~ "C2)"), y = expression(-log[10] ~ "(BH-adjusted" ~ italic(p) * ")")) +
  labs(tag = "(a)") +
  theme_fig + theme(legend.position = "top", legend.direction = "horizontal", legend.margin = margin(0, 0, -4, 0),
                    legend.key = element_blank(), legend.key.size = unit(2.5, "mm"))

f6 <- go_set(c("Tab9_DEG_GO_BP.csv", "Tab9_DEG_GO_CC.csv", "Tab9_DEG_GO_MF.csv", "Tab9_DEG_KEGG.csv"), horiz = TRUE)
go6 <- (f6[[1]] | f6[[2]]) / (f6[[3]] | f6[[4]]) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")
fig6 <- wrap_elements(full = volcano, ignore_tag = TRUE) / go6 + plot_layout(heights = c(6.4, 13.4)) +
  plot_annotation(tag_levels = list(c("(b)", "(c)", "(d)", "(e)")))
save_both(fig6, "Figure6_draft", h_in = 20.3 / 2.54)

# ------------------------------------------------------------------ heatmap (proposed Figure S1)
cat("Heatmap (proposed Figure S1)\n")
expr <- fread(file.path(E, "data/expr_ESCA_tumor.csv"))
cl <- fread(file.path(E, "tables/Tab7_cluster_assignment.csv"))
m <- as.matrix(expr[, -1]); rownames(m) <- expr[[1]]
m <- m[, cl$sampleID]
sig <- deg[sig == TRUE]
top_up <- sig[direction == "C1>C2"][order(-log2FC_C1minusC2)][1:25]
top_dn <- sig[direction == "C2>C1"][order(log2FC_C1minusC2)][1:25]
genes <- c(top_up$gene, top_dn$gene)
cl[, cluster := factor(cluster, levels = c("C1", "C2"))]
cl[, histology := droplevels(factor(histology, levels = c("ESCC", "EAC", "Other")))]
ord <- order(cl$cluster, cl$histology)
hz <- t(scale(t(m[genes, ord])))
anno <- HeatmapAnnotation(Cluster = cl$cluster[ord], Histology = cl$histology[ord],
                          col = list(Cluster = c(C1 = "#1B9E77", C2 = "#D95F02"),
                                     Histology = c(ESCC = "#7570B3", EAC = "#E7298A", Other = "#66A61E")[levels(cl$histology)]),
                          annotation_name_side = "left",
                          annotation_name_gp = gpar(fontsize = 7.5, fontfamily = "Arial"),
                          simple_anno_size = unit(3, "mm"),
                          annotation_legend_param = list(title_gp = gpar(fontsize = 7.5, fontfamily = "Arial"),
                                                         labels_gp = gpar(fontsize = 7, fontfamily = "Arial")))
split <- factor(ifelse(genes %in% top_up$gene, "C1 > C2 (up)", "C2 > C1 (down)"), levels = c("C1 > C2 (up)", "C2 > C1 (down)"))
ht <- Heatmap(hz, name = "Row z-score", col = colorRamp2(c(-2, 0, 2), c("#0571B0", "white", "#CA0020")),
              top_annotation = anno, cluster_columns = FALSE, cluster_rows = TRUE, row_split = split,
              show_column_names = FALSE, row_names_side = "left",
              row_names_gp = gpar(fontsize = 7, fontface = "italic"), row_names_max_width = unit(2.2, "cm"),
              row_title_gp = gpar(fontsize = 7.5, fontface = "bold", fontfamily = "Arial"),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 7.5, fontfamily = "Arial"),
                                          labels_gp = gpar(fontsize = 7, fontfamily = "Arial")))
for (dev in c("pdf", "png")) {
  f <- file.path(OUT, paste0("FigureS1_heatmap_draft.", dev))
  if (dev == "pdf") cairo_pdf(f, width = W_IN, height = 16 / 2.54) else
    agg_png(f, width = W_IN, height = 16 / 2.54, units = "in", res = RES)
  draw(ht, merge_legend = TRUE); invisible(dev.off())
}
cat(sprintf("  FigureS1_heatmap_draft: 16.5 x 16.0 cm, %d genes\n", length(genes)))

###############################################################################
## Stage 4: TCGA-ESCA Mutation OncoPrint (18 Primary FMG)
##
## Reference:
## Liu et al. 2024 (Heliyon) Figure 2A convention -- per-sample alteration
## landscape of FMG with histology / clinical annotation.
##
## Inputs:
## data/ESCA_mc3.csv (TCGA MC3 mutation calls)
## data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv (TCGA clinical matrix)
## tables/Tab1_2_FMG_list_primary.csv (18 Primary FMG)
##
## Sample filter:
## clinical sample_type == "Primary Tumor"
## intersect with samples that appear in MC3
##
## Functional consequence filter:
## Missense_Mutation, Nonsense_Mutation, Frame_Shift_Del, Frame_Shift_Ins,
## Splice_Site, In_Frame_Del, In_Frame_Ins, Translation_Start_Site,
## Nonstop_Mutation
## (Silent, 3'UTR, 5'UTR, Intron, Flank, RNA -> excluded)
##
## Outputs (tables/):
## Tab4_TCGA_mutation_frequency.csv
## columns: gene, n_mut, N_total, prevalence, percent,
## n_mut_ESCC, N_ESCC, percent_ESCC,
## n_mut_EAC, N_EAC, percent_EAC
##
## Outputs (figures/, PDF + PNG @ 600 dpi -- CLAUDE.md fig ):
## Fig2A_TCGA_oncoprint.{pdf,png}
##
## Random seed: 1262118388 (Liu et al. 2024)
##
## Notes:
## - Histology mapping:
## "Esophagus Squamous Cell Carcinoma" -> "ESCC"
## "Esophagus Adenocarcinoma, NOS" -> "EAC"
## else -> "Other"
## - oncoPrint is saved with cairo_pdf() / ragg::agg_png() (not ggsave).
## - Layout (option A; ComplexHeatmap needs percentages and row names on different sides):
## top : anno_oncoprint_barplot (per-sample mutation count)
## bottom : Histology annotation (ESCC vs EAC)
## right : percent + anno_oncoprint_barplot (per-gene n samples)
## left : gene name
## Note: ComplexHeatmap 2.22.0 does not allow pct_side == row_names_side
## (Error: Percent values and row names should be at different side)
###############################################################################

# 0. Working directory + output directory
# run from the repository root (see README.md)

INPUT_MC3 <- "data/ESCA_mc3.csv"
INPUT_CLI <- "data/TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv"
INPUT_FMG <- "tables/Tab1_2_FMG_list_primary.csv"

TBLDIR <- "tables"
FIGDIR <- "figures"
LOGDIR <- "logs"

for (d in c(TBLDIR, FIGDIR, LOGDIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

FIG_STEM <- "Fig2A_TCGA_oncoprint"
TBL_OUT <- file.path(TBLDIR, "Tab4_TCGA_mutation_frequency.csv")

FUNCTIONAL_EFFECTS <- c(
  "Missense_Mutation", "Nonsense_Mutation",
  "Frame_Shift_Del", "Frame_Shift_Ins",
  "Splice_Site",
  "In_Frame_Del", "In_Frame_Ins",
  "Translation_Start_Site", "Nonstop_Mutation"
)

# 0.5. Reproducibility
set.seed(1262118388)

# 0.6. Logging start
log_file <- file.path(LOGDIR, "04_diagnostics.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

cat(sprintf("=== Stage 4 Pipeline start (%s) ===\n\n", Sys.time()))

# 1. load
suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
  library(ragg)
})

cat(sprintf("ComplexHeatmap version : %s\n", packageVersion("ComplexHeatmap")))
cat(sprintf("circlize version : %s\n\n", packageVersion("circlize")))

# 2. Load input data
cat(">>> STEP 2: Load inputs\n")
mc3 <- fread(INPUT_MC3)
cli <- fread(INPUT_CLI)
fmg <- fread(INPUT_FMG)
fmg_genes <- fmg$gene

cat(sprintf(" MC3 rows : %d\n", nrow(mc3)))
cat(sprintf(" Clinical rows : %d\n", nrow(cli)))
cat(sprintf(" FMG count : %d\n\n", length(fmg_genes)))

# 3. Sample filter: Primary Tumor & MC3 intersection
cat(">>> STEP 3: Sample filtering\n")
primary_tumor_samples <- cli[sample_type == "Primary Tumor", sampleID]
cat(sprintf(" Clinical Primary Tumor : %d\n", length(primary_tumor_samples)))

mc3_samples <- unique(mc3$sample)
cat(sprintf(" MC3 unique samples : %d\n", length(mc3_samples)))

samples_keep <- intersect(primary_tumor_samples, mc3_samples)
N_TOTAL <- length(samples_keep)
cat(sprintf(" intersection (analysis) : %d\n", N_TOTAL))

# Histology mapping
cli_keep <- cli[sampleID %in% samples_keep, .(sampleID, histological_type)]
cli_keep[, histology := fifelse(
  histological_type == "Esophagus Squamous Cell Carcinoma", "ESCC",
  fifelse(histological_type == "Esophagus Adenocarcinoma, NOS", "EAC", "Other")
)]
cat(" Histology breakdown:\n")
print(cli_keep[, .N, by=histology][order(-N)])
cat("\n")

# 4. Mutation matrix
cat(">>> STEP 4: Build mutation matrix (sample x gene)\n")
mc3_filt <- mc3[sample %in% samples_keep & gene %in% fmg_genes & effect %in% FUNCTIONAL_EFFECTS]
cat(sprintf(" Filtered MC3 rows : %d\n", nrow(mc3_filt)))

# OncoPrint mapping
cat_map <- function(eff) {
  out <- character(length(eff))
  out[eff == "Missense_Mutation"] <- "Missense"
  out[eff == "Nonsense_Mutation"] <- "Nonsense"
  out[eff %in% c("Frame_Shift_Del","Frame_Shift_Ins")] <- "Frameshift"
  out[eff %in% c("In_Frame_Del","In_Frame_Ins")] <- "Inframe_indel"
  out[eff == "Splice_Site"] <- "Splice"
  out[eff %in% c("Translation_Start_Site","Nonstop_Mutation")] <- "Other_functional"
  out
}
mc3_filt[, alt_type := cat_map(effect)]
cat(" Alteration categories:\n")
print(mc3_filt[, .N, by=alt_type][order(-N)])
cat("\n")

# Per (sample, gene) -> alteration string ";"-joined
mat_long <- mc3_filt[, .(alt = paste(unique(alt_type), collapse = ";")),
                     by = .(sample, gene)]

# 18 FMG x N_TOTAL matrix (string), rows = gene, cols = sample
mat <- matrix("", nrow = length(fmg_genes), ncol = N_TOTAL,
              dimnames = list(fmg_genes, samples_keep))
for (i in seq_len(nrow(mat_long))) {
  mat[mat_long$gene[i], mat_long$sample[i]] <- mat_long$alt[i]
}

# 5. tablecompose ( +
cat(">>> STEP 5: Write mutation frequency table\n")
escc_samples <- cli_keep[histology == "ESCC", sampleID]
eac_samples <- cli_keep[histology == "EAC", sampleID]
N_ESCC <- length(escc_samples)
N_EAC <- length(eac_samples)

freq_tbl <- data.table(
  gene = fmg_genes,
  n_mut = vapply(fmg_genes, function(g) sum(mat[g, ] != ""), integer(1)),
  N_total = N_TOTAL,
  n_mut_ESCC = vapply(fmg_genes, function(g) sum(mat[g, escc_samples] != ""), integer(1)),
  N_ESCC = N_ESCC,
  n_mut_EAC = vapply(fmg_genes, function(g) sum(mat[g, eac_samples] != ""), integer(1)),
  N_EAC = N_EAC
)
freq_tbl[, prevalence := n_mut / N_total]
freq_tbl[, percent := round(100 * prevalence, 2)]
freq_tbl[, percent_ESCC := round(100 * n_mut_ESCC / N_ESCC, 2)]
freq_tbl[, percent_EAC := round(100 * n_mut_EAC / N_EAC, 2)]
setorder(freq_tbl, -n_mut)

setcolorder(freq_tbl, c("gene","n_mut","N_total","prevalence","percent",
                        "n_mut_ESCC","N_ESCC","percent_ESCC",
                        "n_mut_EAC","N_EAC","percent_EAC"))
fwrite(freq_tbl, TBL_OUT)
cat(sprintf(" [OK] %s\n", TBL_OUT))
cat(" Top 5 by total mutation count:\n")
print(head(freq_tbl, 5))
cat("\n")

# 6. OncoPrint
cat(">>> STEP 6: Build oncoPrint\n")

# color (ColorBrewer Set1, color-blind-safe)
# : previous closuregeneratealter_funR lazy-eval
# cellOther_functional, brown)apply .
# colours are hard-coded in each function body (avoids closure capture)
col_named_vec <- c(
  Missense = "#377EB8", # blue - type, Nonsense = "#E41A1C", # red - truncating, “
  Frameshift = "#984EA3", # purple - truncating
  Inframe_indel = "#FF7F00", # orange - in-frame indel
  Splice = "#4DAF4A", # green - splice
  Other_functional = "#A65628" # brown - rare (TSS, Nonstop)
)

alter_fun_list <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.92, h*0.92, gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Missense = function(x, y, w, h) {
    grid.rect(x, y, w*0.92, h*0.92, gp = gpar(fill = "#377EB8", col = NA))
  },
  Nonsense = function(x, y, w, h) {
    grid.rect(x, y, w*0.92, h*0.92, gp = gpar(fill = "#E41A1C", col = NA))
  },
  Frameshift = function(x, y, w, h) {
    grid.rect(x, y, w*0.92, h*0.92, gp = gpar(fill = "#984EA3", col = NA))
  },
  Inframe_indel = function(x, y, w, h) {
    grid.rect(x, y, w*0.92, h*0.92, gp = gpar(fill = "#FF7F00", col = NA))
  },
  Splice = function(x, y, w, h) {
    grid.rect(x, y, w*0.92, h*0.92, gp = gpar(fill = "#4DAF4A", col = NA))
  },
  Other_functional = function(x, y, w, h) {
    grid.rect(x, y, w*0.92, h*0.92, gp = gpar(fill = "#A65628", col = NA))
  }
)

# Histology vector (column-aligned)
hist_vec <- cli_keep[match(samples_keep, sampleID), histology]
hist_colors <- c("ESCC" = "#1B9E77", "EAC" = "#D95F02", "Other" = "#7570B3")

# Layout (as in the earlier Figure 2A):
# top : per-sample mutation count barplot
# bottom : Histology
# right : per-gene n samples barplot
# left : gene name + percent
top_anno <- HeatmapAnnotation(
  cbar = anno_oncoprint_barplot(),
  show_annotation_name = FALSE
)
bottom_anno <- HeatmapAnnotation(
  Histology = hist_vec,
  col = list(Histology = hist_colors),
  annotation_name_side = "left",
  annotation_legend_param = list(Histology = list(title = "Histology"))
)
right_anno <- rowAnnotation(
  rbar = anno_oncoprint_barplot(width = unit(2, "cm")),
  show_annotation_name = FALSE
)

ht <- oncoPrint(
  mat,
  alter_fun = alter_fun_list,
  alter_fun_is_vectorized = FALSE,
  col = col_named_vec,
  top_annotation = top_anno,
  bottom_annotation = bottom_anno,
  right_annotation = right_anno,
  row_order = freq_tbl$gene,
  remove_empty_columns = FALSE,
  remove_empty_rows = FALSE,
  pct_side = "right",
  row_names_side = "left",
  column_title = NULL,
  heatmap_legend_param = list(title = "Alteration",
                              at = names(col_named_vec),
                              labels = names(col_named_vec))
)

# 7. Save (PDF cairo + PNG ragg @1200 dpi)
cat(">>> STEP 7: Save oncoprint\n")
pdf_file <- file.path(FIGDIR, paste0(FIG_STEM, ".pdf"))
png_file <- file.path(FIGDIR, paste0(FIG_STEM, ".png"))

W <- 12
H <- 5

cairo_pdf(pdf_file, width = W, height = H)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()
cat(sprintf(" [OK] %s\n", pdf_file))

ragg::agg_png(png_file, width = W, height = H, units = "in", res = 600)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()
cat(sprintf(" [OK] %s\n\n", png_file))

# 8. Summary Report
cat("========================================================================\n")
cat(" Stage 4 TCGA-ESCA OncoPrint Summary Report \n")
cat("========================================================================\n")
cat(sprintf(" Total Primary Tumor (analysis) : %d\n", N_TOTAL))
cat(sprintf(" ESCC : %d\n", N_ESCC))
cat(sprintf(" EAC : %d\n", N_EAC))
cat(sprintf(" Other : %d\n", N_TOTAL - N_ESCC - N_EAC))
cat(sprintf(" FMG analyzed : %d\n", length(fmg_genes)))
cat(sprintf(" Filtered functional MC3 rows : %d\n\n", nrow(mc3_filt)))

cat(" --- 18 FMG mutation frequency (sorted by n_mut) ---\n")
print(freq_tbl)
cat("\n========================================================================\n")

# 9. sessionInfo
sess_file <- file.path(LOGDIR, "sessionInfo/sessionInfo_stage04.txt")
sink(sess_file)
cat(sprintf("Session info captured: %s\n\n", Sys.time()))
print(sessionInfo())
sink()
cat(sprintf("\n>>> sessionInfo saved: %s\n", sess_file))

cat(sprintf("\n=== Stage 4 Pipeline done (%s) ===\n", Sys.time()))

# end
sink(type = "message")
sink()
close(log_con)

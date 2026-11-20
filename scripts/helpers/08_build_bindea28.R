###############################################################################
##  Build canonical Bindea 28 immune signature for Stage 8C ssGSEA
##
##  Sources:
##    1. IOBR R package v2.2.0 signature_collection (27 *_Bindea_et_al keys)
##       - faithful to Bindea 2013 supp Table 1 per content validation
##       - missing only NK CD56bright cells
##    2. NK CD56bright cells (6 genes) extracted directly from
##       Bindea 2013 *Immunity* (PMID 24138885) supp Table 1
##       (cache/bindea_supp.txt, mmc1.pdf pages around Table S1 (3/6))
##
##  Output:
##    data/Bindea_28_signatures.tsv  (long: cell_type, gene_symbol)
##
##  Validation:
##    - 28 canonical cell types (Bindea 2013 supp methods + Table S1)
##    - ~577-585 unique genes (Bindea 2013 published: 577)
###############################################################################

# run from the repository root (see README.md)
suppressPackageStartupMessages({
  library(IOBR)
  library(data.table)
})

OUT_TSV <- "data/Bindea_28_signatures.tsv"

# ── 1. IOBR Bindea signatures (27 cell types) ──────────────────────────────
e <- new.env()
data("signature_collection", package = "IOBR", envir = e)
sc <- e$signature_collection

bin_keys <- grep("Bindea_et_al$", names(sc), value = TRUE)
cat(sprintf("IOBR Bindea signatures: %d cell types\n", length(bin_keys)))

# Map IOBR keys to canonical cell type names (Bindea 2013 supp methods convention)
iobr_to_canonical <- c(
  "B_cells_Bindea_et_al"            = "B cells",
  "T_cells_Bindea_et_al"            = "T cells",
  "T_helper_cells_Bindea_et_al"     = "T helper cells",
  "Tcm_Bindea_et_al"                = "Tcm",
  "Tem_Bindea_et_al"                = "Tem",
  "Th1_cells_Bindea_et_al"          = "Th1 cells",
  "Th2_cells_Bindea_et_al"          = "Th2 cells",
  "TFH_Bindea_et_al"                = "TFH",
  "Th17_cells_Bindea_et_al"         = "Th17 cells",
  "Treg_Bindea_et_al"               = "TReg",
  "CD8_T_cells_Bindea_et_al"        = "CD8 T cells",
  "Tgd_Bindea_et_al"                = "Tgd",
  "Cytotoxic_cells_Bindea_et_al"    = "Cytotoxic cells",
  "NK_cells_Bindea_et_al"           = "NK cells",
  "NK_CD56dim_cells_Bindea_et_al"   = "NK CD56dim cells",
  "DC_Bindea_et_al"                 = "DC",
  "iDC_Bindea_et_al"                = "iDC",
  "aDC_Bindea_et_al"                = "aDC",
  "pDC_Bindea_et_al"                = "pDC",
  "Eosinophils_Bindea_et_al"        = "Eosinophils",
  "Macrophages_Bindea_et_al"        = "Macrophages",
  "Mast_cells_Bindea_et_al"         = "Mast cells",
  "Neutrophils_Bindea_et_al"        = "Neutrophils",
  "SW480_cancer_cells_Bindea_et_al" = "SW480 cancer cells",
  "Normal_mucosa_Bindea_et_al"      = "Normal mucosa",
  "Blood_vessels_Bindea_et_al"      = "Blood vessels",
  "Lymph_vessels_Bindea_et_al"      = "Lymph vessels"
)
stopifnot(setequal(names(iobr_to_canonical), bin_keys))

# Build long-format from IOBR
iobr_records <- lapply(names(iobr_to_canonical), function(k) {
  data.table(
    cell_type   = iobr_to_canonical[[k]],
    gene_symbol = sc[[k]],
    source      = "IOBR_signature_collection_v2.2.0"
  )
})
iobr_dt <- rbindlist(iobr_records)
cat(sprintf("IOBR-derived rows: %d (cell_type, gene) pairs\n", nrow(iobr_dt)))

# ── 2. NK CD56bright cells (6 genes from Bindea 2013 PDF Table S1) ─────────
# Extracted from cache/bindea_supp.txt parsing of mmc1.pdf
# (NK CD56bright cells section in Table S1 (3/6), page 6 of supp PDF)
nk_cd56bright <- data.table(
  cell_type   = "NK CD56bright cells",
  gene_symbol = c("NCR1", "KIR3DL2", "KIR3DS1", "KIR2DL3", "GZMB", "SPON2"),
  source      = "Bindea_2013_PMID24138885_supp_TableS1"
)
cat(sprintf("NK CD56bright cells (PDF-extracted): %d genes\n", nrow(nk_cd56bright)))

# ── 3. Merge + validate ─────────────────────────────────────────────────────
final_dt <- rbindlist(list(iobr_dt, nk_cd56bright))

# Canonical 28 cell types (Bindea 2013 supp methods)
canonical_28 <- c(
  "B cells", "T cells", "T helper cells",
  "Tcm", "Tem",
  "Th1 cells", "Th2 cells",
  "TFH", "Th17 cells", "TReg",
  "CD8 T cells", "Tgd", "Cytotoxic cells",
  "NK cells", "NK CD56dim cells", "NK CD56bright cells",
  "DC", "iDC", "aDC", "pDC",
  "Eosinophils", "Macrophages", "Mast cells", "Neutrophils",
  "SW480 cancer cells", "Normal mucosa",
  "Blood vessels", "Lymph vessels"
)

obs_cell_types <- unique(final_dt$cell_type)
cat(sprintf("\n=== Validation ===\n"))
cat(sprintf("Cell types: %d / %d canonical\n", length(obs_cell_types), length(canonical_28)))
cat(sprintf("Missing canonical: %s\n",
            ifelse(length(setdiff(canonical_28, obs_cell_types)) == 0, "NONE",
                   paste(setdiff(canonical_28, obs_cell_types), collapse = ", "))))
cat(sprintf("Extra (non-canonical): %s\n",
            ifelse(length(setdiff(obs_cell_types, canonical_28)) == 0, "NONE",
                   paste(setdiff(obs_cell_types, canonical_28), collapse = ", "))))

cat(sprintf("\nTotal unique genes: %d (Bindea 2013 published: 577)\n",
            length(unique(final_dt$gene_symbol))))

cat("\n=== Cell type x gene count (sorted by count desc) ===\n")
ct_summary <- final_dt[, .(n_genes = uniqueN(gene_symbol)), by = cell_type][order(-n_genes)]
print(ct_summary)

# ── 4. Order: canonical 28 cell types order, then alphabetical gene ─────────
final_dt[, cell_type := factor(cell_type, levels = canonical_28)]
setorder(final_dt, cell_type, gene_symbol)
final_dt[, cell_type := as.character(cell_type)]

# ── 5. Save ────────────────────────────────────────────────────────────────
fwrite(final_dt, OUT_TSV, sep = "\t")
cat(sprintf("\n[OK] Saved: %s\n", OUT_TSV))
cat(sprintf("     %d rows  (%d unique genes across %d cell types)\n",
            nrow(final_dt), length(unique(final_dt$gene_symbol)),
            length(unique(final_dt$cell_type))))

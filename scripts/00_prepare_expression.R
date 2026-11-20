###############################################################################
##  00_prepare_expression.R  —  expression matrices read by stages 04c and 06-14
##
##  Extracts the samples the analysis needs from the UCSC Xena "TCGA TARGET GTEx" expression matrix (Toil
##  re-processing) and writes them as separate files. This is the preparation script of the original analysis
##  (April 2026); for this release its comments and messages were translated, the working-directory line was removed
##  and the files are read from and written to data/. Its three outputs used later are byte-identical to the
##  analysis files.
##
##  Inputs (data/; see README.md for the Xena datasets):
##    1) TcgaTargetGtex_RSEM_Hugo_norm_count            (tab-separated, about 8 GB decompressed)
##    2) TcgaTargetGTEX_phenotype.csv                   (phenotype table of the same Xena hub, saved as CSV)
##    3) TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv    (TCGA-ESCA clinical matrix, saved as CSV)
##
##  Outputs (data/):
##    [1] expr_ESCA_tumor.csv              TCGA-ESCA primary tumors (sample IDs ending in -01) found in the
##                                         matrix: 181 in the analysis
##    [2] expr_ESCA_normal_TCGA.csv        TCGA-ESCA solid tissue normal samples (-11); not read by later stages
##    [3] expr_GTEx_esophagus.csv          GTEx esophagus samples (detailed_category contains "Esophagus"): 655
##    [4] GTEx_esophagus_subtype_map.csv   sample, subtype (Mucosa 273 / Muscularis 245 / Gastroesophageal
##                                         Junction 137), sex
##
##  Values are log2(norm_count + 1) as distributed by Xena; no further transformation.
###############################################################################

# ── 0. Packages ──────────────────────────────────────────────────────────────
if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table", repos = "https://cloud.r-project.org")

library(data.table)

# run from the repository root (see README.md)
DATA <- "data"

# ── 1. Sample IDs ────────────────────────────────────────────────────────────
cat(strrep("=", 70), "\n")
cat(">>> STEP 1: sample IDs\n")
cat(strrep("=", 70), "\n\n")

# ─── 1-1. TCGA-ESCA samples (clinical matrix) ──────────────────────────────
clin <- fread(file.path(DATA, "TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv"))

# Primary tumor (-01)
esca_tumor_ids <- clin[grepl("-01$", sampleID), sampleID]
cat(sprintf("  TCGA ESCA Tumor (-01):      %d samples\n", length(esca_tumor_ids)))

# Solid tissue normal (-11)
esca_normal_ids <- clin[grepl("-11$", sampleID), sampleID]
cat(sprintf("  TCGA ESCA Normal (-11):     %d samples\n", length(esca_normal_ids)))

# ─── 1-2. GTEx esophagus samples (phenotype table) ─────────────────────────
pheno <- fread(file.path(DATA, "TcgaTargetGTEX_phenotype.csv"))

gtex_eso <- pheno[
  grepl("Esophagus", detailed_category, ignore.case = TRUE) &
    `_study` == "GTEX"
]

cat(sprintf("  GTEx Esophagus Normal:      %d samples\n", nrow(gtex_eso)))
cat("    Subtypes:\n")
print(gtex_eso[, .N, by = detailed_category][order(-N)])

gtex_eso_ids <- gtex_eso$sample

# ── 2. Expression matrix ─────────────────────────────────────────────────────
cat("\n", strrep("=", 70), "\n")
cat(">>> STEP 2: reading the expression matrix (this takes a while)\n")
cat(strrep("=", 70), "\n\n")

rsem_file <- file.path(DATA, "TcgaTargetGtex_RSEM_Hugo_norm_count")

if (!file.exists(rsem_file)) {
  stop("File not found: ", rsem_file, "\n",
       "Run from the repository root; working directory: ", getwd())
}

cat(sprintf("  File size: %.1f MB\n", file.size(rsem_file) / 1e6))
cat("  Reading...\n")

t_start <- Sys.time()
expr_all <- fread(rsem_file, header = TRUE)
t_elapsed <- difftime(Sys.time(), t_start, units = "secs")

gene_col <- names(expr_all)[1]
all_samples <- setdiff(names(expr_all), gene_col)

cat(sprintf("  Done (%.0f s)\n", as.numeric(t_elapsed)))
cat(sprintf("  All data: %d genes × %d samples\n", nrow(expr_all), length(all_samples)))
cat(sprintf("  Gene column: '%s'\n", gene_col))
cat(sprintf("  Sample IDs, e.g.: %s\n", paste(head(all_samples, 3), collapse = ", ")))

# ── 3. Matching sample IDs ───────────────────────────────────────────────────
cat("\n", strrep("=", 70), "\n")
cat(">>> STEP 3: matching sample IDs\n")
cat(strrep("=", 70), "\n\n")

# Xena files may use IDs with '-' replaced by '.'; the format is detected automatically

match_samples <- function(target_ids, available_ids, label) {
  # 1: direct match
  m1 <- intersect(target_ids, available_ids)
  if (length(m1) >= length(target_ids) * 0.5) {
    cat(sprintf("  [%s] Direct match: %d / %d\n", label, length(m1), length(target_ids)))
    return(list(matched = m1, conversion = "none"))
  }

  # 2: '-' -> '.'
  target_dots <- gsub("-", ".", target_ids)
  m2 <- intersect(target_dots, available_ids)
  if (length(m2) > length(m1)) {
    cat(sprintf("  [%s] Dot-converted match: %d / %d\n", label, length(m2), length(target_ids)))
    return(list(matched = m2, conversion = "dash_to_dot",
                map = setNames(target_ids, target_dots[target_dots %in% available_ids])))
  }

  # 3: '.' -> '-'
  target_dashes <- gsub("\\.", "-", target_ids)
  m3 <- intersect(target_dashes, available_ids)
  if (length(m3) > length(m1)) {
    cat(sprintf("  [%s] Dash-converted match: %d / %d\n", label, length(m3), length(target_ids)))
    return(list(matched = m3, conversion = "dot_to_dash",
                map = setNames(target_ids, target_dashes[target_dashes %in% available_ids])))
  }

  cat(sprintf("  [%s] Best match: %d / %d (direct)\n", label, length(m1), length(target_ids)))
  return(list(matched = m1, conversion = "none"))
}

tumor_match  <- match_samples(esca_tumor_ids,  all_samples, "ESCA Tumor")
normal_match <- match_samples(esca_normal_ids,  all_samples, "ESCA Normal")
gtex_match   <- match_samples(gtex_eso_ids,     all_samples, "GTEx Eso")

# ── 4. Extract and save ──────────────────────────────────────────────────────
cat("\n", strrep("=", 70), "\n")
cat(">>> STEP 4: extracting and saving the expression data\n")
cat(strrep("=", 70), "\n\n")

# ─── 4-1. ESCA tumor ─────────────────────────────────────────────────────────
if (length(tumor_match$matched) > 0) {
  cols_tumor <- c(gene_col, tumor_match$matched)
  expr_tumor <- expr_all[, ..cols_tumor]

  # back to the TCGA IDs if a dash-to-dot conversion was needed
  if (tumor_match$conversion == "dash_to_dot") {
    old_names <- names(expr_tumor)[-1]
    new_names <- gsub("\\.", "-", old_names)
    setnames(expr_tumor, old_names, new_names)
  }

  out_file_1 <- file.path(DATA, "expr_ESCA_tumor.csv")
  fwrite(expr_tumor, out_file_1)
  cat(sprintf("  [1] %s\n", out_file_1))
  cat(sprintf("      %d genes × %d samples (%.1f MB)\n",
              nrow(expr_tumor), ncol(expr_tumor) - 1,
              file.size(out_file_1) / 1e6))
} else {
  cat("  ESCA Tumor: no matching samples, file not written\n")
}

# ─── 4-2. ESCA solid tissue normal ───────────────────────────────────────────
if (length(normal_match$matched) > 0) {
  cols_normal <- c(gene_col, normal_match$matched)
  expr_normal <- expr_all[, ..cols_normal]

  if (normal_match$conversion == "dash_to_dot") {
    old_names <- names(expr_normal)[-1]
    new_names <- gsub("\\.", "-", old_names)
    setnames(expr_normal, old_names, new_names)
  }

  out_file_2 <- file.path(DATA, "expr_ESCA_normal_TCGA.csv")
  fwrite(expr_normal, out_file_2)
  cat(sprintf("  [2] %s\n", out_file_2))
  cat(sprintf("      %d genes × %d samples (%.1f MB)\n",
              nrow(expr_normal), ncol(expr_normal) - 1,
              file.size(out_file_2) / 1e6))
} else {
  cat("  ESCA Normal: no matching samples, file not written\n")
}

# ─── 4-3. GTEx esophagus ─────────────────────────────────────────────────────
if (length(gtex_match$matched) > 0) {
  cols_gtex <- c(gene_col, gtex_match$matched)
  expr_gtex <- expr_all[, ..cols_gtex]

  out_file_3 <- file.path(DATA, "expr_GTEx_esophagus.csv")
  fwrite(expr_gtex, out_file_3)
  cat(sprintf("  [3] %s\n", out_file_3))
  cat(sprintf("      %d genes × %d samples (%.1f MB)\n",
              nrow(expr_gtex), ncol(expr_gtex) - 1,
              file.size(out_file_3) / 1e6))

  # ─── GTEx subtype map (mucosa / muscularis / gastroesophageal junction) ───
  gtex_subtype <- gtex_eso[sample %in% gtex_eso_ids,
                           .(sample, subtype = detailed_category,
                             gender = `_gender`)]

  # the sample ID as written in the expression file
  if (gtex_match$conversion == "dash_to_dot") {
    gtex_subtype[, sample_expr := gsub("-", ".", sample)]
  } else {
    gtex_subtype[, sample_expr := sample]
  }

  out_file_4 <- file.path(DATA, "GTEx_esophagus_subtype_map.csv")
  fwrite(gtex_subtype, out_file_4)
  cat(sprintf("  [4] %s (subtype map)\n", out_file_4))
  cat("      Subtypes:\n")
  print(gtex_subtype[, .N, by = subtype][order(-N)])

} else {
  cat("  GTEx Esophagus: no matching samples, file not written\n")
}

# ── 5. Summary ───────────────────────────────────────────────────────────────
cat("\n", strrep("=", 70), "\n")
cat(">>> SUMMARY\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("  [1] expr_ESCA_tumor.csv              %3d samples (primary tumor)\n", length(tumor_match$matched)))
cat(sprintf("  [2] expr_ESCA_normal_TCGA.csv        %3d samples (solid tissue normal)\n", length(normal_match$matched)))
cat(sprintf("  [3] expr_GTEx_esophagus.csv          %3d samples (GTEx esophagus)\n", length(gtex_match$matched)))
cat("  [4] GTEx_esophagus_subtype_map.csv   (subtype map)\n")
cat("\n  Values: log2(norm_count + 1) as distributed by Xena; no further transformation.\n")

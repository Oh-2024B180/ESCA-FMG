# Mutation-guided transcriptomic stratification of esophageal cancer — analysis code

Analysis code for

> Oh J, Byun J, Lee H. *Mutation-Guided Transcriptomic Stratification of Esophageal Cancer: Molecular Subtypes and
> a Candidate Prognostic Signature.* Submitted to *Current Issues in Molecular Biology* (MDPI).

The repository contains code only. No data are included (see **Data**).

## Contents

| Folder | What it does |
|---|---|
| `scripts/` | The analysis pipeline, stages 00–14, in the order of the Methods |
| `scripts/helpers/` | Steps that belong to a stage but were run as separate scripts |
| `scripts/figures/` | Figures 2–11 and Figure S1 of the article, drawn at their printed size from the stage outputs |
| `sessionInfo/` | The R session information that stages 01–14 wrote when the analysis was run (R version, attached and loaded packages with their versions) |

Each script states its inputs and outputs in its header. Random seed: 1262118388.

### Scripts and the parts of the article they produce

| Script | Article |
|---|---|
| `00_prepare_expression.R` | Expression matrices of the TCGA-ESCA primary tumors (181) and the GTEx esophagus samples, taken from the UCSC Xena TCGA TARGET GTEx (Toil) data (Methods 2.3); writes `data/expr_ESCA_tumor.csv`, `data/expr_GTEx_esophagus.csv` and `data/GTEx_esophagus_subtype_map.csv` |
| `01_02_FMG_analysis.R` | Variant filtering and the frequently mutated gene (FMG) panel of the CODA cohort (Methods 2.2 and 2.4; Table 2, CODA columns). Needs the restricted CODA data; without it, `tables/Tab1_2_FMG_list_primary.csv` (column `gene`) can be written from the 18 genes of Table 2 |
| `helpers/01_02_hotspot_check.py` | What the variant filters did to TP53 and to known hotspot variants (Discussion, second paragraph); prints only |
| `03_FMG_GO_KEGG_enrichment.R` | GO and KEGG enrichment of the FMG (Results 3.1; Figure 2) |
| `04_TCGA_mutation_oncoprint.R` | FMG mutations in TCGA-ESCA (Results 3.1; Figure 3a; Table 2, TCGA columns, from `tables/Tab4_TCGA_mutation_frequency.csv`) |
| `04c_FMG_mutation_vs_expression.R` | Expression of each FMG in tumors with and without its mutation (Results 3.1, last paragraph; Table S1) |
| `05_TCGA_CNV_analysis.R` | Copy-number gain and loss (Results 3.1; Figure 3b) |
| `06_expression_tumor_vs_normal.R` | TCGA-ESCA tumors versus GTEx mucosa (Results 3.1; Figure 3c) |
| `07_consensus_clustering.R` | Consensus clustering on FMG expression (Results 3.2; Figure 4) |
| `helpers/07_make_Fig3F_BH_FMG_boxplot.R` | Wilcoxon tests of FMG expression between the clusters (`tables/Tab7C_*`; Figure 4g, whose script checks its p values against this table) |
| `08_subtype_immune_stromal.R`, `helpers/08_build_bindea28.R`, `helpers/08_rerun_cibersort_QNfalse.R`, `helpers/08_import_cibersortx_result.R` | ESTIMATE, CIBERSORTx, ssGSEA, immune checkpoint and HLA comparisons between clusters (Results 3.3; Figure 5) |
| `09_DEG_between_clusters.R`, `helpers/09_split_enrichment.R` | Differential expression between clusters and its enrichment (Results 3.4; Figure 6; Table S2, from `tables/Tab9_DEG_significant_C1_vs_C2.csv`) |
| `10_LASSO_Cox_model.R` | Univariate Cox screening, LASSO-Cox risk score and risk groups (Results 3.5; Figure 7; Table 3, univariate and LASSO columns; Table 4, hazard ratio and log-rank p of the training cohort) |
| `10b_bootstrap_optimism.R` | Optimism-corrected C-index and selection stability over 1,000 bootstrap resamples (Results 3.5; Table 4, optimism-corrected row) |
| `helpers/11_fix_cindex_timeROC_tables.R` | Discrimination tables (C-index and time-dependent AUC with 95% CI) of the training and external cohorts |
| `11_external_validation_GSE53625.R` | GSE53625 processing and probe annotation (Methods 2.9) |
| `11b_external_on_training_scale.R`, `11c_external_figures.R` | External evaluation with the score on the training scale (Results 3.6; Figure 8; Table 4, external row and the apparent C-index with 95% CI of the training cohort) |
| `12_RS_clinical_immune_association.R` | Risk score versus clinical variables and 102 immune features (Results 3.7; Figure 9; Table S3, from `tables/Tab12_RS_immune_spearman.csv` and `tables/Tab12_RS_immune_wilcoxon.csv`) |
| `12_RS_histology_alignment.R` | Risk score by histological type (Results 3.7; Table 3, median expression by histology and "Weighted toward") |
| `12b_cibersortx_confidence_sensitivity.R` | Risk score – CIBERSORTx correlations without the non-confident sample (Methods 2.6, Results 3.7) |
| `12c_within_histology_discrimination.R` | Discrimination of the risk score within ESCC and within EAC (Methods 2.10 and 2.14, Results 3.7; Table 4, rows "Within histology") |
| `13_nomogram_construction.R`, `helpers/13_compute_nomogram_timeROC.R`, `helpers/13_fix_Fig9_nomogram_calibration.R` | Multivariable Cox model, nomogram, C-index comparison, calibration and time-dependent AUC (Results 3.8; Figure 10; Table 4, nomogram, risk-score-alone and clinical-variables-alone rows; Table 5) |
| `14_treatment_response_prediction.R` | Immunophenoscore comparison (Results 3.9; Figure 11b) and the first TIDE and oncoPredict runs (Methods 2.12–2.13) |
| `helpers/14a_hgnc_symbol_audit.py`, `14a_tide_rerun.py`, `14a_tide_rerun.R` | HGNC symbol audit of the expression genes, TIDE with TIDEpy and the risk-group comparison reported in the article (Results 3.9; Figure 11a) |
| `14b_oncopredict_predictions.R`, `14c_oncopredict_vignette.R` | oncoPredict predictions following the package vignette (`14b`, written to `tables/Tab14Cb_oncoPredict_IC50_vignette.csv`) and the risk-group comparison (`14c`; Results 3.10; Figure 11c; Table S4, from `tables/Tab14C_oncoPredict_drug_test.csv`) |
| `figures/*.R` | Figures 2–11 as printed (16.5 cm wide) and Figure S1 (`redraw_fig2_fig6.R`), from the tables written by the stages above |

Figure 1 (study design) and Table 1 (datasets) were drawn and compiled by hand and have no script.

Note on stage 08: stage 08 and `helpers/08_rerun_cibersort_QNfalse.R` first ran CIBERSORT locally (`data/CIBERSORT.R`,
`data/LM22.txt`); both runs failed the deconvolution check (mean correlation below 0.03). The second run also writes the
linear-scale mixture (`cache/tcga_mixture_linear_for_cibersort.tsv`) that was submitted to the CIBERSORTx web tool
(B-mode). `helpers/08_import_cibersortx_result.R` reads the CIBERSORTx result and overwrites the CIBERSORT tables of
stage 08 (`tables/Tab8B_*`) that the later stages use. The order therefore matters: stage 08,
`helpers/08_rerun_cibersort_QNfalse.R`, the CIBERSORTx web job on the mixture file it wrote, and then
`helpers/08_import_cibersortx_result.R`; run the other way round, the failed local fractions would reach stages 12
and 12b.

Note on stage 13 and Table 5: `tables/Tab13_univariate_Cox_clinical.csv` stores sex and histology as male vs. female and
ESCC vs. EAC; Table 5 gives the reverse contrasts (female vs. male, EAC vs. ESCC), whose hazard ratios and confidence
limits are the reciprocals of the stored ones. The multivariable table already uses the Table 5 contrasts.

The histological types of the CODA cohort given in Methods 2.2 are not computed by any script here.

Notes on `scripts/figures/`:

* `rerun_ccp_fig4.R` and `rerun_lasso_fig7.R` run the consensus-clustering and LASSO steps of stages 07 and 10 again with
  the same inputs and seed to recover values for Figure 4a–c and Figure 7a–b, and stop unless the result equals the
  saved stage output.
* `v3_fig10C_IPS.R` computes the per-sample immunophenoscore (IOBR `IPS_calculation`, the values tested in stage 14) and
  caches it in `cache/v3_fig10C_ips_long.rds` for Figure 11b.
* Figure 10b keeps the data area of the calibration plot drawn by `helpers/13_fix_Fig9_nomogram_calibration.R` (bootstrap
  calibration) and draws new axes; the image of that area is not part of this repository.

## Data (not included)

| Data | Files in `data/` | Source | Used by |
|---|---|---|---|
| CODA cohort: tumor-only somatic variant calls of 98 Korean patients | `variants_eso_final_masked.csv` | Clinical & Omics Data Archive, Korea Disease Control and Prevention Agency (https://coda.nih.go.kr); restricted, available subject to approval (approval number in the article). CODA distributes one VCF per patient; within the CODA analysis environment the variants with FILTER = PASS and a non-reference genotype were extracted into this table (Methods 2.2), of which `01_02` requires the columns `sample_id_masked`, `chrom`, `pos`, `ref` and `alt` | `01_02` |
| Optional list of variants to keep as recurrent sites | `recurrent_site_rescue_list.csv` | `01_02` writes an empty template (`hgvs_id`, `rescue_reason`) when the file is absent | `01_02` |
| TCGA-ESCA MC3 somatic mutations, GISTIC2 thresholded copy number, clinical and survival data | `ESCA_mc3.csv`, `Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes`, `TCGA.ESCA.sampleMap_ESCA_clinicalMatrix.csv`, `survival_ESCA_survival.csv` | UCSC Xena TCGA hub (https://xenabrowser.net), datasets `mc3/ESCA_mc3.txt`, `TCGA.ESCA.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes`, `TCGA.ESCA.sampleMap/ESCA_clinicalMatrix` and `survival/ESCA_survival.txt`; the mutation, clinical and survival tables were saved as comma-separated files under the names given here | `00`, stages 04–14 |
| UCSC Xena TCGA TARGET GTEx expression (log2(norm_count + 1)) and its phenotype table | `TcgaTargetGtex_RSEM_Hugo_norm_count`, `TcgaTargetGTEX_phenotype.csv` | UCSC Xena Toil hub (https://xenabrowser.net), datasets `TcgaTargetGtex_RSEM_Hugo_norm_count` and `TcgaTargetGTEX_phenotype` (the phenotype table saved as a comma-separated file) | `00` |
| TCGA-ESCA tumor expression (181 tumors) and GTEx esophagus expression | `expr_ESCA_tumor.csv`, `expr_GTEx_esophagus.csv`, `GTEx_esophagus_subtype_map.csv` | written by `00_prepare_expression.R` from the two rows above | stages 04c and 06–14 |
| GSE53625 (Agilent lncRNA + mRNA microarray, GPL18109; overall survival) | `GSE53625_series_matrix.txt.gz` | NCBI Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo) | `11`, `11b`, `11c` |
| CIBERSORT R script and LM22 signature matrix | `CIBERSORT.R`, `LM22.txt` | CIBERSORT website (https://cibersort.stanford.edu; registration required) | stage `08` (needed for it to run) and `helpers/08_rerun_cibersort_QNfalse.R` |
| CIBERSORTx results (LM22, B-mode) | `cibersortx_LM22_TCGA_ESCA_Bmode.csv` | CIBERSORTx web tool (https://cibersortx.stanford.edu) | `helpers/08_import_cibersortx_result.R`, `12b` |
| GDSC2 expression and drug-response reference data | `GDSC2_Expr (RMA Normalized and Log Transformed).rds`, `GDSC2_Res.rds` | prepared by the oncoPredict developers and distributed through their Open Science Framework project (https://osf.io/c6tfx/), as given in the oncoPredict README; not part of the R package itself | `14`, `14b` |
| HGNC complete set (approved symbols, aliases, previous symbols, locus types) | `hgnc_complete_set.txt` | HUGO Gene Nomenclature Committee (https://www.genenames.org); the article used the file dated 2026-04-19, and a later release can classify a few symbols differently | `helpers/14a_hgnc_symbol_audit.py` |
| An earlier TIDEpy run (April 2026): its input and output | `esca_tide_output_hgnc_self.txt`, `tide_2026-04/` | not distributed; `helpers/14a_hgnc_symbol_audit.py` writes the gene-symbol audit (`tide_2026-04/gene_symbol_audit.csv`) that `14a_tide_rerun.py` needs | `14` (first TIDE comparison, skipped without the file), `14a_tide_rerun.py` (comparisons with the April run, skipped without the files) |

Written by the scripts and read by later ones: `Bindea_28_signatures.tsv` (`helpers/08_build_bindea28.R`, from the IOBR
signature collection and Table S1 of Bindea et al. 2013), `esca_tide_input_hgnc_centered.txt` and
`esca_tide_output_hgnc_centered.txt` (`14a_tide_rerun.py`).

## Running

Run every script from the repository root. Inputs are read from `data/`; outputs go to `tables/`, `figures/`, `logs/` and
`cache/`, which are ignored by git. Create the folders first:

```bash
mkdir -p data tables figures/final figures/remake logs/sessionInfo cache
```

Run `00_prepare_expression.R` first; it reads the 8 GB Xena expression matrix once and writes the smaller expression
files the later stages read. Then run the stages in numerical order and each helper after the stage whose number it
carries, except
`helpers/08_build_bindea28.R`, which writes `data/Bindea_28_signatures.tsv` for stage 08 and runs before it. Stage 14
runs in this order: `14_treatment_response_prediction.R`, `helpers/14a_hgnc_symbol_audit.py`, `14a_tide_rerun.py`,
`14a_tide_rerun.R`, `14b_oncopredict_predictions.R`, `14c_oncopredict_vignette.R`. Stage 14 writes first versions of
`tables/Tab14A_*` and `tables/Tab14C_*`; `14a_tide_rerun.R` and `14c` overwrite them with the values reported in the
article, and `14b` compares its first run with the stage 14 predictions, so it runs before `14c`. The figure scripts run
last: `rerun_ccp_fig4.R`, `rerun_lasso_fig7.R` and `v3_fig10C_IPS.R` before the `redraw_*` scripts, which read their
results.

The scripts have not been re-run end to end in this repository layout.

## Software

R: stages 01–11 ran in R 4.4.3 and stages 12–14 in R 4.4.1. The session information each stage wrote at the time
(R version, attached and loaded packages with their versions) is in `sessionInfo/`, one file per stage. The scripts
added after those runs (04c, 10b, 11b, 11c, 12b, 12c, 14a, 14b, 14c and the figure scripts) ran in R 4.4.1;
`00_prepare_expression.R` needs only data.table.

Python 3.14.4 for `14a_tide_rerun.py`, with pandas 3.0.2, numpy 2.4.4 and TIDEpy 1.3.9; the two helper scripts in
Python need pandas only.

Package versions as stated in the article: TxDb.Hsapiens.UCSC.hg19.knownGene 3.2.2, org.Hs.eg.db 3.20.0,
clusterProfiler 4.14.6, ComplexHeatmap 2.22.0, ConsensusClusterPlus 1.70.0, estimate 1.0.13, GSVA 2.0.7, IOBR 2.2.0
(signature collection) and 0.99.0 (immunophenoscore), limma 3.62.2, biomaRt 2.62.1, glmnet 4.1-10, survival 3.8-6,
timeROC 0.4.1 (0.4 for the nomogram), AnnoProbe 0.1.7, rms 8.1-1, oncoPredict 1.2, TIDEpy 1.3.9.

Packages from outside CRAN and Bioconductor: estimate from R-Forge
(`install.packages("estimate", repos = "http://r-forge.r-project.org", dependencies = TRUE)`, as given by its authors),
IOBR from GitHub (https://github.com/IOBR/IOBR), TIDEpy from GitHub (https://github.com/jingxinfu/TIDEpy), and the
CIBERSORT R script from the CIBERSORT website (see **Data**).

Figure scripts: R 4.4.1 with ggplot2 4.0.2, patchwork 1.3.2, ggrepel 0.9.8, ComplexHeatmap 2.20.0, ConsensusClusterPlus
1.68.0 and glmnet 4.1.10.

## License

MIT (see `LICENSE`).

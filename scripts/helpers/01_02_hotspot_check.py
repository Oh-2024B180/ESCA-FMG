###############################################################################
## Stage 1/2 helper: what the CODA variant filters did to TP53 and to known hotspots (Discussion, second paragraph)
##
## Prints which TP53 / PIK3CA hotspot variants are in the stage 01_02 QC table and whether any filter removed them,
## the variants removed as possible germline (gnomAD AF 0.001-0.01), and the raw calls (SNVs and indels) in the
## TP53 locus. Writes nothing.
## Inputs  tables/variant_qc_flags_v10_1.csv, tables/
##         mutation_matrix_primary_v10_1.csv (stage 01_02 outputs), data/variants_eso_final_masked.csv
##         (restricted CODA data)
###############################################################################
import io, sys
import pandas as pd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_rows", 200)
q = pd.read_csv("tables/variant_qc_flags_v10_1.csv", low_memory=False)
print("variants:", len(q))
for c in ("flag_common_germline", "flag_possible_germline", "flag_consequence_missing", "flag_functional_primary", "flag_primary_drop_variant"):
    print("  %-28s %s" % (c, q[c].value_counts(dropna=False).to_dict()))
print("\nexclusion_reason_primary:", q.exclusion_reason_primary.value_counts(dropna=False).head(12).to_dict())

tp = q[q.genename == "TP53"].copy()
cols = ["hgvs_id", "picked_consequence", "all_consequence", "gnomad_af_max", "n_samples", "median_tumor_af",
        "flag_common_germline", "flag_possible_germline", "flag_consequence_missing", "flag_functional_primary", "exclusion_reason_primary"]
print("\nTP53 variants in the QC table:", len(tp), "| samples carrying them:", int(tp.n_samples.sum()))
print(tp[cols].sort_values("n_samples", ascending=False).to_string())

# well-known hotspots (hg19)
HOT = {("chr17", 7578406): "TP53 R175H", ("chr17", 7577538): "TP53 R248Q", ("chr17", 7577539): "TP53 R248W",
       ("chr17", 7577120): "TP53 R273H", ("chr17", 7577121): "TP53 R273C", ("chr17", 7577548): "TP53 G245S",
       ("chr17", 7577094): "TP53 R282W", ("chr17", 7578190): "TP53 Y220C", ("chr3", 178936091): "PIK3CA E545K",
       ("chr3", 178936082): "PIK3CA E542K", ("chr3", 178952085): "PIK3CA H1047R"}
q["key"] = list(zip(q.chrom, q.pos))
hot = q[q.key.isin(HOT.keys())].copy(); hot["hotspot"] = hot.key.map(HOT)
print("\nknown hotspots present:\n", hot[["hotspot", "hgvs_id", "n_samples", "gnomad_af_max", "flag_possible_germline", "flag_common_germline", "exclusion_reason_primary"]].to_string())

band = q[(q.flag_possible_germline == True)]
print("\nvariants removed as possible germline (0.001-0.01):", len(band), "| in the 18 FMG:",
      band[band.genename.isin(["TP53", "ATM", "FANCD2", "LRP1B", "NOTCH1", "RUNX1", "NOTCH3", "PKHD1", "BRCA2", "NTRK1", "PIK3CA", "ASXL1", "ERCC6", "FAT3", "FGFR2", "KMT2D", "ROS1", "SLX4"])].genename.value_counts().to_dict())
print("possible-germline band in common cancer genes:",
      band[band.genename.isin(["TP53", "PIK3CA", "NOTCH1", "NFE2L2", "CDKN2A", "KMT2D", "FAT1", "ZNF750", "ARID1A", "SMAD4", "ERBB2", "KRAS", "PTEN"])][["genename", "hgvs_id", "n_samples", "gnomad_af_max"]].to_string())

# raw variant table: TP53 indels and all TP53 records before filtering
v = pd.read_csv("data/variants_eso_final_masked.csv", usecols=["sample_id_masked", "chrom", "pos", "ref", "alt", "AF", "filter"], low_memory=False)
t = v[(v.chrom == "chr17") & (v.pos.between(7565097, 7590856))]
t = t.assign(kind=["indel" if len(r) != len(a) else "snv" for r, a in zip(t.ref, t.alt)])
print("\nraw calls in the TP53 locus (hg19 chr17:7,565,097-7,590,856):", len(t), "| samples:", t.sample_id_masked.nunique(),
      "| by kind:", t.kind.value_counts().to_dict(), "| samples with an indel:", t[t.kind == "indel"].sample_id_masked.nunique())
coding = t[t.pos.between(7572927, 7579912)]
print("   within the coding span (7,572,927-7,579,912):", len(coding), "calls in", coding.sample_id_masked.nunique(), "samples; indels in",
      coding[coding.kind == "indel"].sample_id_masked.nunique(), "samples; filter values:", coding["filter"].value_counts().head(5).to_dict())
mm = pd.read_csv("tables/mutation_matrix_primary_v10_1.csv")
print("mutation matrix columns:", list(mm.columns)[:5], "| shape", mm.shape)

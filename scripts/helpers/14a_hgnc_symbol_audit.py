###############################################################################
## Stage 14A helper: HGNC symbol audit of the expression genes (the "keep" list used by scripts/14a_tide_rerun.py)
##
## Inputs  data/expr_ESCA_tumor.csv (first column = gene symbol)
##         data/hgnc_complete_set.txt (HGNC complete set)
##
## Outputs data/tide_2026-04/gene_list_from_expr.csv, gene_list_summary.txt,
##         data/tide_2026-04/gene_symbol_audit.csv (read by 14a_tide_rerun.py), gene_symbol_keep_list.csv,
##         gene_symbol_drop_list.csv, gene_symbol_audit_summary.txt
###############################################################################
from pathlib import Path
import pandas as pd

expr_file = Path("data/expr_ESCA_tumor.csv")
hgnc_file = Path("data/hgnc_complete_set.txt")
OUT = Path("data/tide_2026-04")
OUT.mkdir(parents=True, exist_ok=True)
for f in (expr_file, hgnc_file):
    if not f.exists():
        raise FileNotFoundError(f"input file not found: {f}")

# ---------------------------------------------------------------- 1. gene list of the expression file
df = pd.read_csv(expr_file)
if df.shape[1] < 2:
    raise ValueError("expr_ESCA_tumor.csv needs a gene column and sample columns")
first_col = df.columns[0]
gene_clean = df[first_col].copy().astype(str).str.strip()
gene_clean = gene_clean[gene_clean.notna() & (gene_clean != "") & (gene_clean.str.lower() != "nan")]
n_total_rows, n_nonempty, n_unique = len(df), len(gene_clean), gene_clean.nunique(dropna=True)
gene_unique = pd.Series(gene_clean.unique(), name="gene").to_frame()
gene_unique["gene_length"] = gene_unique["gene"].str.len()
gene_unique["has_space"] = gene_unique["gene"].str.contains(r"\s", regex=True)
gene_unique["starts_with_LOC"] = gene_unique["gene"].str.upper().str.startswith("LOC")
gene_unique["starts_with_LINC"] = gene_unique["gene"].str.upper().str.startswith("LINC")
gene_unique["starts_with_MIR"] = gene_unique["gene"].str.upper().str.startswith("MIR")
gene_unique["contains_dash"] = gene_unique["gene"].str.contains("-", regex=False)
gene_unique.to_csv(OUT / "gene_list_from_expr.csv", index=False)
(OUT / "gene_list_summary.txt").write_text("\n".join([
    "gene list extraction summary", f"input file: {expr_file}", f"gene column used: {first_col}",
    f"total rows in expression file: {n_total_rows}", f"non-empty gene entries: {n_nonempty}",
    f"unique genes: {n_unique}", f"duplicate gene entries removed: {n_nonempty - n_unique}"]), encoding="utf-8")
print("genes: %d rows, %d unique" % (n_total_rows, n_unique))

# ---------------------------------------------------------------- 2. HGNC audit
gene_df = gene_unique
hgnc = pd.read_csv(hgnc_file, sep="\t", dtype=str, low_memory=False)
for col in ("symbol", "status"):
    if col not in hgnc.columns:
        raise ValueError(f"HGNC file lacks the column {col}")
alias_col = "alias_symbol" if "alias_symbol" in hgnc.columns else None
prev_col = "prev_symbol" if "prev_symbol" in hgnc.columns else None
locus_type_col = "locus_type" if "locus_type" in hgnc.columns else None
locus_group_col = "locus_group" if "locus_group" in hgnc.columns else None
hgnc = hgnc.fillna("")
approved_hgnc = hgnc[hgnc["status"].str.upper().eq("APPROVED")].copy()


def norm(x):
    return str(x).strip().upper()


approved_hgnc["symbol_norm"] = approved_hgnc["symbol"].map(norm)
symbol_map = {}
for _, row in approved_hgnc.iterrows():
    symbol_map.setdefault(row["symbol_norm"], []).append(row)


def split_pipe_values(x):
    if pd.isna(x):
        return []
    x = str(x).strip()
    if x == "":
        return []
    return [v.strip() for v in x.split("|") if v.strip()]


alias_map, prev_map = {}, {}
for _, row in approved_hgnc.iterrows():
    approved_symbol = row["symbol"]
    if alias_col:
        for alias in split_pipe_values(row[alias_col]):
            alias_map.setdefault(norm(alias), set()).add(approved_symbol)
    if prev_col:
        for prev in split_pipe_values(row[prev_col]):
            prev_map.setdefault(norm(prev), set()).add(approved_symbol)
all_symbol_status = {}
for _, row in hgnc.iterrows():
    all_symbol_status.setdefault(norm(row["symbol"]), set()).add(row["status"])

audit_rows = []
for gene in gene_df["gene"].astype(str):
    g = gene.strip()
    g_norm = norm(g)
    match_type, approved_symbol, status, locus_type, locus_group, action, note = "unmapped", "", "", "", "", "drop", ""
    if g_norm in symbol_map and len(symbol_map[g_norm]) == 1:            # exact approved symbol
        row = symbol_map[g_norm][0]
        approved_symbol, status = row["symbol"], row["status"]
        locus_type = row[locus_type_col] if locus_type_col else ""
        locus_group = row[locus_group_col] if locus_group_col else ""
        match_type = "exact_approved"
    elif g_norm in all_symbol_status:                                     # exact non-approved symbol
        status = "|".join(sorted(all_symbol_status[g_norm]))
        match_type = "exact_nonapproved"
        note = "symbol exists in HGNC but not as current approved symbol"
    elif g_norm in alias_map:                                             # alias rescue
        targets = sorted(alias_map[g_norm])
        if len(targets) == 1:
            approved_symbol = targets[0]
            row = approved_hgnc[approved_hgnc["symbol"] == approved_symbol].iloc[0]
            status = row["status"]
            locus_type = row[locus_type_col] if locus_type_col else ""
            locus_group = row[locus_group_col] if locus_group_col else ""
            match_type = "alias_rescued"
        else:
            match_type, approved_symbol = "alias_ambiguous", "|".join(targets)
            note = "alias maps to multiple approved symbols"
    elif g_norm in prev_map:                                              # previous symbol rescue
        targets = sorted(prev_map[g_norm])
        if len(targets) == 1:
            approved_symbol = targets[0]
            row = approved_hgnc[approved_hgnc["symbol"] == approved_symbol].iloc[0]
            status = row["status"]
            locus_type = row[locus_type_col] if locus_type_col else ""
            locus_group = row[locus_group_col] if locus_group_col else ""
            match_type = "previous_rescued"
        else:
            match_type, approved_symbol = "previous_ambiguous", "|".join(targets)
            note = "previous symbol maps to multiple approved symbols"
    # conservative TIDE input rule: keep only approved symbols of protein-coding genes
    if match_type in ["exact_approved", "alias_rescued", "previous_rescued"]:
        if str(status).upper() != "APPROVED":
            action = "drop"
            note = f"{note}; non-approved status".strip("; ")
        elif locus_type and locus_type != "gene with protein product":
            action = "drop_noncoding"
        else:
            action = "keep"
    elif match_type in ["alias_ambiguous", "previous_ambiguous"]:
        action = "drop_ambiguous"
    elif match_type == "exact_nonapproved":
        action = "drop_nonapproved"
    else:
        action = "drop_unmapped"
    audit_rows.append({"gene_original": g, "gene_norm": g_norm, "match_type": match_type,
                       "approved_symbol": approved_symbol, "status": status, "locus_type": locus_type,
                       "locus_group": locus_group, "action": action, "note": note})

audit_df = pd.DataFrame(audit_rows)
audit_df.to_csv(OUT / "gene_symbol_audit.csv", index=False)
keep_df = audit_df[audit_df["action"] == "keep"].copy()
drop_df = audit_df[audit_df["action"] != "keep"].copy()
keep_df.to_csv(OUT / "gene_symbol_keep_list.csv", index=False)
drop_df.to_csv(OUT / "gene_symbol_drop_list.csv", index=False)
lines = ["HGNC audit summary", f"total input genes: {len(audit_df)}", "", "match_type counts"]
lines += [f"{k}: {v}" for k, v in audit_df["match_type"].value_counts(dropna=False).items()]
lines += ["", "action counts"] + [f"{k}: {v}" for k, v in audit_df["action"].value_counts(dropna=False).items()]
(OUT / "gene_symbol_audit_summary.txt").write_text("\n".join(lines), encoding="utf-8")
print("audit: %d genes, keep %d" % (len(audit_df), len(keep_df)))

###############################################################################
## Stage 14A re-run, part 1 of 2 (2026-10-19): TIDE with a correctly transformed input
##
## Why: the April run (2026-04-19, according to the command history of the analysis computer) was
##   tidepy esca_tide_input_hgnc_selfsubtract.txt -o esca_tide_output_hgnc_self.txt -c Other --ignore_norm
## with TIDEpy 1.3.9, and its input had been made by make_tide_input_hgnc.py as log2(x + 1) of values that were
## already log2(norm_count + 1) (UCSC Xena), then centred on the gene-wise mean. The other April output, from the
## un-centred HGNC matrix without --ignore_norm, is identical, because TIDEpy itself applies log2(x + 1) and centring
## when more than 80% of the values are positive. Both April outputs therefore rest on a doubly logged matrix.
##
## What this script does:
##   1. rebuilds the April HGNC matrix (log2(norm_count + 1), one row per approved symbol, the highest-variance row
##      kept for duplicated symbols) from data/expr_ESCA_tumor.csv and the April symbol audit, and checks that it is
##      identical to the April file;
##   2. checks that TIDEpy here reproduces the April output from the April input exactly (same options);
##   3. centres the HGNC matrix on the gene-wise mean across the 181 tumours WITHOUT a second log, which is what TIDE
##      asks for ("normalized toward ... the all sample average"), and runs TIDEpy with the April options:
##      cancer = "Other", no previous immunotherapy, vthres = 0, ignore_norm (the input is already normalized).
##
## Inputs  data/expr_ESCA_tumor.csv; the April folder (read only):
##         data/tide_2026-04/{gene_symbol_audit.csv, esca_tide_input_hgnc.txt,
##         esca_tide_input_hgnc_selfsubtract.txt, esca_tide_output_hgnc_self.txt}
##         gene_symbol_audit.csv is required (written by scripts/helpers/14a_hgnc_symbol_audit.py); the other three
##         April files are used only for the checks in steps 1-3, which are skipped when the files are absent.
## Outputs data/esca_tide_input_hgnc_centered.txt, data/esca_tide_output_hgnc_centered.txt,
##         logs/14a_tide_rerun_py.log
## Part 2: scripts/14a_tide_rerun.R (risk-group tests as in stage 14A, Figure 11a panel).
###############################################################################
import io, os, sys
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
import tidepy
from tidepy.pred import TIDE

ROOT = ""  # repository root = working directory (see README.md)
APRIL = "data/tide_2026-04/"  # files of the April 2026 TIDE run (not distributed; see README.md)
OUT_IN = ROOT + "data/esca_tide_input_hgnc_centered.txt"
OUT_TIDE = ROOT + "data/esca_tide_output_hgnc_centered.txt"
LOG = ROOT + "logs/14a_tide_rerun_py.log"
OPTS = dict(cancer="Other", pretreat=False, vthres=0.0, ignore_norm=True)   # the April options


class Tee(io.TextIOBase):
    def __init__(self, *s): self.s = s
    def write(self, t):
        for x in self.s: x.write(t)
        return len(t)
    def flush(self):
        for x in self.s: x.flush()


log = open(LOG, "w", encoding="utf-8")
console = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.stdout = Tee(console, log)
print("=== Stage 14A TIDE re-run (part 1, Python %s, pandas %s, numpy %s, tidepy %s) ===" % (
    sys.version.split()[0], pd.__version__, np.__version__, tidepy.__version__))
print("TIDE options:", OPTS)

# ---------------------------------------------------------------- 1. the April HGNC matrix, rebuilt
expr = pd.read_csv(ROOT + "data/expr_ESCA_tumor.csv")
expr = expr.rename(columns={expr.columns[0]: "gene_original"})
audit = pd.read_csv(APRIL + "gene_symbol_audit.csv")
keep = audit[audit["action"] == "keep"].copy()
keep["approved_symbol"] = keep["approved_symbol"].fillna("").astype(str).str.strip()
keep["gene_original"] = keep["gene_original"].astype(str).str.strip()
keep["gene_final"] = np.where(keep["approved_symbol"] == "", keep["gene_original"], keep["approved_symbol"])
keep = keep[["gene_original", "gene_final"]].drop_duplicates()
expr["gene_original"] = expr["gene_original"].astype(str).str.strip()
m = pd.merge(expr, keep, on="gene_original", how="inner")
samples = [c for c in m.columns if c not in ("gene_original", "gene_final")]
m[samples] = m[samples].apply(pd.to_numeric, errors="coerce")
m = m.dropna(subset=samples, how="all").copy()
m["row_var"] = m[samples].var(axis=1, skipna=True)
m = m.sort_values(["gene_final", "row_var"], ascending=[True, False])
hgnc = m.drop_duplicates(subset="gene_final", keep="first").set_index("gene_final")[samples]
hgnc.index.name = "gene"
if os.path.exists(APRIL + "esca_tide_input_hgnc.txt"):   # April files are not distributed: the checks run only if present
    april = pd.read_csv(APRIL + "esca_tide_input_hgnc.txt", sep="\t", index_col=0)
    assert list(april.columns) == samples and april.index.equals(hgnc.index)
    d1 = float(np.nanmax(np.abs(april.to_numpy() - hgnc.to_numpy())))
    print("1. HGNC matrix rebuilt: %d genes x %d tumours; max |difference| from the April file = %.1e" % (
        hgnc.shape[0], hgnc.shape[1], d1))
    assert d1 < 1e-9
else:
    print("1. HGNC matrix rebuilt: %d genes x %d tumours (April file not present: comparison skipped)" % hgnc.shape)
print("   value range %.2f to %.2f (already log2(norm_count + 1)); %.1f%% of values > 0" % (
    np.nanmin(hgnc.to_numpy()), np.nanmax(hgnc.to_numpy()), 100 * (hgnc.to_numpy() > 0).mean()))


# ---------------------------------------------------------------- 2. April output reproduced?
def run(df):
    return TIDE(df.copy(), **OPTS)


def compare(a, b):
    a = a.loc[b.index]
    num = [c for c in b.columns if b[c].dtype != bool]
    dev = float((a[num].astype(float) - b[num].astype(float)).abs().max().max())
    same = all((a[c] == b[c]).all() for c in b.columns if b[c].dtype == bool)
    return dev, same


HAVE_OLD = all(os.path.exists(APRIL + f) for f in ("esca_tide_input_hgnc_selfsubtract.txt", "esca_tide_output_hgnc_self.txt"))
if HAVE_OLD:
    old_in = pd.read_csv(APRIL + "esca_tide_input_hgnc_selfsubtract.txt", sep="\t", index_col=0)
    old_out = pd.read_csv(APRIL + "esca_tide_output_hgnc_self.txt", sep="\t", index_col=0)
    again = run(old_in)
    dev, same = compare(again, old_out)
    print("2. April input re-run here: max |difference| from the April output = %.1e, flags identical: %s" % (dev, same))
    assert dev < 1e-9 and same
else:
    print("2. April input and output not present: check skipped")

# ---------------------------------------------------------------- 3. corrected input and new run
cen = hgnc.sub(hgnc.mean(axis=1, skipna=True), axis=0)
pos = 100 * (cen.to_numpy() > 0).mean()
print("3. corrected input: gene-wise centred log2(norm_count + 1), no second log; %.1f%% of values > 0 "
      "(TIDEpy would not re-normalize it either way)" % pos)
if HAVE_OLD:
    check_old = np.log2(hgnc + 1.0)
    check_old = check_old.sub(check_old.mean(axis=1), axis=0)
    assert float(np.nanmax(np.abs(check_old.to_numpy() - old_in.loc[hgnc.index, samples].to_numpy()))) < 1e-9
cen.to_csv(OUT_IN, sep="\t", float_format="%.10g")
new = run(cen)
new.to_csv(OUT_TIDE, sep="\t")
print("   wrote", OUT_IN, "and", OUT_TIDE, "(%d tumours x %d columns)" % new.shape)

# ---------------------------------------------------------------- 4. old vs new, per metric
if HAVE_OLD:                       # only with the April files (not distributed), like checks 2 and 3
    new = new.loc[old_out.index]
    print("\n4. April (double log) vs corrected input, per metric (Spearman rho across the 181 tumours; medians):")
    for c in ["TIDE", "Dysfunction", "Exclusion", "MDSC", "CAF", "TAM M2", "IFNG", "MSI Score", "CD274", "CD8", "CTL"]:
        rho = pd.Series(old_out[c]).rank().corr(pd.Series(new[c]).rank())
        print("   %-12s rho %.3f   median %8.4f -> %8.4f" % (c, rho, old_out[c].median(), new[c].median()))
    for c in ["Responder", "No benefits", "CTL.flag"]:
        print("   %-12s TRUE %3d -> %3d   (changed in %d tumours)" % (
            c, int(old_out[c].sum()), int(new[c].sum()), int((old_out[c] != new[c]).sum())))
else:
    print("4. April output not present: comparison skipped")
sys.stdout.flush()
sys.stdout = console
log.close()

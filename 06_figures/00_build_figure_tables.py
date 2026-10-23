"""Build the small per-line tables used by Figs 1A, 2A and 2B from scripted inputs.

Inputs
  samples.tsv            sample sheet (repo root)           sample_id, wgs_id, MYCN_amp, ecDNA_group, ...
  MYCN_CHEK1_FPKM.csv    05_differential_expression/03_export_FPKM_tables.R   Sample, ecDNA_group, MYCN, CHEK1
  MYCN_CN.tsv            02_ecDNA/03_extract_MYCN_CN.py     sample_id, MYCN_CN, source, feature, segment
  master_in.csv          existing master.csv: WGS / ecDNA / chromothripsis columns are carried over unchanged

Outputs (written to --out)
  MYCN_FPKM.txt          Sample / MYCN / MYCN_amp                 -> fig2B_MYCN_FPKM_boxplot.py
  FPKM_CN_MYCN.txt       Sample / GE / CN / MYCNamp               -> fig2A_MYCN_CN_vs_FPKM.py
  master.csv             master_in.csv with MYCN_CN, MYCN_FPKM, CHEK1_FPKM replaced -> fig1A_oncoprint.ipynb

Every expression and copy-number value in the three outputs therefore traces to a file written by a
pipeline script. Omitting --cn keeps the MYCN_CN column of master_in.csv (manually recorded, 2024) and
says so loudly; that is for comparison only and not for the published tables.

Usage
  python 00_build_figure_tables.py --samples samples.tsv --fpkm data/MYCN_CHEK1_FPKM.csv \
         --cn data/MYCN_CN.tsv --master-in data/master.csv --out data/
"""
import argparse, sys
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--samples", required=True)
ap.add_argument("--fpkm", required=True, help="MYCN_CHEK1_FPKM.csv from 03_export_FPKM_tables.R")
ap.add_argument("--cn", default=None, help="MYCN_CN.tsv from 02_ecDNA/03_extract_MYCN_CN.py")
ap.add_argument("--master-in", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

samples = pd.read_csv(a.samples, sep="\t")
fpkm = pd.read_csv(a.fpkm)
master = pd.read_csv(a.master_in)
assert len(samples) == 21 and len(fpkm) == 21, "expected 21 lines"
assert set(fpkm.Sample) == set(samples.sample_id), "sample IDs differ between samples.tsv and FPKM table"

# master.csv uses WGS-style names (SKN-BE2, SMS-KCNR, SJ-NB10, SKNFI ...): map through samples.tsv
id_map = dict(zip(samples.wgs_id, samples.sample_id))
missing = [s for s in master.Sample if s not in id_map]
assert not missing, f"master.csv samples not in samples.tsv wgs_id column: {missing}"
master["_sid"] = master.Sample.map(id_map)

# ── MYCN copy number ──────────────────────────────────────────────────────────
if a.cn:
    cn = pd.read_csv(a.cn, sep="\t")
    assert set(cn.sample_id) == set(samples.sample_id), "MYCN_CN.tsv does not cover the 21 lines"
    cn_map = dict(zip(cn.sample_id, cn.MYCN_CN))
    old = dict(zip(master._sid, master.MYCN_CN))
    changed = {s: (old[s], cn_map[s]) for s in cn_map if abs(float(old[s]) - float(cn_map[s])) > 0.05}
    print(f"MYCN_CN from {a.cn}: sources {cn.source.value_counts().to_dict()}", file=sys.stderr)
    if changed:
        print("MYCN_CN values that differ from master_in.csv (old -> new):", file=sys.stderr)
        for s, (o, n) in changed.items():
            print(f"  {s}: {o:g} -> {n:g}", file=sys.stderr)
    else:
        print("MYCN_CN identical to master_in.csv for all 21 lines", file=sys.stderr)
else:
    cn_map = dict(zip(master._sid, master.MYCN_CN))
    print("WARNING: no --cn given; MYCN_CN carried over from master_in.csv (manually recorded values). "
          "Not for the published tables.", file=sys.stderr)

# ── assemble ──────────────────────────────────────────────────────────────────
f = fpkm.set_index("Sample")
tbl = samples[["sample_id", "MYCN_amp", "ecDNA_group"]].copy()
tbl["MYCN"]  = tbl.sample_id.map(f.MYCN).round(2)
tbl["CHEK1"] = tbl.sample_id.map(f.CHEK1)
tbl["CN"]    = tbl.sample_id.map(cn_map)
tbl["amp"]   = tbl.MYCN_amp.eq("amp")
assert tbl.CN.notna().all(), "missing MYCN_CN for some lines"

TF = {True: "TRUE", False: "FALSE"}                      # same spelling as the original hand-made tables
fmt_cn = lambda v: f"{float(v):g}"                       # 273 not 273.0; 2.1 stays 2.1

# Fig 2B input
pd.DataFrame({"Sample": tbl.sample_id, "MYCN": tbl.MYCN, "MYCN_amp": tbl.amp.map(TF)}) \
  .to_csv(f"{a.out}/MYCN_FPKM.txt", sep="\t", index=False)

# Fig 2A input (amplified lines first, as in the original file)
t2 = tbl.sort_values("amp", ascending=False, kind="stable")
pd.DataFrame({"Sample": t2.sample_id, "GE": t2.MYCN, "CN": t2.CN.map(fmt_cn), "MYCNamp": t2.amp.map(TF)}) \
  .to_csv(f"{a.out}/FPKM_CN_MYCN.txt", sep="\t", index=False)

# Fig 1A input
master["MYCN_CN"]    = master._sid.map(cn_map)
master["MYCN_FPKM"]  = master._sid.map(dict(zip(tbl.sample_id, tbl.MYCN)))
master["CHEK1_FPKM"] = master._sid.map(dict(zip(tbl.sample_id, tbl.CHEK1)))
master.drop(columns="_sid").to_csv(f"{a.out}/master.csv", index=False)

print(f"wrote {a.out}/MYCN_FPKM.txt, FPKM_CN_MYCN.txt, master.csv")
print(f"MYCN FPKM  ecDNA+: median {tbl.MYCN[tbl.amp].median():.1f} (min {tbl.MYCN[tbl.amp].min():.1f}, max {tbl.MYCN[tbl.amp].max():.1f})"
      f" | ecDNA-: median {tbl.MYCN[~tbl.amp].median():.1f} (max {tbl.MYCN[~tbl.amp].max():.1f})")
print(f"MYCN CN    ecDNA+: median {tbl.CN[tbl.amp].astype(float).median():.0f} (min {tbl.CN[tbl.amp].astype(float).min():.0f}, max {tbl.CN[tbl.amp].astype(float).max():.0f})")

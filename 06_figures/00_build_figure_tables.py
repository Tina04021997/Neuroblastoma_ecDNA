"""Build the small per-line tables used by Figs 1A, 2A and 2B from scripted inputs.

Inputs
  samples.tsv            sample sheet (repo root)           sample_id, wgs_id, MYCN_amp, ecDNA_group, ...
  MYCN_CHEK1_FPKM.csv    from 03_export_FPKM_tables.R       Sample, ecDNA_group, MYCN, CHEK1  (StringTie FPKM)
  MYCN_CN.tsv (optional) frozen MYCN copy number per line   sample_id, MYCN_CN   (CNVkit; see data/README.md)
  master_in.csv          existing master.csv                WGS/ecDNA/chromothripsis columns are carried over

Outputs (written to --out)
  MYCN_FPKM.txt          Sample / MYCN / MYCN_amp                 -> fig2B_MYCN_FPKM_boxplot.py
  FPKM_CN_MYCN.txt       Sample / GE / CN / MYCNamp               -> fig2A_MYCN_CN_vs_FPKM.py
  master.csv             master_in.csv with MYCN_FPKM, CHEK1_FPKM (and MYCN_CN) replaced -> fig1A_oncoprint.ipynb

Only the expression (and, if MYCN_CN.tsv is given, copy-number) columns are rebuilt; every other column of
master.csv is copied unchanged. Sample order of MYCN_FPKM.txt / FPKM_CN_MYCN.txt follows samples.tsv.
"""
import argparse, sys
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--samples", required=True)
ap.add_argument("--fpkm", required=True, help="MYCN_CHEK1_FPKM.csv from 03_export_FPKM_tables.R")
ap.add_argument("--cn", default=None, help="MYCN_CN.tsv (sample_id, MYCN_CN); if omitted, CN is taken from master_in.csv")
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

if a.cn:
    cn = pd.read_csv(a.cn, sep="\t")
    assert set(cn.sample_id) == set(samples.sample_id)
    cn_map = dict(zip(cn.sample_id, cn.MYCN_CN))
else:
    cn_map = dict(zip(master._sid, master.MYCN_CN))
    print("NOTE: MYCN_CN carried over from master_in.csv (no MYCN_CN.tsv given)", file=sys.stderr)

f = fpkm.set_index("Sample")
tbl = samples[["sample_id", "MYCN_amp", "ecDNA_group"]].copy()
tbl["MYCN"]  = tbl.sample_id.map(f.MYCN).round(2)
tbl["CHEK1"] = tbl.sample_id.map(f.CHEK1)
tbl["CN"]    = tbl.sample_id.map(cn_map)
tbl["amp"]   = tbl.MYCN_amp.eq("amp")

TF = {True: "TRUE", False: "FALSE"}                      # same spelling as the original hand-made tables
fmt_cn = lambda v: f"{v:g}"                              # 273 not 273.0; 2.1 stays 2.1

# Fig 2B input
pd.DataFrame({"Sample": tbl.sample_id, "MYCN": tbl.MYCN, "MYCN_amp": tbl.amp.map(TF)}) \
  .to_csv(f"{a.out}/MYCN_FPKM.txt", sep="\t", index=False)

# Fig 2A input (amplified lines first, as in the original file)
t2 = tbl.sort_values("amp", ascending=False, kind="stable")
pd.DataFrame({"Sample": t2.sample_id, "GE": t2.MYCN, "CN": t2.CN.map(fmt_cn), "MYCNamp": t2.amp.map(TF)}) \
  .to_csv(f"{a.out}/FPKM_CN_MYCN.txt", sep="\t", index=False)

# Fig 1A input
master["MYCN_FPKM"]  = master._sid.map(dict(zip(tbl.sample_id, tbl.MYCN)))
master["CHEK1_FPKM"] = master._sid.map(dict(zip(tbl.sample_id, tbl.CHEK1)))
master["MYCN_CN"]    = master._sid.map(cn_map)
master.drop(columns="_sid").to_csv(f"{a.out}/master.csv", index=False)

print(f"wrote {a.out}/MYCN_FPKM.txt, FPKM_CN_MYCN.txt, master.csv")
print(f"MYCN FPKM  ecDNA+: median {tbl.MYCN[tbl.amp].median():.1f} (min {tbl.MYCN[tbl.amp].min():.1f}, max {tbl.MYCN[tbl.amp].max():.1f})"
      f" | ecDNA-: median {tbl.MYCN[~tbl.amp].median():.1f} (max {tbl.MYCN[~tbl.amp].max():.1f})")

#!/usr/bin/env python3
"""Cross-check sample IDs, MYCN status, ecDNA group/count, MYCN copy number and MYCN / CHEK1 FPKM across every table that
feeds Figures 1A, 2A and 2B (audit point 12). Exits non-zero on any disagreement.

  python 06_figures/01_check_figure_tables.py --samples samples.tsv --data data/
Tables: samples.tsv, data/master.csv (wgs_id names), data/MYCN_CN.tsv, data/FPKM_CN_MYCN.txt (Fig 2A),
        data/MYCN_FPKM.txt (Fig 2B), data/MYCN_CHEK1_FPKM.csv (03_export_FPKM_tables.R output).
"""
import argparse, sys
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--samples", required=True)
ap.add_argument("--data", required=True)
a = ap.parse_args(); d = a.data.rstrip("/")

S  = pd.read_csv(a.samples, sep="\t").set_index("sample_id")
m  = pd.read_csv(f"{d}/master.csv")
K  = pd.read_csv(f"{d}/MYCN_CN.tsv", sep="\t").set_index("sample_id")
B  = pd.read_csv(f"{d}/FPKM_CN_MYCN.txt", sep="\t").set_index("Sample")
C  = pd.read_csv(f"{d}/MYCN_FPKM.txt", sep="\t").set_index("Sample")
D  = pd.read_csv(f"{d}/MYCN_CHEK1_FPKM.csv").set_index("Sample")
m["sid"] = m.Sample.map(dict(zip(S.wgs_id, S.index)))
if m.sid.isna().any():
    sys.exit(f"master.csv samples not in samples.tsv wgs_id: {m.Sample[m.sid.isna()].tolist()}")
A = m.set_index("sid")

ids = sorted(S.index)
for name, t in [("master.csv", A), ("MYCN_CN.tsv", K), ("FPKM_CN_MYCN.txt", B), ("MYCN_FPKM.txt", C), ("MYCN_CHEK1_FPKM.csv", D)]:
    if set(t.index) != set(ids):
        sys.exit(f"{name}: sample set differs from samples.tsv: {sorted(set(t.index) ^ set(ids))}")
A, K, B, C, D, S = [t.loc[ids] for t in (A, K, B, C, D, S)]
amp = S.MYCN_amp == "amp"
checks = [
    (f"{len(ids)} sample IDs identical across samples.tsv and 5 data tables", len(ids) == 21),
    ("MYCN status: master.csv == samples.tsv",            ((A.MYCN_status == 1) == amp).all()),
    ("MYCN status: FPKM_CN_MYCN.txt (Fig 2A) == samples.tsv", (B.MYCNamp.astype(bool) == amp).all()),
    ("MYCN status: MYCN_FPKM.txt (Fig 2B) == samples.tsv", (C.MYCN_amp.astype(bool) == amp).all()),
    ("ecDNA group: MYCN_CHEK1_FPKM.csv == samples.tsv",    (D.ecDNA_group == S.ecDNA_group).all()),
    ("ecDNA count: master.csv == samples.tsv",             (A.ecDNA_counts == S.ecDNA_count).all()),
    ("MYCN CN: master.csv == MYCN_CN.tsv",                 (A.MYCN_CN.round(1) == K.MYCN_CN.round(1)).all()),
    ("MYCN CN: FPKM_CN_MYCN.txt (Fig 2A) == MYCN_CN.tsv",  (B.CN.round(1) == K.MYCN_CN.round(1)).all()),
    ("MYCN FPKM: master.csv == Fig 2A == Fig 2B == MYCN_CHEK1_FPKM.csv",
        ((A.MYCN_FPKM.round(2) == B.GE.round(2)) & (A.MYCN_FPKM.round(2) == C.MYCN.round(2)) & (A.MYCN_FPKM.round(2) == D.MYCN.round(2))).all()),
    ("CHEK1 FPKM: master.csv == MYCN_CHEK1_FPKM.csv",     (A.CHEK1_FPKM.round(4) == D.CHEK1.round(4)).all()),
    ("cell_line (figure labels, Fig 1A/1B/2A): 21 unique non-empty names in samples.tsv",
        S.cell_line.notna().all() and (S.cell_line.str.strip() != "").all() and S.cell_line.nunique() == 21),
]
bad = 0
for label, ok in checks:
    print(("OK   " if ok else "FAIL ") + label); bad += (not ok)
if bad:
    sys.exit(f"{bad} cross-check(s) failed")
print(f"OK: all {len(checks)} cross-checks agree for {len(ids)} lines")

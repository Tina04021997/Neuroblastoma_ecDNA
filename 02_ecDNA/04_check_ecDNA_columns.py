#!/usr/bin/env python3
"""Cross-check the curated ecDNA columns of data/master.csv against AmpliconClassifier output of the final runs.

  ecDNA_counts        vs  <dna_id>_ecDNA_counts.tsv            (AmpliconClassifier; header-only file = 0)
                      and the sum of 'ecDNA_amplicons' in <dna_id>_amplicon_classification_profiles.tsv
  oncogene_on_ecDNA   vs  canonical oncogenes (is_canonical_oncogene == True) located on ecDNA features
                          (feature starts with 'ecDNA') in <dna_id>_gene_list.tsv; lines without a gene list
                          (no amplicon) must have an empty oncogene list

Exits non-zero on any mismatch (audit point 13). Usage:
  python 02_ecDNA/04_check_ecDNA_columns.py --master data/master.csv --samples samples.tsv --dir data/ampliconclassifier
"""
import argparse, csv, os, sys

ap = argparse.ArgumentParser()
ap.add_argument("--master", required=True)
ap.add_argument("--samples", required=True)
ap.add_argument("--dir", required=True, help="folder with <dna_id>_ecDNA_counts.tsv, *_amplicon_classification_profiles.tsv, *_gene_list.tsv")
a = ap.parse_args()

samples = list(csv.DictReader(open(a.samples), delimiter="\t"))
wgs2dna = {r["wgs_id"]: r["dna_id"] for r in samples}
master = list(csv.DictReader(open(a.master)))
assert len(master) == 21 and len(samples) == 21, "expected 21 lines"

def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

errors, rows = [], []
for r in master:
    wgs = r["Sample"]
    if wgs not in wgs2dna:
        errors.append(f"{wgs}: not in samples.tsv wgs_id column"); continue
    d = wgs2dna[wgs]
    # 1. ecDNA count
    cnt_file = os.path.join(a.dir, f"{d}_ecDNA_counts.tsv")
    if not os.path.isfile(cnt_file):
        errors.append(f"{wgs}: missing {cnt_file}"); continue
    cnt_rows = read_tsv(cnt_file)
    ac_count = int(cnt_rows[0]["ecDNA_count"]) if cnt_rows else 0
    prof_file = os.path.join(a.dir, f"{d}_amplicon_classification_profiles.tsv")
    prof = read_tsv(prof_file) if os.path.isfile(prof_file) else []
    prof_sum = sum(int(p["ecDNA_amplicons"]) for p in prof)
    if prof and prof_sum != ac_count:
        errors.append(f"{wgs}: ecDNA_counts.tsv = {ac_count} but classification profile sums to {prof_sum}")
    master_count = int(r["ecDNA_counts"])
    if master_count != ac_count:
        errors.append(f"{wgs}: master.csv ecDNA_counts = {master_count}, AmpliconClassifier = {ac_count}")
    # 2. oncogenes on ecDNA
    gl_file = os.path.join(a.dir, f"{d}_gene_list.tsv")
    if os.path.isfile(gl_file):
        gl = read_tsv(gl_file)
        ac_genes = sorted({g["gene"] for g in gl if g["is_canonical_oncogene"] == "True" and g["feature"].startswith("ecDNA")})
    else:
        if ac_count > 0:
            errors.append(f"{wgs}: ecDNA present but no {gl_file}")
        ac_genes = []
    master_genes = sorted({g.strip() for g in r["oncogene_on_ecDNA"].split(",") if g.strip()})
    if master_genes != ac_genes:
        errors.append(f"{wgs}: master.csv oncogene_on_ecDNA = {master_genes}, AmpliconClassifier = {ac_genes}")
    rows.append((wgs, d, master_count, ac_count, ",".join(master_genes) or "-", ",".join(ac_genes) or "-"))

print(f"{'Sample':10s} {'dna_id':8s} {'master':>6s} {'AC':>3s}  master_oncogenes | AC_oncogenes")
for row in rows:
    print(f"{row[0]:10s} {row[1]:8s} {row[2]:6d} {row[3]:3d}  {row[4]} | {row[5]}")
if errors:
    print("\nMISMATCHES:"); print("\n".join("  " + e for e in errors)); sys.exit(1)
print(f"\nOK: ecDNA_counts and oncogene_on_ecDNA agree with AmpliconClassifier for all {len(rows)} lines")

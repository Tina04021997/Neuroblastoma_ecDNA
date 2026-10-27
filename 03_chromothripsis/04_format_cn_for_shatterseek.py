#!/usr/bin/env python3
"""Convert CNVkit calls from AmpliconSuite into the 4-column CN table ShatterSeek expects.

Input : <sample>*_CNV_CALLS.bed written by PrepareAA/CNVkit on the duplicate-marked WGS BAM
        (e.g. CHP134_tumor_mkdp_CNV_CALLS.bed, KAN_mkdp_CNV_CALLS.bed); *_pre_filtered.bed copies are skipped.
Output: <sample>_cn.txt  -> chromosome (no "chr"), start, end, CN (rounded to integer); chr1-22 and chrX only.
        Also consumed by 06_figures/fig1B_chr2_copy_number.py.

Usage:
  python 04_format_cn_for_shatterseek.py --in /path/to/cnvkit_beds --out /path/to/shatterseek/CNs
  (--in may be repeated, e.g. one folder per sequencing batch; sub-folders are searched)
"""
import argparse, glob, os, re, sys

ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="indirs", required=True, action="append", help="folder(s) with *_CNV_CALLS.bed")
ap.add_argument("--out", required=True)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

CHROMS = {str(i) for i in range(1, 23)} | {"X"}
SUFFIX = re.compile(r"[._](tumor_)?(mkdp|cs\.rmdup)_CNV_CALLS\.bed$")   # -> sample id is what precedes it

beds = []
for d in a.indirs:
    beds += glob.glob(os.path.join(d, "**", "*_CNV_CALLS.bed"), recursive=True)
beds = sorted(b for b in set(beds) if "pre_filtered" not in os.path.basename(b))
if not beds:
    sys.exit("ERROR: no *_CNV_CALLS.bed files found under " + ", ".join(a.indirs))

seen = {}
for bed in beds:
    name = os.path.basename(bed)
    if not SUFFIX.search(name):
        print(f"WARNING: skipping {name}: cannot derive sample id (expected <sample>[_tumor]_mkdp_CNV_CALLS.bed)", file=sys.stderr)
        continue
    sample = SUFFIX.sub("", name)
    if sample in seen:
        sys.exit(f"ERROR: two CNV_CALLS files for {sample}: {seen[sample]} and {bed}")
    seen[sample] = bed
    out = os.path.join(a.out, f"{sample}_cn.txt")
    n = 0
    with open(bed) as fi, open(out, "w") as fo:
        fo.write("chromosome\tstart\tend\tCN\n")
        for line in fi:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                continue
            chrom = c[0].replace("chr", "")
            if chrom in CHROMS:
                fo.write("\t".join([chrom, c[1], c[2], str(round(float(c[4])))]) + "\n")
                n += 1
    print(f"{sample}: {n} segments -> {out}")

print(f"done: {len(seen)} samples")

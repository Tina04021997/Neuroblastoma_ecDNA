#!/usr/bin/env python3
"""Extract the MYCN copy number used in Figs 1A/2A from AmpliconSuite output.

For each sample listed in samples.tsv, in order of preference:
  1. AmpliconClassifier <sample>_gene_list.tsv: gene_cn of MYCN (max over amplicons)  -> source "AC"
     (= AmpliconArchitect's PredictedCopyCount of the MYCN-containing amplicon segment)
  2. AmpliconArchitect <sample>_amplicon*_graph.txt: highest PredictedCopyCount among sequence
     edges overlapping the MYCN locus                                                  -> source "AA"
  3. No amplicon covers MYCN (MYCN non-amplified lines): the CNVkit segment overlapping MYCN
     in <sample>.cs.rmdup_CNV_CALLS.bed (column 5), rounded to 1 decimal              -> source "CNVkit"

Output: MYCN_CN.tsv  (sample_id, MYCN_CN, source, feature, segment)

Usage:
  python 03_extract_MYCN_CN.py --aa-dir /path/to/ecDNA_analysis --samples samples.tsv --out data/MYCN_CN.tsv
AmpliconSuite writes one folder per sample under --aa-dir (01_run_ampliconsuite.sh); the graph files are
<aa-dir>/<sample>/<sample>_amplicon<N>_graph.txt and the CNVkit calls <aa-dir>/<sample>/<sample>.cs.rmdup_CNV_CALLS.bed.
"""
import argparse, glob, os, re, sys

MYCN = ("chr2", 15_940_550, 15_947_007)   # GENCODE v44 gene span, hg38

ap = argparse.ArgumentParser()
ap.add_argument("--aa-dir", required=True)
ap.add_argument("--samples", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

def overlaps(chrom, s, e):
    return chrom == MYCN[0] and s < MYCN[2] and e > MYCN[1]

def aa_cn(sample):
    best = None
    for g in sorted(glob.glob(os.path.join(a.aa_dir, sample, f"{sample}_amplicon*_graph.txt"))):
        for line in open(g):
            if not line.startswith("sequence"):
                continue
            _, p1, p2, cn = line.rstrip("\n").split("\t")[:4]
            c1, s = re.match(r"(chr\w+):(\d+)", p1).groups()
            c2, e = re.match(r"(chr\w+):(\d+)", p2).groups()
            if c1 == c2 and overlaps(c1, int(s), int(e)):
                cand = (float(cn), f"{os.path.basename(g)}:{c1}:{s}-{e}")
                if best is None or cand[0] > best[0]:
                    best = cand
    return best

def ac_cn(sample):
    best = None
    for gl in glob.glob(os.path.join(a.aa_dir, sample, "**", f"{sample}_gene_list.tsv"), recursive=True):
        for line in open(gl):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5 and f[3] == "MYCN":
                cand = (float(f[4]), f[2], f"{os.path.basename(gl)}:{f[1]}")
                if best is None or cand[0] > best[0]:
                    best = cand
    return best

def cnvkit_cn(sample):
    beds = glob.glob(os.path.join(a.aa_dir, sample, f"{sample}*_CNV_CALLS.bed"))
    beds = [b for b in beds if "pre_filtered" not in b] or beds
    for b in beds:
        for line in open(b):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5 and overlaps(f[0], int(f[1]), int(f[2])):
                return (round(float(f[4]), 1), f"{os.path.basename(b)}:{f[0]}:{f[1]}-{f[2]}")
    return None

samples = [l.split("\t")[0] for l in open(a.samples).read().splitlines()[1:] if l.strip()]
rows, missing = [], []
for s in samples:
    r = ac_cn(s)
    if r is not None:
        cn, feat, seg, src = round(r[0]), r[1], r[2], "AC"
    else:
        r = aa_cn(s)
        if r is not None:
            cn, feat, seg, src = round(r[0]), "amplicon", r[1], "AA"
        else:
            r = cnvkit_cn(s)
            if r is None:
                missing.append(s); continue
            cn, feat, seg, src = r[0], "none", r[1], "CNVkit"
    rows.append((s, cn, src, feat, seg))

with open(a.out, "w") as fh:
    fh.write("sample_id\tMYCN_CN\tsource\tfeature\tsegment\n")
    for s, cn, src, feat, seg in rows:
        fh.write(f"{s}\t{cn:g}\t{src}\t{feat}\t{seg}\n")
        print(f"{s}\t{cn:g}\t{src}\t{feat}\t{seg}")
if missing:
    print("ERROR: no MYCN copy number found for:", ", ".join(missing), file=sys.stderr); sys.exit(1)
print(f"wrote {a.out} ({len(rows)} lines)")

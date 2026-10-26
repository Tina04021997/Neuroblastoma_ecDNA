#!/usr/bin/env python3
"""Extract the MYCN copy number used in Figs 1A/2A from AmpliconSuite output.

Inputs (anywhere below --dir, searched recursively; a flat folder of copied files is fine):
  <sample>_gene_list.tsv        AmpliconClassifier per-gene output   -> MYCN gene_cn (max over amplicons)
  <sample>*_CNV_CALLS.bed       CNVkit calls (from PrepareAA)        -> CN of the segment overlapping MYCN,
                                                                        used only when no amplicon covers MYCN
Rule, checked against samples.tsv:
  MYCN_amp == "amp"     -> must have a MYCN row in its gene list (source "AC"); otherwise the script fails
  MYCN_amp == "no_amp"  -> CNVkit segment value, 1 decimal (source "CNVkit"); a MYCN gene-list row would be a
                           warning (it would mean an amplicon covers MYCN in a line labelled non-amplified)

Output: MYCN_CN.tsv  (sample_id, MYCN_CN, source, feature, file)  — read by 06_figures/00_build_figure_tables.py --cn

Usage:
  python 03_extract_MYCN_CN.py --dir /path/to/ac_files --samples samples.tsv --id-column dna_id --out data/MYCN_CN.tsv
  (--dir may be repeated; --id-column names the samples.tsv column matching the AmpliconSuite file names,
   because three lines carry different IDs in the WGS data: NB10=NBro, SKNBE2=SKNEB2, SKNrI=SKNYrI)
"""
import argparse, glob, os, sys

MYCN = ("chr2", 15_940_550, 15_947_007)   # GENCODE v44 gene span, hg38

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True, action="append", help="folder(s) holding the gene lists / CNV_CALLS beds")
ap.add_argument("--samples", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--id-column", default="sample_id",
                help="samples.tsv column holding the name used in the AmpliconSuite file names (e.g. dna_id)")
a = ap.parse_args()

def find(pattern):
    hits = []
    for root in a.dir:
        hits += glob.glob(os.path.join(root, "**", pattern), recursive=True)
    return sorted(set(hits))

def overlaps(chrom, s, e):
    return chrom == MYCN[0] and s < MYCN[2] and e > MYCN[1]

def ac_cn(sample):
    """(gene_cn, feature, file) for MYCN, max over amplicons; None if no gene list or no MYCN row"""
    best = None
    for gl in find(f"{sample}_gene_list.tsv"):
        for line in open(gl):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5 and f[3] == "MYCN":
                cand = (float(f[4]), f[2], f"{os.path.basename(gl)} ({f[1]})")
                if best is None or cand[0] > best[0]:
                    best = cand
    return best

def cnvkit_cn(sample):
    beds = [b for b in find(f"{sample}*_CNV_CALLS.bed") if "pre_filtered" not in b] or find(f"{sample}*_CNV_CALLS.bed")
    for b in beds:
        for line in open(b):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5 and overlaps(f[0], int(f[1]), int(f[2])):
                return (round(float(f[4]), 1), f"{os.path.basename(b)} {f[0]}:{f[1]}-{f[2]}")
    return None

# samples.tsv: sample_id in column 1, MYCN_amp somewhere in the header
hdr, *lines = [l.split("\t") for l in open(a.samples).read().splitlines() if l.strip()]
i_amp = hdr.index("MYCN_amp"); i_id = hdr.index(a.id_column)
samples = [(l[0], l[i_id], l[i_amp]) for l in lines]   # (sample_id, file-name id, amp)

rows, errors = [], []
for s, fid, amp in samples:
    ac = ac_cn(fid)
    if amp == "amp":
        if ac is None:
            errors.append(f"{s}: MYCN-amplified but no MYCN row in any {fid}_gene_list.tsv under --dir"); continue
        rows.append((s, round(ac[0]), "AC", ac[1], ac[2]))
    else:
        if ac is not None:
            print(f"WARNING: {s} is labelled no_amp but its gene list has MYCN (gene_cn {ac[0]:.1f}, {ac[1]}); "
                  f"using CNVkit value as labelled", file=sys.stderr)
        ck = cnvkit_cn(fid)
        if ck is None:
            errors.append(f"{s}: non-amplified but no {fid}*_CNV_CALLS.bed with a segment overlapping MYCN under --dir"); continue
        rows.append((s, ck[0], "CNVkit", "none", ck[1]))

with open(a.out, "w") as fh:
    fh.write("sample_id\tMYCN_CN\tsource\tfeature\tfile\n")
    for r in rows:
        fh.write("\t".join([r[0], f"{r[1]:g}", r[2], r[3], r[4]]) + "\n")
        print("\t".join([r[0], f"{r[1]:g}", r[2], r[3], r[4]]))
if errors:
    print("\n".join("ERROR: " + e for e in errors), file=sys.stderr); sys.exit(1)
print(f"wrote {a.out} ({len(rows)} lines)")

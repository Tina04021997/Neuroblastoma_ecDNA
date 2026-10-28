# REANALYSIS: featureCounts read-level → fragment-level counting (audit point 4)

Date: 2026-10-19. Original outputs are preserved unchanged; this file documents a re-analysis and its effect on the reported numbers.

## What changed and why

The 2024 run used `featureCounts v2.0.6` with `-p` only. Since Subread 2.0.2, `-p` declares paired-end input but counts each mate as a separate read unless `--countReadPairs` is given. MultiQC confirmed read-level counting: featureCounts Assigned ÷ STAR aligned pairs = 1.46–1.67 for all 21 lines (a fragment count cannot exceed 1.0). The re-count uses `-p --countReadPairs -B -C`, i.e. one count per properly paired, non-chimeric fragment.

Everything downstream (`08_merge_counts.sh` steps 1–5 with the gene-ID guard, `01_deseq2_ecDNA.R` with the point-2 fixes, `02_gsea_and_RS_genes.R`) was run unchanged on both matrices. The read-level re-run reproduces the original export exactly (498/374 DE genes, identical padj, identical Hallmark NES), which validates the corrected scripts before comparing.

## Per-sample totals

| Sample | Original (reads) | Re-count (fragments) | Ratio |
|---|---:|---:|---:|
| CHP134 | 37,532,784 | 19,413,659 | 1.933 |
| CHP212 | 40,803,255 | 21,172,293 | 1.927 |
| IMR32 | 37,611,931 | 19,431,473 | 1.936 |
| Kelly | 36,616,087 | 19,005,054 | 1.927 |
| LAN6 | 48,288,779 | 25,043,627 | 1.928 |
| NBLS | 40,401,191 | 20,889,793 | 1.934 |
| NGP | 35,231,173 | 18,194,442 | 1.936 |
| SKNAS | 41,611,253 | 21,522,789 | 1.933 |
| SKNSH | 33,007,583 | 17,113,585 | 1.929 |
| TR14 | 37,481,483 | 19,385,592 | 1.933 |
| KAN | 33,695,268 | 17,400,649 | 1.936 |
| KCNR | 34,501,146 | 17,736,927 | 1.945 |
| KANR | 24,512,984 | 12,577,678 | 1.949 |
| KPNYN | 22,659,735 | 11,682,824 | 1.940 |
| LAN1 | 25,847,000 | 13,255,519 | 1.950 |
| LANS | 22,905,669 | 11,801,257 | 1.941 |
| NB10 | 29,547,855 | 15,200,853 | 1.944 |
| NB1643 | 23,050,139 | 11,832,374 | 1.948 |
| NB69 | 29,292,912 | 15,107,273 | 1.939 |
| SKNBE2 | 32,663,124 | 16,778,404 | 1.947 |
| SKNrI | 28,016,046 | 14,391,170 | 1.947 |

## Differential expression (ecDNA+ vs ecDNA−, FDR < 0.05, |log2FC| > 1)

| Quantity | Original (reads) | Fragments |
|---|---:|---:|
| Genes exported (with symbol, deduplicated) | 35,278 | 35,278 |
| DE genes up | 498 | 486 |
| DE genes down | 374 | 404 |
| DE genes total | 872 | 890 |
| Overlap of DE sets | — | 816 (94% of original retained) |
| Pearson r of log2FC (all shared genes) | — | 0.958 |
| Median log2FC shift | — | +0.002 (5–95th pct −0.58 … +0.48) |
| MYCN log2FC / padj | 6.57 / 1.4e-17 | 6.56 / 3.1e-19 |

DE genes lost in fragment analysis (56): ABHD1, ACP3, AGBL3, AJUBA-DT, AK2P1, C21orf62, C6, CA4, CENPP, CKM, CLEC17A, DCAF8L2, FAM237A, FLNC-AS1, GAGE12H, GIP, GRAMD4P3, HTR4, ITGB1BP2, KLK14, LBX1, LINC01589, LINC01625, LOC100129291, LOC100287792, LOC101927164, LOC102546299, LOC613206, MAJIN, MBOAT4, MICA, NCAPD3, OR1F1, PCAT4, PCDH11X, PDE11A, PGM5-AS1, PQBP1, PRICKLE4, RN7SL456P, RPL21P93, RPL3P2, RPS15AP12, RPS6KA4, SH2D4A, SLC1A6, SNORA12, SNORD34, SNORD65, SPINT1-AS1, ST13P19, TAS2R5, TBL1Y, TLR6, TMEM132E-DT, USP18

DE genes gained in fragment analysis (74): ADD3-AS1, ALDH18A1, AMDHD1, APOL6, ARX, B3GNT4, B3GNT7, BAALC, BHLHE22, CTSS, CXADR, DACT1, DCLK3, DLK2, DNAJC5G, DUSP22, EML6, EPPK1, F11R, FANK1, FBXL18, FLJ20021, GABRA2, GABRG2, GSKIP, HSD17B14, KRT15, LINC01886, LOC105376805, MAGEA9B, MAP3K15, MARCHF8, MIR25, MTND1P23, ND5, ND6, NECAB2, NME1-NME2, OPN1SW, PDE1C, PDK4, PIAS3, PIP5KL1, POLR2J4, PON2, PPIC, PPP1R17, PSPN, PTPRJ, RASGRP2, RNF39, RNLS, RPL13AP26, SAMD9, SH3BP5, SHFL, SLC5A12, SLC6A6, SMIM32, SNHG28, SNORA40, SP5, SPATA17, SYCE1, TAC1, TAF1D, TGFB2, THPO, TMEM102, VWA5A, YBX3, ZNF233, ZNF575, ZNF92

Interpretation: fold changes are essentially unchanged (ratio ≈ 1.93 cancels between groups); the DE set changes at the margin because read-level counting roughly doubles counts and understates dispersion. Fewer genes are flagged as outliers (1,035 vs 1,439) and more are filtered as low-count (49% vs 44%).

## Replication-stress panel (29 genes)

| Gene | log2FC orig | padj orig | log2FC frag | padj frag |
|---|---:|---:|---:|---:|
| ATR | -0.128 | 8.15e-01 | -0.141 | 7.86e-01 |
| ATRIP | +0.117 | 8.71e-01 | +0.101 | 8.88e-01 |
| TOPBP1 | -0.041 | 9.50e-01 | -0.040 | 9.50e-01 |
| CHEK1 | +1.111 | 3.41e-02 | +1.113 | 3.10e-02 |
| CLSPN | +0.081 | 9.47e-01 | +0.087 | 9.40e-01 |
| TIMELESS | +0.423 | 4.89e-01 | +0.417 | 4.72e-01 |
| TIPIN | +0.160 | 8.26e-01 | +0.166 | 8.07e-01 |
| RAD17 | -0.278 | 2.99e-01 | -0.279 | 2.80e-01 |
| WEE1 | +0.721 | 1.42e-01 | +0.721 | 1.34e-01 |
| CHEK2 | +0.520 | 3.16e-01 | +0.510 | 3.06e-01 |
| RPA1 | +0.679 | 8.28e-02 | +0.680 | 7.70e-02 |
| RPA2 | +0.118 | 8.89e-01 | +0.120 | 8.81e-01 |
| RPA3 | +0.147 | 8.65e-01 | +0.143 | 8.60e-01 |
| RRM1 | +0.953 | 3.57e-02 | +0.953 | 3.13e-02 |
| RRM2 | +0.447 | 5.87e-01 | +0.448 | 5.64e-01 |
| MCM2 | +0.674 | 5.45e-01 | +0.686 | 5.04e-01 |
| MCM7 | +1.086 | 6.51e-04 | +1.084 | 6.20e-04 |
| CDC45 | +0.949 | 5.44e-02 | +0.955 | 5.09e-02 |
| PCNA | +0.689 | 3.00e-01 | +0.691 | 2.79e-01 |
| RFC1 | +0.230 | 3.60e-01 | +0.227 | 3.53e-01 |
| POLD1 | +0.247 | 7.56e-01 | +0.251 | 7.37e-01 |
| POLE | +0.172 | 8.47e-01 | +0.160 | 8.50e-01 |
| CDC6 | +0.943 | 1.13e-01 | +0.941 | 1.01e-01 |
| FEN1 | -0.183 | 7.96e-01 | -0.183 | 7.87e-01 |
| MRE11 | +0.803 | 1.48e-01 | +0.803 | 1.37e-01 |
| RAD51 | +0.476 | 4.62e-01 | +0.475 | 4.39e-01 |
| BRCA1 | +0.553 | 3.93e-01 | +0.554 | 3.69e-01 |
| BRCA2 | +0.173 | 7.83e-01 | +0.174 | 7.72e-01 |
| FANCD2 | +0.993 | 3.93e-02 | +0.980 | 3.94e-02 |

Mean log2FC of the panel: 0.443 (original) vs 0.441 (fragments). Genes passing FDR < 0.05 are the same in both: CHEK1, RRM1, MCM7, FANCD2. Fig. 2E colouring (FDR < 0.2) and stars are unchanged for every gene.

## Hallmark GSEA (fgsea multilevel, seed 42, msigdbr 26.1.0)

Significant pathways (FDR < 0.05): 18 original, 22 fragments; all 18 original pathways remain significant. Newly significant: HALLMARK_ESTROGEN_RESPONSE_LATE, HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY, HALLMARK_UNFOLDED_PROTEIN_RESPONSE, HALLMARK_UV_RESPONSE_DN.

| Pathway | NES orig | FDR orig | NES frag | FDR frag |
|---|---:|---:|---:|---:|
| MYC_TARGETS_V1 | 2.440 | 1.26e-16 | 2.420 | 7.13e-17 |
| E2F_TARGETS | 2.430 | 1.26e-16 | 2.410 | 1.82e-16 |
| G2M_CHECKPOINT | 2.316 | 1.46e-13 | 2.269 | 5.89e-13 |
| INTERFERON_GAMMA_RESPONSE | -2.084 | 1.26e-07 | -2.057 | 3.56e-08 |
| INTERFERON_ALPHA_RESPONSE | -2.048 | 5.51e-06 | -2.051 | 2.31e-06 |
| INFLAMMATORY_RESPONSE | -1.797 | 1.49e-04 | -1.775 | 1.59e-04 |
| COMPLEMENT | -1.783 | 1.92e-04 | -1.752 | 1.11e-04 |
| OXIDATIVE_PHOSPHORYLATION | 1.766 | 3.49e-04 | 1.757 | 1.59e-04 |
| COAGULATION | -1.752 | 1.69e-03 | -1.732 | 1.35e-03 |
| PROTEIN_SECRETION | -1.666 | 4.00e-03 | -1.634 | 4.14e-03 |
| BILE_ACID_METABOLISM | -1.659 | 3.71e-03 | -1.614 | 5.03e-03 |
| MTORC1_SIGNALING | 1.657 | 2.48e-03 | 1.673 | 8.71e-04 |
| TNFA_SIGNALING_VIA_NFKB | -1.605 | 3.69e-03 | -1.583 | 1.74e-03 |
| MYC_TARGETS_V2 | 1.521 | 3.47e-02 | 1.560 | 2.61e-02 |
| MITOTIC_SPINDLE | 1.507 | 1.27e-02 | 1.485 | 1.29e-02 |
| REACTIVE_OXYGEN_SPECIES_PATHWAY | 1.498 | 7.80e-02 | 1.528 | 2.78e-02 |
| IL6_JAK_STAT3_SIGNALING | -1.498 | 3.47e-02 | -1.517 | 4.37e-02 |
| EPITHELIAL_MESENCHYMAL_TRANSITION | -1.406 | 3.47e-02 | -1.394 | 2.61e-02 |
| UNFOLDED_PROTEIN_RESPONSE | 1.405 | 5.45e-02 | 1.417 | 4.37e-02 |
| PI3K_AKT_MTOR_SIGNALING | -1.385 | 5.45e-02 | -1.392 | 6.93e-02 |

The top-20 Hallmark pathways and their direction are unchanged; NES differ by < 0.05 for every pathway.

## Files

- `data/original_2024_readcounts/`: original 21-line read-level matrix, the DESeq2 re-run from the corrected script (identical to the published export) and Hallmark results.
- `data/deseq2_gsea/`: fragment-level matrix (`-p --countReadPairs -B -C`), DESeq2 results, Hallmark results, RS panel, archived gene-set membership.

## Decision

**Adopted: fragment-level results (decided 2026-10-19).** The manuscript (v7) reports the
fragment-level DESeq2 and GSEA results; the read-level outputs are retained in
`original_reads/` for the record and are not used in any figure or number.

Re-run performed by the authors on 2026-10-19 from `counts_21lines.txt`
(62,700 genes x 21 lines; column order = `heatmap_condition.csv` = `list.txt` order:
CHP134 CHP212 IMR32 Kelly LAN6 NBLS NGP SKNAS SKNSH TR14 | KAN KCNR KANR KPNYN LAN1 LANS
NB10 NB1643 NB69 SKNBE2 SKNrI; batch 1 = first 10, batch 2 = last 11). The authors'
re-run of `01_deseq2_ecDNA.R` and `02_gsea_and_RS_genes.R` is byte-identical to
`fragments/` (max |delta log2FC| = 0, max |delta padj| = 0, identical NES).

Manuscript changes (tracked in v7):
- Results / Fig 2C legend: 498 up / 374 down / 16,943 ns -> **486 / 404 / 16,245**.
- Results top hits: *SP8* and *SOX8* removed. Under read-level counting they appeared at
  log2FC ~ 23 (padj 1e-17/1e-18); under fragment counting they fall to log2FC ~ 8
  (padj 0.015 / 0.022) -- still significant but no longer top hits. *EYA2*, *MNX1*,
  *GAGE1*, *MYCNOS*, *MYCN* remain top hits.
- Hallmark: 22 sets at FDR < 0.05 (all 18 original sets retained; the five named in
  Results -- E2F, G2M, mTORC1, OXPHOS, mitotic spindle -- and the five named for ecDNA-
  -- IFN-gamma, IFN-alpha, inflammatory, complement, coagulation -- are unchanged).
- RS panel: CHEK1 (*), MCM7 (***), FANCD2 (*), RRM1 (*) significant in both; CDC45
  padj 0.054 -> 0.051, CDC6 0.113 -> 0.101. Discussion wording "most significantly
  induced" -> "most strongly induced".
- Methods: featureCounts sentence now states fragments were counted
  (`-p --countReadPairs -B -C`, Subread 2.0.6) and how the 21-sample matrix is built.

Figures 2C-E regenerated from the fragment-level outputs (now `data/deseq2_gsea/`) and placed in the figure file.

## Audit points 1–4 (Priority 1) — fixes and verification record (2026-10-19 to 2026-10-27)

**Point 1 — DELLY VCF filter.** Original `02_filter_delly_vcf.sh` used `awk '/^chr(…)\s.*\sPASS\s/'`; `\s` is a
GNU-awk extension (on BSD awk the same command keeps 0 of 12,738 PASS records of CHP134). Replaced by explicit
tab-delimited column tests (`$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X)$/ && $7 == "PASS"`), with a stop if a sample keeps no
record. Verification (2026-10-27): applied to the 21 DELLY 1.7.2 VCFs (fileDate 2026-02-13), the fixed filter
reproduces byte for byte the 21 `<sample>_sv.txt` files that were fed to ShatterSeek in the 2026 run (CHP134: 10,435 of
12,738 PASS records kept; the 2,303 dropped are on unplaced/decoy contigs, chrM and chrY; 7,516 non-PASS records
excluded), so variant counts, chromosome distribution and the ShatterSeek calls (17 positive lines, 30 events) are
unchanged by the fix.

**Point 2 — DESeq2 count-matrix input.** Original `01_deseq2_ecDNA.R` read `merged_counts.txt` (versioned Ensembl IDs)
and took groups from a separate `heatmap_condition.csv`. It now reads the version-stripped `counts_21lines.txt` written
by `08_merge_counts.sh --join`, strips any residual version suffix defensively, takes groups from `samples.tsv`, and
stops unless > 50 % of IDs map to a symbol. Mapping yield on the 21-line matrix: 35,377 of 62,700 Ensembl IDs
(56.4 %; the remainder are unannotated/novel GENCODE IDs with no org.Hs.eg.db symbol), 99 duplicate-symbol rows dropped
(first occurrence kept, GTF order). Re-run on the original read-level matrix reproduces the published 498 / 374 DE genes
with identical log2FC and padj; the fragment-level matrix gives 486 / 404 (point 4, table below).

**Point 3 — gene-ID-safe merge.** Original step 4 pasted count columns positionally. `08_merge_counts.sh` now runs
`08b_featurecounts_sanity_check.sh` itself before pasting; 08b writes `reference_column.txt` and compares every sample's
gene-ID column with it (content and order), and the merge aborts on any mismatch. `--join` additionally requires
identical gene-ID columns between the two batch matrices and no versioned IDs. Verification (2026-10-27, synthetic
featureCounts files): well-formed inputs merge; one sample with shuffled rows → "MISMATCH … gene IDs differ", exit 1,
nothing merged; one sample with a missing gene → exit 1.

**Point 4 — paired-end featureCounts.** Documented in the first section of this file: `-p` alone counted reads in
Subread 2.0.6; re-counted with `-p --countReadPairs -B -C`; all downstream numbers re-derived and compared below.

**Point 12 — cross-check across figures.** `06_figures/01_check_figure_tables.py` (also run by `make_figures.sh`)
asserts that sample IDs, MYCN status, ecDNA group and count, MYCN copy number and MYCN / CHEK1 FPKM agree across
`samples.tsv`, `master.csv` (Fig 1A), `FPKM_CN_MYCN.txt` (Fig 2A), `MYCN_FPKM.txt` (Fig 2B), `MYCN_CN.tsv` and
`MYCN_CHEK1_FPKM.csv`: 10 of 10 checks pass for all 21 lines (2026-10-27); a flipped label is caught (exit 1).

## Audit point 5 — duplicate handling and StringTie modes (2026-10-23)

**Chains actually executed**

| Batch | Lines | RNA root (TSCC) | StringTie input | StringTie mode | featureCounts input |
|---|---|---|---|---|---|
| 1 | 10 (CHP134 … TR14) | `…/CureBound/Pilot/RNA/20230717/RNA_analysis` | `markDP/<s>.bam` (Picard MarkDuplicates 3.0.0, flag only) | reference-guided assembly (no `-e`) | raw STAR PASS2 BAM |
| 2 | 11 (KAN … SKNrI) | `…/CureBound/Second/RNA/rna_analysis` | raw STAR PASS2 BAM (no Picard step) | reference-only (`-e`) | raw STAR PASS2 BAM |

GTF for both: `…/Databases/STAR/GRCh38/gencode.v44.basic.annotation.gtf`. `--ignoreDup` never used; Picard
`REMOVE_DUPLICATES` left at default (false). No quantification in the paper is duplicate-filtered.

**Duplicate-flag test (CHP212, batch 1).** StringTie 2.2.1 was re-run with the batch-1 command on the raw STAR BAM
and compared with the existing output made from the Picard-flagged BAM: 103,574 rows each; all annotated genes
(`diff` of sorted ENSG rows = 0) and all 40,816 novel loci (ID column ignored, diff = 0) identical in coverage,
FPKM and TPM. Only the `STRG.N` numbering order differed (thread scheduling). Conclusion: StringTie ignores flag
0x400; the Picard step had no effect on batch-1 FPKMs. Picard metrics, CHP212: 25,410,660 read pairs examined,
PERCENT_DUPLICATION 0.687.

**StringTie mode difference.** Batch 1 gene.abundance files contain 40,816 `STRG.` novel genes; batch 2 files
contain none. FPKM denominators therefore differ between batches. Decision: not re-run; disclosed in Methods (v7)
and to be recorded per sample in `samples.tsv` (`stringtie_mode`).

**FPKM table provenance.** The MYCN FPKM values in `master.csv`, `MYCN_FPKM.txt` and `FPKM_CN_MYCN.txt`
(e.g. Kelly 526.89) do not match any StringTie output that still exists on TSCC or locally (Kelly 493.44 in the
only Kelly abundance file; CHEK1 values in `master.csv` *do* match the existing files). Decision: rebuild all three
tables from the existing `gene.abundance.txt` files with `03_export_FPKM_tables.R` and regenerate Figs 1A/2A/2B.

**Picard MarkDuplicates metrics, batch 1 (RNA; Picard 3.0.0; flag only, nothing removed).**
Batch 2 was not passed through Picard, so no metrics exist for it.

| Line | Read pairs examined | PERCENT_DUPLICATION |
|---|---|---|
| CHP134 | 23,634,124 | 0.689 |
| CHP212 | 25,410,660 | 0.687 |
| IMR32 | 23,933,616 | 0.757 |
| Kelly | 23,143,120 | 0.606 |
| LAN6 | 30,714,451 | 0.723 |
| NBLS | 25,267,184 | 0.622 |
| NGP | 22,072,373 | 0.586 |
| SKNAS | 25,485,855 | 0.720 |
| SKNSH | 22,354,703 | 0.838 |
| TR14 | 23,668,345 | 0.685 |

Mean 0.691 (range 0.586–0.838). High rates are expected for deep poly(A) RNA-seq dominated by highly expressed
genes and are not evidence of library failure; duplicates were retained in all quantifications (see Methods).

## Audit points 11/12 — FPKM table rebuild (2026-10-23)


**FPKM table rebuild — done 2026-10-23 (audit points 11/12).** `03_export_FPKM_tables.R` (21 × `gene.abundance.txt`
in `Pilot_Second/RS/`, groups from `samples.tsv`) → `MYCN_CHEK1_FPKM.csv` → `00_build_figure_tables.py` →
`MYCN_FPKM.txt`, `FPKM_CN_MYCN.txt`, `master.csv` (`master_sheet/2026_rerun/`; CN column carried over from the 2024
sheet). Only `MYCN_FPKM` changed (CHEK1 already matched): ecDNA+ median 426.6 (208.2–941.1), ecDNA− median 2.0
(max 19.1). Perfect group separation before and after → Fig 2B Mann–Whitney p unchanged. Figs 1A, 2A, 2B
regenerated from the new tables (`MYCN_ecDNA/2026_rerun/`). Results sentences ("FPKM exceeding 800", "median over
400 FPKM") remain correct; no manuscript text change required. The 2024 tables (Kelly MYCN 526.89 etc.) came from a
StringTie run that no longer exists and are superseded.

## MYCN copy-number provenance — Figs 1A/2A (2026-10-26)

**Finding.** The `MYCN_CN` column of `master.csv` (hand-entered 2024) matched the AmpliconClassifier `gene_cn` of a
*superseded* AmpliconSuite run (PrepareAA `--fastqs … --cnvkit_segmentation hmm-tumor`, old /oasis TSCC), not the final
runs (`Pilot/DNA/ecDNA_analysis_cbs/`, `Second/DNA/ecDNA_analysis/`: `--bam <sample>_tumor_mkdp.bam`, CNVkit CBS,
`--cngain 4.5 --cnsize_min 50000 --downsample -1`). CHP134 gave 273 in the old run and 275.7 in the final run.

**What the number is.** For MYCN-amplified lines it is the AmpliconClassifier `gene_cn` of MYCN (= AmpliconArchitect
`PredictedCopyCount` of the MYCN-containing segment, max over amplicons); for the six non-amplified lines it is the
CNVkit segment overlapping MYCN. It is not a CNVkit call for amplified lines, which the Fig 2A legend previously implied.

**Script.** `02_ecDNA/03_extract_MYCN_CN.py --dir <AmpliconSuite output> --samples samples.tsv --id-column dna_id`
(three lines carry different WGS IDs: NB10 = NBro, SKNBE2 = SKNEB2, SKNrI = SKNYrI) writes `data/MYCN_CN.tsv`;
`06_figures/00_build_figure_tables.py --cn` rebuilds `FPKM_CN_MYCN.txt` and `master.csv`.

**Result (final run vs 2024 table).** CHP134 273→276, CHP212 103→104, NGP 282→285, TR14 394→403; LAN6 2→2.9, NBLS 2.1→3,
SKNAS 2→1.9, SKNSH 1.9→2.9, NB69 2→2.7, SKNrI 2→2.2; 11 amplified lines unchanged. ecDNA+ MYCN CN median 268 (17–505).
Non-amplified values of 2.7–3.0 reflect the 2p copy state of the CNVkit segment. Fig 2A and Fig 1A re-plotted from the
rebuilt tables; no visible change. Manuscript: Methods provenance sentence and Fig 2A/1B legends (v8), NGP example
282→285 and "reflecting 2p copy state" (v9). Original tables kept in `FPKM_2026_rerun/before_CN_update/`.

**ecDNA columns cross-checked (2026-10-27).** The `ecDNA_counts` and `oncogene_on_ecDNA` columns of `master.csv` were
compared with the AmpliconClassifier output of the final runs (batch 1 CBS run, files dated 2024-02-25; batch 2 run, files
dated 2024-02-12), now archived in `data/ampliconclassifier/`: `<dna_id>_ecDNA_counts.tsv` (ecDNA per line),
`*_amplicon_classification_profiles.tsv` (class and ecDNA count per amplicon) and `*_gene_list.tsv` (canonical oncogenes on
ecDNA features). `02_ecDNA/04_check_ecDNA_columns.py` asserts both columns and exits non-zero on a mismatch; all 21 lines
agree (TR14 5, NGP 2, LAN1 2, twelve amplified lines 1, six non-amplified lines 0; oncogene sets identical). Note that the
AmpliconArchitect `#Amplicons` value counts amplicons of every class, so it is not the ecDNA count (TR14: 4 amplicons,
5 ecDNA; SKN-BE2: 2 amplicons, 1 ecDNA; NB69: 1 linear amplicon, 0 ecDNA); the notebook no longer warns on that
difference. `data/ecDNA_summary/` was confirmed identical to the final-run summaries for both batches. The Fig 1B
chromosome-2 profile was not re-derived (user decision 2026-10-26).

## Audit point 6 — software versions (2026-10-26)

**Reconciliation.** The AmpliconSuite run records (`<sample>_run_metadata.json`; batch 2 KAN launched 2024-02-09,
batch 1 CHP134 re-run 2024-02-19) give CNVkit **0.9.10**, AmpliconSuite-pipeline 1.2.2, AmpliconClassifier 1.1.1;
the AA log gives AmpliconArchitect 1.3.r7, Python 3.10.13. The README value 0.9.10 was correct; the manuscript value
0.9.12 came from the current conda environment (`ampsuite`, exported 2026-10-26: ampliconsuite 1.4.0, cnvkit 0.9.12),
which was updated after the analysis. Manuscript corrected 0.9.12→0.9.10 (v10, tracked change, comment 85).

**Other versions.** DELLY 1.7.2 static binary (VCF headers, run 2026-02-13); ShatterSeek 1.1 in env ASCAT (R 4.2.2);
CycleViz 0.2.2 (`_version.py`); RNA tools from `rna-seq.yml` (STAR 2.7.10b, StringTie 2.2.1, Subread 2.0.6,
Picard 3.0.0, Trim Galore 0.6.10, FastQC 0.12.1); WGS from `evc_main.yml` (bwa 0.7.17, Picard 2.18.27, GATK 4.1.4.1);
local R 4.2.2 with DESeq2 1.36.0, fgsea 1.22.0, msigdbr 26.1.0, org.Hs.eg.db 3.15.0. All match the Methods.

**Files.** `envs/{rna-seq,evc_main,ampsuite,ASCAT}.yml`, `envs/R_sessionInfo_local.txt`, `envs/versions.md`,
`envs/run_records/`. Figure-script Python package versions are not pinned (printed at run time).

## ecDNA summary files for Fig 1A — replaced for batch 1 (2026-10-26)

**Finding.** The six batch 1 AmpliconArchitect summary files (`CHP134, CHP212, IMR32, Kelly, NGP, TR14`) that the Fig 1A
notebook reads (`data/ecDNA_summary/<sample>_summary.txt`) were from the superseded hmm-tumor run; the other 15 were already
from the final runs. They were replaced by `<sample>_AA_results/<sample>_summary.txt` from `Pilot/DNA/ecDNA_analysis_cbs/`
(old copies kept in `master_sheet/2026_rerun/ecDNA_summary_old_hmm_batch1/`).

**Differences (old → new).** Amplicon numbers and MYCN-amplicon oncogene content unchanged for all six lines; average
amplified copy counts shift by a few percent (e.g. CHP134 191.0→180.3, Kelly 221.2→202.1, NGP amplicon1 107.2→105.1).
TR14 amplicon 2 no longer carries SNX29 (new run: no oncogene on that amplicon); the other TR14 amplicons are unchanged
(amplicon 1 AGAP2/AGRN/CDK4/MYCN, amplicon 4 MDM2/YEATS4). Gene order within the lists differs (alphabetical in the new run).

**Effect.** `master.csv` `ecDNA_counts` / `oncogene_on_ecDNA` validate against the new summaries (notebook ran through);
Fig 1A regenerated 2026-10-26 16:43: identical except that the SNX29 row is gone. SNX29 is not mentioned in the manuscript
text or the Fig 1 legend, so no text change is needed. `data/ecDNA_summary/` in the repository now holds the final-run files.

## Audit points 8, 13, 15 — sample sheet, fail-fast scripts, config.sh (2026-10-26)

**Point 8.** `05_differential_expression/01_deseq2_ecDNA.R` now reads the sample groups from `samples.tsv` (`MYCN_amp`,
levels `no_amp` < `amp`) instead of a separate `heatmap_condition.csv`, checks that every count column is present, that
every label is valid, that the group sizes are exactly 6 `no_amp` / 15 `amp`, and that `MYCN_amp` and `ecDNA_group`
coincide for all 21 lines (`stopifnot`; the design uses `MYCN_amp`, the Methods describe ecDNA+ vs ecDNA−, and in this
cohort the two contrasts are the same). It writes `deseq2_sessionInfo.txt`. Re-run locally on `counts_21lines.txt` with
the assertions in place (2026-10-27): `deseq2_results_ecDNA.csv` byte-identical to the 2026-10 re-run (486 up / 404 down).
**Points 13/15.** All 15 shell/PBS scripts now `set -euo pipefail`, source `config.sh` (template `config.example.sh`,
git-ignored copy holds site paths), check their inputs with `require_file`, take sample lists from `samples.tsv` or the
per-batch `$RNA_LIST`, and log tool versions (`log_version`). No script contains an absolute path. `06_stringtie.sh`
(PBS, batch 1) was replaced by `06_stringtie_batch1.sh` and a new `06_stringtie_batch2.sh` reproducing the batch 2
SLURM job (`-e`, unflagged STAR BAM). `08_merge_counts.sh --join` builds the 21-line matrix from `$BATCH1_COUNTS` /
`$BATCH2_COUNTS`. Resource requests and tool flags are unchanged from the executed jobs.

## Required final validation — before / after (2026-10-26)

"Before" = original 2024 analysis (read-level featureCounts, hand-entered tables; `data/original_2024_readcounts/`).
"After" = frozen inputs in `data/` and `data/deseq2_gsea/` (fragment-level counts, scripted tables).

| Item | Before (2024) | After (2026) | Note |
|---|---|---|---|
| Cell lines | 21 | 21 | identical IDs (`samples.tsv`) |
| Group membership | 15 MYCN-amp/ecDNA+ : 6 non-amp/ecDNA− | 15 : 6 | identical lines; `MYCN_amp` and `ecDNA_group` coincide for all 21 (point 8) |
| DE genes (FDR < 0.05, \|log2FC\| > 1) | 498 up / 374 down (17,815 tested) | 486 up / 404 down (17,135 tested) | fragment counting (point 4); log2FC Pearson r = 0.958 over 35,278 shared genes |
| MYCN expression, ecDNA+ vs ecDNA− | median 413.8 FPKM; one-sided Mann-Whitney p = 1.8e-05 | median 426.6 FPKM; p = 1.8e-05 | FPKM table rebuilt from StringTie output; ranks unchanged |
| Replication-stress panel | 29 genes ("n = 30" string in script/legend) | 29 genes (`data/deseq2_gsea/rs_genes.txt`) | point 9; `02_gsea_and_RS_genes.R` and `03_export_FPKM_tables.R` now `stopifnot` 29 unique symbols and that all 29 are present in the DESeq2 results (re-run 2026-10-27: outputs byte-identical) |
| Hallmark GSEA | 18 of 50 sets FDR < 0.05 | 22 of 50 | NES Pearson r = 0.9998; same top sets (MYC targets V1, E2F, G2M up; IFN-γ/α down); 4 sets newly below FDR 0.05 (estrogen response late, ROS pathway, UPR, UV response down); none lost |
| Chromothripsis calls per line | 17 positive (30 high-confidence events) | 17 positive (30 events) | `Chromothripsis_status` / `_counts` of `master.csv` reproduced exactly by `06_call_high_confidence_chromothripsis.ipynb` run on the 21 ShatterSeek summaries of the 2026 DELLY/ShatterSeek run (2026-10-27); the earlier "15 positive" in this row was a transcription error |
| ecDNA counts / SV burden per line | — | — | identical to 2024 values; ecDNA counts and oncogenes on ecDNA asserted against final-run AmpliconClassifier output (`04_check_ecDNA_columns.py`, 21/21 agree) |
| MYCN copy number | hand-entered from superseded run | scripted from final run | 10 lines differ by ≤ 2 % (amplified) or < 1 copy (non-amplified) |

## Code review before release (2026-10-27)

An independent read of the day's changes, plus functional tests, found and fixed the following (all verified by
re-running; no result changed):

- `.gitignore`: a trailing comment on the `figures/` line made the pattern inert; fixed, and the notebook's temporary
  `06_figures/neuroblastoma_heatmap.*` is ignored too.
- `06_figures/make_figures.sh`: an R failure could be masked by `|| true` and a previous run's files re-used; the notebook
  step could fail silently. Every step now logs to `<out_dir>/logs/`, stops on failure with the log tail, clears the GSEA
  folder before re-running, and the notebook kernel is selectable (`JUPYTER_KERNEL`).
- All 13 job scripts: `for s in $(sample_ids); do` did not stop on a missing sample list (a failed substitution inside a
  `for` word list is not caught by `set -e`); replaced by `ids=$(sample_ids) || exit 1`, which is.
- `config.example.sh`: unquoted `<partition>` placeholders made the template unparsable; quoted, and a guard now stops any
  script sourcing a `config.sh` in which any value is still a placeholder.
- `03_export_FPKM_tables.R`: takes its four paths as arguments (README usage corrected), stops on placeholders, no longer
  duplicates `MYCN_FPKM.txt` (single source: `00_build_figure_tables.py`), and sorts ecDNA+ rows first independent of the
  locale (`data/CHEK1_FPKM.csv` re-written in that order; same rows).
- `02_gsea_and_RS_genes.R`, `fig2C_volcano.py`: stop with a usage message instead of failing on the placeholder path.
- `fig1B_chr2_copy_number.py`: exits non-zero on a bad input instead of printing the error and returning 0; locus comment
  corrected to hg38.
- Functional tests run: `04_check_ecDNA_columns.py` exits 1 on a corrupted `master.csv`; the DESeq2 group assertions fail
  on a wrong label and on a `MYCN_amp` / `ecDNA_group` disagreement; `00_build_figure_tables.py` and
  `03_export_FPKM_tables.R` reproduce the frozen tables byte for byte; `make_figures.sh` produces all panels and the
  re-computed Hallmark table is identical to the frozen copy.

## Audit point 15 — single reproducible workflow (2026-10-27)

- **No placeholder or hard-coded path remains in any script.** Cluster scripts take everything from `config.sh`
  (template `config.example.sh`, which refuses to be used unedited); local scripts take their paths as command-line
  arguments and stop with a usage message if a default placeholder is left. Converted today in addition to the figure
  scripts: `01_deseq2_ecDNA.R <deseq_dir> <samples.tsv>`, `03_format_sv_for_shatterseek.py --in --out` (also fails on
  a malformed record instead of raising an IndexError), `05_run_shatterseek.R <sample> <sv_dir> <cn_dir> <results_dir>
  [chromosomes]` (the five hand-copied per-chromosome plot blocks became a loop over the same chromosomes 1, 2, 7, 3, 12;
  `<sample>_summary.txt` is now written to `<results_dir>` instead of the working directory), and the chromothripsis
  notebook reads `$SHATTERSEEK_RESULTS`. `01_deseq2_ecDNA.R` re-run through its arguments: results byte-identical.
  Testing `03_format_sv_for_shatterseek.py` on the 2026 DELLY files showed that the repository copy differed from the
  script actually executed on the cluster (`process_sv.py`): the executed version sets `chr2 = chr` when a record has
  no `CHR2` field (DELLY 1.x omits it for intra-chromosomal DEL/DUP/INV), whereas the repository copy would have raised
  an IndexError on the first such record. The repository script now carries the executed logic; its output for all 21
  lines is byte-identical to the `sv_<sample>.txt` files used in the 2026 ShatterSeek run.
  `06_call_high_confidence_chromothripsis.ipynb` executed on the 21 `<sample>_summary.txt` files of that run reproduces
  `Chromothripsis_status` and `Chromothripsis_counts` of `master.csv` for all 21 lines (17 positive, 30 events); the
  notebook now stops if the results folder is unset or empty. A second review also found that the three WGS job
  generators (`01_run_ampliconsuite.sh`, `01_delly_call.sh`, `02_filter_delly_vcf.sh`) defaulted to the `sample_id`
  column of `samples.tsv`, whereas the WGS files are named by `dna_id` (three lines differ); they now use a `wgs_ids`
  helper that reads the `dna_id` column by header name.
- **Inputs, commands, environments and output locations** are tabulated per stage in the README *Workflow* section;
  run-time versions with evidence in `envs/versions.md`; conda exports in `envs/*.yml`; local R `sessionInfo` and the
  Python package set used for figure regeneration (`envs/figures_python.txt`).
- **Data access** is documented in the README *Data availability* section: raw sequencing data and per-sample pipeline
  outputs stay on the cluster; everything needed for stages 5–6 (DESeq2 → figures) is archived under `data/`, and those
  stages were run cleanly from the repository alone on 2026-10-27 (`make_figures.sh`).

## Audit point 11 — frozen figure inputs and a documented regeneration command (2026-10-27)

All inputs behind Figures 1 and 2 are in `data/`: `master.csv`, `FPKM_CN_MYCN.txt`, `MYCN_FPKM.txt`, `MYCN_CN.tsv`,
`MYCN_CHEK1_FPKM.csv`, the AmpliconArchitect summaries (`ecDNA_summary/`), the AmpliconClassifier output
(`ampliconclassifier/`), the chromosome-2 copy-number tables for NGP and SKNFI (`chr2_cn/`, from
`04_format_cn_for_shatterseek.py` on the final CNVkit calls; NGP MYCN segment CN = 238), the fragment-level count matrix
and the final DESeq2 / GSEA results (`deseq2_gsea/`), plus the original 2024 read-level counts and results
(`original_2024_readcounts/`). `heatmap_condition.csv` is superseded by `samples.tsv`.

`06_figures/make_figures.sh` regenerates every panel from these files with one command. The figure scripts and
`02_gsea_and_RS_genes.R` now take their input and output paths as arguments (defaults unchanged), so no placeholder
editing is needed. Verified 2026-10-27 (Python 3.9.13, pandas 2.2.3, matplotlib 3.9.4, scipy 1.13.1, adjustText 0.8;
R 4.2.2, fgsea 1.22.0, msigdbr 26.1.0): all nine panel files were produced; Fig 2C reports 486 / 404 / 16,245; the
re-computed `hallmark_GSEA_results.csv` is byte-identical to the frozen copy; the RS panel is 29 / 29 with CHEK1, MCM7,
FANCD2 and RRM1 below FDR 0.05; the ecDNA-column check passes 21 / 21. Panels match the manuscript figure file
(`Figure_rerun.pptx`, checked panel by panel on 2026-10-27).

## Audit point 10 — Fig 2E FDR encoding (2026-10-27)

Fig 2E uses two encodings: bar colour at FDR < 0.2 (red higher in ecDNA+, blue higher in ecDNA−, grey otherwise) and
asterisks at FDR < 0.05 / 0.01 / 0.001. The authors chose to keep the dual encoding (audit option "keep only if
intentional and clearly described"): the 0.2 colouring is a deliberately lenient, exploratory cut that shows the direction
of change across the whole 29-gene panel and is not a significance claim; significance is carried by the asterisks only.
Manuscript v11 states this explicitly in the Methods and in the Fig 2E legend (tracked changes, comments 86–87), and the
comment in `02_gsea_and_RS_genes.R` was corrected to say the same. No value, colour or star changed: 8 genes are below
FDR 0.2 and CHEK1 (*), MCM7 (***), FANCD2 (*) and RRM1 (*) are below 0.05, identical in the original and re-run results.

## Audit point 14 — Figure 3–4 analysis code: not addressed (2026-10-26)

The functional-assay analyses (R2 survival grouping, DepMap processing, dose-response / AUC, caspase and immunoblot
quantification) were not performed by the maintainer of this repository and no scripts for them are included.
The README states that Figures 3 and 4 are not generated from this repository. This point remains open.

## Audit point 7 — Hallmark gene sets pinned and archived (2026-10-26)

`02_gsea_and_RS_genes.R` uses `msigdbr(collection = "H")`, stops unless msigdbr >= 26.1.0 (MSigDB 2026.1.Hs), writes
`gsea_session_versions.txt`, and now also writes the gene-set membership it used
(`hallmark_genesets_msigdbr_<version>.tsv`: gs_name, gene_symbol, ensembl_gene). The archived file for the published run,
`data/deseq2_gsea/hallmark_genesets_msigdbr_26.1.0.tsv` (50 sets, 7,331 set-gene rows), was generated with the same
msigdbr 26.1.0 installation that produced `hallmark_GSEA_results.csv`.

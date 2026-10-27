# ecDNA-MYCN-neuroblastoma

Analysis and figure code for:

> Yang T, *et al.* **Extrachromosomal *MYCN* Amplification Defines a Replication Stress-Addicted Neuroblastoma Subset with Selective Vulnerability to CHK1 Inhibition.**

This repository covers the whole-genome sequencing (WGS) and RNA-seq analyses of 21 neuroblastoma cell lines behind **Figures 1 and 2**:
- ecDNA detection
- chromothripsis calling
- differential expression and gene set enrichment
- figure generation

Figures 3 and 4 (functional assays) are not generated from this repository.

## Figure-to-code map

| Panel | Data tables | Plot |
|---|---|---|
| 1A | `02_ecDNA/03_extract_MYCN_CN.py` + `05_differential_expression/03_export_FPKM_tables.R` → `06_figures/00_build_figure_tables.py` → `data/master.csv`; ecDNA columns validated against `data/ecDNA_summary/` | `06_figures/fig1A_oncoprint.ipynb` |
| 1B (ecDNA cycle, NGP) | AmpliconArchitect cycles/graph (`02_ecDNA/01_run_ampliconsuite.sh`) | `02_ecDNA/02_run_cycleviz.sh` |
| 1B (chr2 copy number) | CNVkit calls → `03_chromothripsis/04_format_cn_for_shatterseek.py` | `06_figures/fig1B_chr2_copy_number.py` |
| 2A | as 1A → `data/FPKM_CN_MYCN.txt` | `06_figures/fig2A_MYCN_CN_vs_FPKM.py` |
| 2B | `03_export_FPKM_tables.R` → `00_build_figure_tables.py` → `data/MYCN_FPKM.txt` | `06_figures/fig2B_MYCN_FPKM_boxplot.py` |
| 2C | `05_differential_expression/01_deseq2_ecDNA.R` | `06_figures/fig2C_volcano.py` |
| 2D, 2E | `05_differential_expression/02_gsea_and_RS_genes.R` | same script |

Every number plotted in Figures 1A, 2A and 2B therefore traces to a file written by a pipeline script; see `data/` and
`REANALYSIS.md` for how the tables were rebuilt.

## Repository structure

```
├── README.md
├── REANALYSIS.md                        # what was re-run during the 2026 code audit, and what changed
├── samples.tsv                          # 21 lines: sample_id, cell_line, rna_batch, MYCN_amp, ecDNA_group, wgs_id, dna_id, stringtie_mode
├── config.example.sh                    # copy to config.sh (git-ignored): paths, scheduler, conda; sourced by every script
├── data/                                # scripted tables behind Figs 1A/2A/2B
│   ├── MYCN_CN.tsv                      # 02_ecDNA/03_extract_MYCN_CN.py
│   ├── MYCN_CHEK1_FPKM.csv, CHEK1_FPKM.csv   # 05_differential_expression/03_export_FPKM_tables.R
│   ├── MYCN_FPKM.txt, FPKM_CN_MYCN.txt, master.csv   # 06_figures/00_build_figure_tables.py
│   └── ecDNA_summary/<sample>_summary.txt   # AmpliconArchitect summaries (final runs); read by fig1A_oncoprint.ipynb
├── envs/                                # conda exports, R sessionInfo, run records, versions.md
├── 01_wgs_alignment/
│   └── bwa_mem_markduplicates.pbs       # BWA-MEM → samtools sort → Picard MarkDuplicates (env evc_main)
├── 02_ecDNA/
│   ├── 01_run_ampliconsuite.sh          # PrepareAA: CNVkit (CBS) → AmpliconArchitect → AmpliconClassifier
│   ├── 02_run_cycleviz.sh               # ecDNA cycle visualization (Fig 1B)
│   └── 03_extract_MYCN_CN.py            # MYCN copy number per line from AmpliconClassifier / CNVkit output
├── 03_chromothripsis/
│   ├── 01_delly_call.sh                 # SV calling (DELLY 1.7.2 static binary)
│   ├── 02_filter_delly_vcf.sh           # PASS calls, chr1–22, X (tab-delimited column checks)
│   ├── 03_format_sv_for_shatterseek.py
│   ├── 04_format_cn_for_shatterseek.py  # CNVkit *_CNV_CALLS.bed → ShatterSeek CN tables (also Fig 1B input)
│   ├── 05_run_shatterseek.R             # env ASCAT
│   └── 06_call_high_confidence_chromothripsis.ipynb
├── 04_rnaseq/
│   ├── 01_fastqc.sh
│   ├── 02_trim_galore.sh
│   ├── 03_sortmerna.sh
│   ├── 04_star_2pass.sh
│   ├── 05_picard_markduplicates.sh      # batch 1 only; flags duplicates for QC, nothing is removed
│   ├── 06_stringtie_batch1.sh           # FPKM, reference-guided assembly mode, on the duplicate-flagged BAM
│   ├── 06_stringtie_batch2.sh           # FPKM, reference-only mode (-e), on the STAR BAM
│   ├── 07_featurecounts.sh              # fragment counts (-p --countReadPairs), duplicates retained
│   ├── 08_merge_counts.sh               # gene × sample matrix, joined on Ensembl gene ID (fails on ID mismatch)
│   └── 08b_featurecounts_sanity_check.sh
├── 05_differential_expression/
│   ├── 01_deseq2_ecDNA.R                # ecDNA+ vs ecDNA− (groups from samples.tsv; reference level ecDNA−)
│   ├── 02_gsea_and_RS_genes.R           # Hallmark GSEA (fgsea, seed 42, MSigDB 2026.1.Hs), 29 RS genes
│   └── 03_export_FPKM_tables.R          # MYCN / CHEK1 FPKM from StringTie gene abundance tables; RS heatmap
└── 06_figures/
    ├── 00_build_figure_tables.py        # assembles data/MYCN_FPKM.txt, FPKM_CN_MYCN.txt, master.csv
    ├── fig1A_oncoprint.ipynb
    ├── fig1B_chr2_copy_number.py
    ├── fig2A_MYCN_CN_vs_FPKM.py
    ├── fig2B_MYCN_FPKM_boxplot.py
    └── fig2C_volcano.py
```

## Usage

1. `cp config.example.sh config.sh` and fill in the paths, scheduler account and conda base. Every shell script sources
   `config.sh`; none contains an absolute path.
2. Run the folders in numeric order. Job-generator scripts write one SLURM/PBS script per sample from `samples.tsv`
   (or the per-batch `$RNA_LIST`) and are submitted with `sbatch`/`qsub`; every script checks its inputs, logs the tool
   version and stops on the first error. RNA-seq steps run once per sequencing batch; `08_merge_counts.sh --join`
   then combines the two batch matrices into `counts_21lines.txt` for DESeq2.
3. Rebuild the figure tables after any upstream change:

```bash
python 02_ecDNA/03_extract_MYCN_CN.py --dir <AmpliconSuite output> --samples samples.tsv --id-column dna_id --out data/MYCN_CN.tsv
Rscript 05_differential_expression/03_export_FPKM_tables.R          # -> data/MYCN_CHEK1_FPKM.csv, CHEK1_FPKM.csv
python 06_figures/00_build_figure_tables.py --samples samples.tsv --fpkm data/MYCN_CHEK1_FPKM.csv \
       --cn data/MYCN_CN.tsv --master-in data/master.csv --out data/
```

Notes
- **STAR 2-pass:** run samples one at a time; each sample rebuilds a shared splice-junction index.
- **Two RNA-seq batches.** Batch 1 (10 lines) was quantified with StringTie in reference-guided assembly mode on the
  Picard-flagged BAM; batch 2 (11 lines) in reference-only mode (`-e`) on the STAR BAM. StringTie ignores the duplicate flag
  (verified: identical gene FPKM for CHP212 with and without flagging), and featureCounts was run without duplicate
  filtering, so duplicates are retained in all quantifications. Both batches use the same STAR index and GENCODE v44 GTF.
- **Sample IDs.** Three lines carry different IDs in the WGS data (`dna_id` column of `samples.tsv`: NB10 = NBro,
  SKNBE2 = SKNEB2, SKNrI = SKNYrI); `wgs_id` holds the display names used in `master.csv` (e.g. SKN-BE2, SMS-KCNR, SJ-NB10).
- **MYCN copy number** (Figs 1A/2A) is the AmpliconClassifier `gene_cn` of MYCN (AmpliconArchitect-predicted copy count of
  the MYCN-containing segment); for the six non-amplified lines it is the CNVkit segment overlapping MYCN.
- **Superseded AmpliconSuite run.** An earlier run (PrepareAA from FASTQ, CNVkit `hmm-tumor` segmentation) was replaced by
  the run documented here (from the duplicate-marked BAM, CBS segmentation). The hand-entered copy numbers in the original
  `master.csv` came from the earlier run and were replaced by scripted values from the final run (`REANALYSIS.md`).
- `master.csv` columns `ecDNA_counts` and `oncogene_on_ecDNA` are curated from AmpliconClassifier output and are checked by
  `fig1A_oncoprint.ipynb` against the AmpliconArchitect summaries in `data/ecDNA_summary/` (the notebook stops on any mismatch);
  `Chromothripsis_*` and `SV_burden` come from `03_chromothripsis/`; the MYCN/CHEK1 columns are rebuilt by the scripts above.

## Software versions

Run-time versions; the evidence for each (conda export, run record or log) is listed in `envs/versions.md`.

| Tool | Version |
|---|---|
| Reference genome (WGS) | GRCh38.d1.vd1 |
| Reference genome / annotation (RNA-seq) | GENCODE GRCh38 primary assembly, GENCODE v44 basic |
| BWA-MEM | 0.7.17-r1188 |
| SAMtools | 1.9 (WGS alignment), 1.6 (AmpliconSuite env) |
| Picard MarkDuplicates | 2.18.27 (WGS), 3.0.0 (RNA-seq batch 1, flag only) |
| AmpliconSuite-pipeline (PrepareAA) | 1.2.2 |
| AmpliconArchitect | 1.3.r7 |
| AmpliconClassifier | 1.1.1 |
| CNVkit | 0.9.10 (CBS segmentation) |
| CycleViz | 0.2.2 |
| DELLY | 1.7.2 |
| ShatterSeek | 1.1 (R 4.2.2) |
| FastQC | 0.12.1 |
| Trim Galore | 0.6.10 |
| SortMeRNA | 4.3.6 |
| STAR | 2.7.10b |
| StringTie | 2.2.1 |
| Subread (featureCounts) | 2.0.6 |
| R | 4.2.2 |
| DESeq2 | 1.36.0 |
| fgsea | 1.22.0 |
| msigdbr (MSigDB) | 26.1.0 (v2026.1.Hs) |
| org.Hs.eg.db | 3.15.0 |

## Reanalysis (2026 code audit)

The RNA-seq count matrix was rebuilt with fragment counting (`featureCounts -p --countReadPairs`) and an ID-keyed merge,
and DESeq2/GSEA were re-run; the ecDNA+ vs ecDNA− signature changed from 498/374 to 486/404 genes (FDR < 0.05,
|log2FC| > 1) with the same Hallmark enrichments. The MYCN FPKM and copy-number tables were rebuilt from pipeline output.
Original outputs are preserved; all changes and the reasons for them are recorded in `REANALYSIS.md`.

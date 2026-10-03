# ecDNA-MYCN-neuroblastoma

Analysis and figure code for:

> Yang T, *et al.* **Extrachromosomal *MYCN* Amplification Defines a Replication Stress-Addicted Neuroblastoma Subset with Selective Vulnerability to CHK1 Inhibition.** *Journal*, Year. DOI: `TODO`

This repository covers the whole-genome sequencing (WGS) and RNA-seq analyses of 21 neuroblastoma cell lines behind **Figures 1 and 2**:
- ecDNA detection
- chromothripsis calling
- differential expression and gene set enrichment
- figure generation

## Figure-to-code map

| Panel | Script |
|---|---|
| 1A | `06_figures/fig1A_oncoprint.ipynb` |
| 1B (ecDNA cycle, NGP) | `02_ecDNA/02_run_cycleviz.sh` |
| 1B (chr2 copy number) | `06_figures/fig1B_chr2_copy_number.py` |
| 2A | `06_figures/fig2A_MYCN_CN_vs_FPKM.py` |
| 2B | `06_figures/fig2B_MYCN_FPKM_boxplot.py` |
| 2C | `05_differential_expression/01_deseq2_ecDNA.R` → `06_figures/fig2C_volcano.py` |
| 2D, 2E | `05_differential_expression/02_gsea_and_RS_genes.R` |

## Repository structure

```
├── 01_wgs_alignment/
│   └── bwa_mem_markduplicates.pbs       # BWA-MEM → samtools sort → Picard MarkDuplicates
├── 02_ecDNA/
│   ├── 01_run_ampliconsuite.sh          # CNVkit → AmpliconArchitect → AmpliconClassifier
│   └── 02_run_cycleviz.sh               # ecDNA cycle visualization
├── 03_chromothripsis/
│   ├── 01_delly_call.sh                 # SV calling
│   ├── 02_filter_delly_vcf.sh           # PASS calls, chr1–22, X
│   ├── 03_format_sv_for_shatterseek.py
│   ├── 04_format_cn_for_shatterseek.py
│   ├── 05_run_shatterseek.R
│   └── 06_call_high_confidence_chromothripsis.ipynb
├── 04_rnaseq/
│   ├── 01_fastqc.sh
│   ├── 02_trim_galore.sh
│   ├── 03_sortmerna.sh
│   ├── 04_star_2pass.sh
│   ├── 05_picard_markduplicates.sh
│   ├── 06_stringtie.sh                  # FPKM
│   ├── 07_featurecounts.sh              # read counts
│   ├── 08_merge_counts.sh               # gene × sample count matrix
│   └── 08b_featurecounts_sanity_check.sh
├── 05_differential_expression/
│   ├── 01_deseq2_ecDNA.R                # ecDNA+ vs ecDNA−
│   ├── 02_gsea_and_RS_genes.R           # Hallmark GSEA, replication stress genes
│   └── 03_export_CHEK1_FPKM.R
└── 06_figures/
```

## Usage

Run the folders in numeric order.

- **Paths:** each script starts with a `Paths: edit before running` block, and scheduler headers use `<partition>`/`<account>` placeholders. Commands and parameters are as used in the paper.
- **Job generators:** most cluster scripts write one job script per sample, which is then submitted with `sbatch` or `qsub`.
- **Single-sample scripts:** these are shown for CHP134 and were run identically for every cell line.
- **STAR 2-pass:** run samples one at a time. Each sample rebuilds a shared splice-junction index.
- **Manual tables:**
  - `master.csv` (Fig 1A), `FPKM_CN_MYCN.txt` (Fig 2A) and `MYCN_FPKM.txt` (Fig 2B) were assembled manually from the AmpliconClassifier, ShatterSeek and StringTie outputs.
  - The 21-sample count matrix was assembled by joining per-sample featureCounts columns on Ensembl gene ID.

## Software versions

| Tool | Version |
|---|---|
| Reference genome (WGS) | GRCh38.d1.vd1 |
| Reference genome / annotation (RNA-seq) | GENCODE GRCh38 primary assembly, GENCODE v44 basic |
| BWA-MEM | ____ |
| SAMtools | ____ |
| Picard | 2.18.27 |
| AmpliconSuite-pipeline | 1.2.2 |
| CNVkit | 0.9.12 |
| AmpliconArchitect | ____ |
| AmpliconClassifier | 1.1.1 |
| CycleViz | 0.2.2 |
| DELLY | 1.7.2 |
| ShatterSeek | 1.1 |
| FastQC | 0.12.1 |
| Trim Galore | 0.6.10 |
| SortMeRNA | 4.3.6 |
| STAR | 2.7.10b |
| StringTie | 2.2.1 |
| Subread (featureCounts) | 2.0.6 |
| R | ____ |
| DESeq2 | 1.36.0 |
| fgsea | 1.22.0 |
| msigdbr (MSigDB) | ____ (v2026.1.Hs) |
| Python | ____ |
| pandas / numpy / matplotlib / scipy / adjustText | ____ |

## Data availability

- **Raw sequencing data:** `TODO` (accession, or as stated in the manuscript).
- **Intermediate files:** available from the corresponding authors upon reasonable request.

## Contact

Ting Yang: `TODO email`; Ludmil B. Alexandrov: L2alexandrov@health.ucsd.edu

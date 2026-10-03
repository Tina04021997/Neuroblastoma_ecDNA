# Extrachromosomal *MYCN* amplification defines a replication stress–addicted neuroblastoma subset

Code accompanying:

> Yang T, *et al.* **Extrachromosomal *MYCN* Amplification Defines a Replication Stress-Addicted Neuroblastoma Subset with Selective Vulnerability to CHK1 Inhibition.** *Journal*, Year. DOI: `TODO`

[![DOI](https://zenodo.org/badge/DOI/TODO.svg)](https://doi.org/TODO) <!-- archive a release on Zenodo before submission -->

This repository contains the bioinformatic analyses and plotting code for **Figures 1 and 2**: WGS and RNA-seq processing of 21 neuroblastoma cell lines, ecDNA detection and reconstruction, chromothripsis calling, differential expression, gene set enrichment, and figure generation.

Figures 3 and 4 were generated outside this repository. Figure 3 uses the R2 and DepMap web portals. Figure 4 contains wet-lab data.

---

## Figure-to-code map

| Panel | What it shows | Script | Input (produced upstream) |
|---|---|---|---|
| **1A** | Oncoprint: MYCN status, ecDNA counts, chromothripsis, MYCN CN/FPKM, oncogenes on ecDNA, CHEK1 FPKM | `06_figures/fig1/fig1A_oncoprint.ipynb` | `master.csv` (manually assembled, see `data/README.md`) + AA `<s>_summary.txt` (step 02) |
| **1B** circle | ecDNA cycle reconstruction, NGP | `02_ecDNA_AmpliconSuite/02_run_cycleviz.sh` | AA `NGP_amplicon1_cycles.txt`, `_graph.txt` |
| **1B** schematic | Linear chr2 with *MYCN* | not code: `TODO` (e.g. BioRender) | — |
| **1B** right | Chr2 copy number, NGP and SKNFI | `06_figures/fig1/fig1B_chr2_copy_number.py` | `<s>_cn.txt` CN segments (step 03, `04_format_cn_for_shatterseek.py`) |
| **2A** | *MYCN* CN vs FPKM | `06_figures/fig2/fig2A_MYCN_CN_vs_FPKM.py` | `FPKM_CN_MYCN.txt` (manually assembled from steps 02 + 04) |
| **2B** | *MYCN* FPKM, ecDNA+ vs ecDNA− | `06_figures/fig2/fig2B_MYCN_FPKM_boxplot.py` | `MYCN_FPKM.txt` (manually assembled from step 04) |
| **2C** | Volcano plot | `05_differential_expression/01_deseq2_ecDNA.R` → `06_figures/fig2/fig2C_volcano.py` | count matrix (step 04) → `deseq2_results_ecDNA.csv` (step 05) |
| **2D** | Hallmark GSEA | `05_differential_expression/02_gsea_and_RS_genes.R` | `deseq2_results_ecDNA.csv` (step 05) |
| **2E** | Replication stress genes | `05_differential_expression/02_gsea_and_RS_genes.R` | `deseq2_results_ecDNA.csv` (step 05) |

---

## Repository structure

```
ecDNA-MYCN-neuroblastoma/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── config/
│   ├── README.md
│   └── samples.tsv                         # 21 cell lines: IDs, display names, groups, batches
├── data/
│   └── README.md                           # data availability + provenance of manually assembled tables
├── 01_wgs_preprocessing/
│   ├── README.md
│   └── 01_bwa_mem_markduplicates.pbs       # BWA-MEM → samtools sort → Picard MarkDuplicates
├── 02_ecDNA_AmpliconSuite/
│   ├── README.md
│   ├── 01a_run_ampliconsuite_pilot.sh      # CNVkit → AA → AC (pilot batch)
│   ├── 01b_run_ampliconsuite_second.sh     # same, second batch (job generator)
│   └── 02_run_cycleviz.sh                  # Fig 1B circle
├── 03_chromothripsis/
│   ├── README.md
│   ├── 01_delly_call.sh
│   ├── 02_filter_delly_vcf.sh
│   ├── 03_format_sv_for_shatterseek.py
│   ├── 04_format_cn_for_shatterseek.py
│   ├── 05_run_shatterseek.R
│   └── 06_call_high_confidence_chromothripsis.ipynb
├── 04_rnaseq/
│   ├── README.md
│   ├── pilot_PBS/                          # 01_fastqc … 07_featurecounts (PBS)
│   └── second_SLURM/                       # 01_fastqc … 08b_sanity_check (SLURM)
├── 05_differential_expression/
│   ├── README.md
│   ├── 01_deseq2_ecDNA.R
│   ├── 02_gsea_and_RS_genes.R
│   └── 03_export_CHEK1_FPKM.R
└── 06_figures/
    ├── README.md
    ├── fig1/  (fig1A_oncoprint.ipynb, fig1B_chr2_copy_number.py)
    └── fig2/  (fig2A_MYCN_CN_vs_FPKM.py, fig2B_MYCN_FPKM_boxplot.py, fig2C_volcano.py)
```

Each numbered folder has its own `README.md` with inputs, outputs, parameters and run order.

---

## Analysis overview

```
WGS FASTQ ─BWA-MEM─► BAM ─Picard─► dedup BAM ─┬─► AmpliconSuite (CNVkit/AA/AC) ─► ecDNA calls, CN segments ──┐
                                              │                                  └─► CycleViz (Fig 1B)        │
                                              └─► DELLY SVs ─┐                                                │
                                         CNVkit CN segments ─┴─► ShatterSeek ─► chromothripsis calls ─────────┤
                                                                                                              ├─► master.csv ─► Fig 1A
RNA FASTQ ─FastQC─TrimGalore─SortMeRNA─STAR 2-pass─┬─► Picard ─► StringTie (FPKM) ────────────────────────────┤
                                                   │                                                          └─► Fig 2A, 2B
                                                   └─► featureCounts ─► DESeq2 ─► Fig 2C ─► fgsea ─► Fig 2D, 2E
```

---

## Software versions

| Step | Tool / resource | Version | Used in | How to check |
|---|---|---|---|---|
| **References** | WGS genome | GRCh38.d1.vd1 (GDC) | 01–03 | — |
| | RNA-seq genome | GENCODE GRCh38 primary assembly, release 44 | 04 | — |
| | Gene annotation | GENCODE v44 basic GTF | 04 | — |
| | AmpliconSuite data repo | GRCh38, `____` | 02 | `$AA_DATA_REPO/GRCh38/` file dates |
| **WGS alignment** | BWA-MEM | `____` | 01 | `bwa 2>&1 \| grep Version` |
| | SAMtools | `____` | 01 | `samtools --version` |
| | Picard | `____` | 01, 04 | `picard MarkDuplicates --version` |
| **ecDNA** | AmpliconSuite-pipeline | 1.2.2 | 02 | `PrepareAA.py --version` |
| | CNVkit | 0.9.12 | 02 | `cnvkit.py version` |
| | AmpliconArchitect | `____` | 02 | first lines of the AA log |
| | AmpliconClassifier | 1.1.1 | 02 | `amplicon_classifier.py --version` |
| | CycleViz | 0.2.2 | 02 | `CycleViz.py -v` |
| **Chromothripsis** | DELLY | 1.7.2 | 03 | `delly` (prints version) |
| | ShatterSeek | 1.1 | 03 | `packageVersion("ShatterSeek")` |
| | GenomicRanges | `____` | 03 | `packageVersion("GenomicRanges")` |
| **RNA-seq** | FastQC | 0.12.1 | 04 | `fastqc --version` |
| | Trim Galore | 0.6.10 | 04 | `trim_galore --version` |
| | cutadapt | `____` | 04 | `cutadapt --version` |
| | SortMeRNA | 4.3.6 (db `smr_v4.3_default_db`) | 04 | `sortmerna --version` |
| | STAR | 2.7.10b | 04 | `STAR --version` |
| | StringTie | 2.2.1 | 04 | `stringtie --version` |
| | Subread (featureCounts) | 2.0.6 | 04 | `featureCounts -v` |
| **R** | R | `____` | 03, 05 | `R --version` |
| | DESeq2 | 1.36.0 | 05 | `packageVersion("DESeq2")` |
| | org.Hs.eg.db | `____` | 05 | `packageVersion("org.Hs.eg.db")` |
| | fgsea | 1.22.0 | 05 | `packageVersion("fgsea")` |
| | msigdbr / MSigDB | `____` / v2026.1.Hs | 05 | `packageVersion("msigdbr")` |
| | tidyverse / ggplot2 | `____` / `____` | 05 | `packageVersion("tidyverse")`, `packageVersion("ggplot2")` |
| | pheatmap, RColorBrewer | `____` | 05 | `packageVersion("pheatmap")` |
| | gridExtra, cowplot | `____` | 03 | `packageVersion("cowplot")` |
| **Python** | Python | `____` | 03, 06 | `python --version` |
| | pandas / numpy | `____` / `____` | 03, 06 | `python -c "import pandas, numpy; print(pandas.__version__, numpy.__version__)"` |
| | matplotlib / scipy | `____` / `____` | 06 | `python -c "import matplotlib, scipy; print(matplotlib.__version__, scipy.__version__)"` |
| | adjustText | `____` | 06 | `pip show adjustText` |
| | Jupyter | `____` | 03, 06 | `jupyter --version` |

All versions should match the Methods section.

---

## Data availability

Raw sequencing data are **not** stored in this repository.

| Data | Location |
|---|---|
| WGS FASTQ (21 cell lines) | `TODO: SRA BioProject PRJNAxxxxxx` |
| RNA-seq FASTQ (21 cell lines) | `TODO: GEO GSExxxxxx / SRA` |
| Intermediate/processed tables (count matrix, CN segments, summary tables) | Available from the corresponding authors upon reasonable request |

---

## Running the analysis

Run the numbered folders in order, following each folder's `README.md`.

- **Steps 01–04** (raw reads → BAMs, ecDNA calls, chromothripsis calls, expression tables) run on an HPC cluster.
- **Steps 05–06** (statistics and figures) run on a laptop from the tables produced upstream.
- **Schedulers:** the scripts were written for the UCSD TSCC cluster. Pilot RNA-seq steps use PBS/Torque; WGS and second-batch steps use SLURM.
- **Job generators:** most `*.sh` files loop over a sample list and write one job script per sample, which is then submitted with `qsub`/`sbatch`.
- **Single-sample files** (e.g. `01_bwa_mem_markduplicates.pbs`, `01a_run_ampliconsuite_pilot.sh`) are the CHP134 instance of a template run identically for every cell line.
- **Batches:** samples were sequenced in two batches (`Pilot`, `Second`; see `config/samples.tsv`) and processed with the same tools and parameters. The one exception is a STAR flag, noted in `04_rnaseq/README.md`.
- **Paths:** every script starts with a `Paths: edit before running` block (`/path/to/...` placeholders), and scheduler headers use `<queue>`, `<partition>` and `<account>`. Edit these for your system; commands and parameters below the block are exactly as run for the paper.

| Step | Folder | Resources per sample |
|---|---|---|
| WGS alignment + MarkDuplicates | `01_wgs_preprocessing` | 28 cores |
| AmpliconSuite | `02_ecDNA_AmpliconSuite` | 32–64 cores, up to 1 TB RAM, ≤250 h |
| DELLY | `03_chromothripsis` | 4 cores, 60 GB, ≤20 h |
| STAR 2-pass | `04_rnaseq` | 32–64 cores |

---

## Sample groups

There are 21 neuroblastoma cell lines: 15 *MYCN*-amplified/ecDNA+ and 6 *MYCN* non-amplified/ecDNA−.

In this cohort, *MYCN* amplification and ecDNA status are identical. The DESeq2 design variable `MYCN_amp` therefore equals the ecDNA+/ecDNA− comparison reported in the paper.

File-level IDs (e.g. `KAN`, `SKNBE2`, `LANS`) map to manuscript names (SMS-KAN, SK-N-BE(2), LAN-5) in `config/samples.tsv`.

---

## Citation

Please cite the paper above and the archived code release (`CITATION.cff`).

## License

`TODO` (e.g. MIT)

## Contact

Ting Yang: `TODO email`; Ludmil B. Alexandrov: L2alexandrov@health.ucsd.edu

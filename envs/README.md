# Software environments

Exports of the conda environments used on TSCC (`conda env export --no-builds`; `ampsuite` via mamba; all 2026-10-26) and the
`sessionInfo()` of the local R installation used for DESeq2 / fgsea / figure tables.
Recreate an environment with `conda env create -f envs/<name>.yml`.

| Environment / file | Used by | Key tools |
|---|---|---|
| `evc_main.yml` | `01_wgs_alignment/*` (WGS alignment, duplicate marking) | bwa 0.7.17, Picard 2.18.27, GATK 4.1.4.1, bcftools 1.9, bedtools 2.30.0, mosdepth 0.3.3, python 3.7; SAMtools 1.9 was in this env at run time and was later moved to its own env (not in the export) |
| `ampsuite.yml` | `02_ecDNA/01_run_ampliconsuite.sh` (AmpliconSuite-pipeline / AmpliconArchitect / AmpliconClassifier, CNVkit) | run time (2024-02): AmpliconSuite-pipeline 1.2.2, AmpliconArchitect 1.3.r7, AmpliconClassifier 1.1.1, CNVkit 0.9.10; export today: ampliconsuite 1.4.0, cnvkit 0.9.12, samtools 1.6, bwa 0.7.17, python 3.10.13 (see note) |
| `ASCAT.yml` | `03_chromothripsis/*` (ShatterSeek, run in this env) | R 4.2.2, ShatterSeek 1.1 (GitHub package, not in the yml), cancerit-allelecount 4.0.2, htslib 1.9 |
| `rna-seq.yml` | `04_rnaseq/*` (Trim Galore 0.6.10, FastQC 0.12.1, STAR 2.7.10b, Picard 3.0.0, StringTie 2.2.1, featureCounts / Subread 2.0.6) | python 3.10, R 4.3.1 (not used for DESeq2) |
| `R_sessionInfo_local.txt` | `05_differential_expression/*` (DESeq2, GSEA, FPKM tables / RS heatmap) | R 4.2.2, DESeq2 1.36.0, fgsea 1.22.0, msigdbr 26.1.0 (MSigDB 2026.1.Hs), org.Hs.eg.db 3.15.0, dplyr 1.2.0, ggplot2 4.0.2 |
| `figures_python.txt` | `06_figures/*.py`, `fig1A_oncoprint.ipynb`, `02_ecDNA/*.py`, `03_chromothripsis/*.py` (local Python steps, `make_figures.sh`) | Python 3.9.13, pandas 2.2.3, numpy 1.24.4, matplotlib 3.9.4, scipy 1.13.1, adjustText 0.8, nbconvert 6.4.4 |

Notes
- DELLY: the structural-variant calls (`03_chromothripsis/01_delly_call.sh`, VCF headers `SVMETHOD=EMBL.DELLYv1.7.2`, run 2026-02-13)
  were produced with the DELLY v1.7.2 static binary (`delly_v1.7.2_linux_x86_64bit`, downloaded from
  https://github.com/dellytools/delly/releases, kept under `$SOFTWARE_DIR`; set `DELLY_BIN` in `config.sh`).
  It is not a conda package; the `ASCAT.yml` environment contains an older `delly=0.7.6`, which was not used.
- `ampsuite.yml` was exported with mamba from `miniforge3/envs/ampsuite` on 2026-10-26 and reflects the environment *today*
  (bioconda `ampliconsuite=1.4.0`, `cnvkit=0.9.12`); it was updated after the analysis. The versions actually used are recorded
  by the pipeline in `<sample>_run_metadata.json` (batch 2 launched 2024-02-09, batch 1 re-run 2024-02-19; excerpts in `run_records/`):
  PrepareAA.py from the git clone `~/AmpliconSuite-pipeline` (PAA_version **1.2.2**), AmpliconArchitect from the env
  `site-packages/ampliconarchitectlib` (**AmpliconArchitect 1.3.r7**, python 3.10.13, from `<sample>_AA_results/<sample>.log`), AmpliconClassifier
  **1.1.1** (`envs/ampsuite/bin/amplicon_classifier.py`), CNVkit **0.9.10** (`-m cbs`, reference
  `GRCh38_cnvkit_filtered_ref.cnn` from the AA data repo). These run-time versions, not the export, are what the README and
  manuscript cite (see REANALYSIS.md).
- Python figure scripts (`06_figures/*.py`): package versions are not pinned inside the scripts; the versions used for the
  2026-10-27 regeneration are listed in `figures_python.txt`. `fig2C_volcano.py` prints its matplotlib / pandas versions at run time.

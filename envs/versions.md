# Software versions used for the published analysis (audit point 6)

Run-time versions, with the evidence file for each. Where the current conda export differs from the
run-time version, the environment was updated after the analysis and the run-time value is the one cited
in the Methods and README.

| Stage (scripts) | Tool | Run-time version | Evidence | Manuscript |
|---|---|---|---|---|
| WGS alignment (`01_wgs_alignment/`) | bwa | 0.7.17 | `envs/evc_main.yml` | not stated |
| | Picard MarkDuplicates | 2.18.27 | `envs/evc_main.yml` | 2.18.27 |
| | SAMtools | 1.9 | in env `evc_main` at run time (`01_wgs_alignment/bwa_mem_markduplicates.pbs` header); later moved to a separate env, so absent from `evc_main.yml` | not stated |
| | GATK | 4.1.4.1 | `envs/evc_main.yml` | not stated |
| ecDNA (`02_ecDNA/01_run_ampliconsuite.sh`) | AmpliconSuite-pipeline (PrepareAA.py, git clone) | 1.2.2 | `<sample>_run_metadata.json`, `PAA_version` | 1.2.2 |
| | AmpliconArchitect (`ampliconarchitectlib`) | 1.3.r7 | `<sample>_AA_results/<sample>.log` | not stated |
| | AmpliconClassifier | 1.1.1 | `<sample>_run_metadata.json`, `AC_version` | 1.1.1 |
| | CNVkit (`-m cbs`, `GRCh38_cnvkit_filtered_ref.cnn`) | **0.9.10** | `<sample>_run_metadata.json`, `cnvkit_version` | 0.9.12 in v6–v9, corrected to 0.9.10 in v10 |
| | samtools (AmpliconSuite env) | 1.6 | `envs/ampsuite.yml` | not stated |
| | Python (AmpliconSuite env) | 3.10.13 | AA log | not stated |
| | CycleViz (`02_ecDNA/02_run_cycleviz.sh`, git clone under `$SOFTWARE_DIR`) | 0.2.2 | `CycleViz/_version.py`, checked 2026-10-26 | not stated |
| SV / chromothripsis (`03_chromothripsis/`) | DELLY (static binary `delly_v1.7.2_linux_x86_64bit`) | 1.7.2 | VCF header `SVMETHOD=EMBL.DELLYv1.7.2`, `fileDate=20260213` | 1.7.2 |
| | ShatterSeek (ASCAT env, R 4.2.2) | 1.1 | `packageVersion("ShatterSeek")` in env ASCAT, 2026-10-26 | 1.1 |
| | R (ASCAT env) | 4.2.2 | `envs/ASCAT.yml` | — |
| RNA-seq (`04_rnaseq/`) | Trim Galore | 0.6.10 | `envs/rna-seq.yml` | 0.6.10 |
| | FastQC | 0.12.1 | `envs/rna-seq.yml` | 0.12.1 |
| | STAR | 2.7.10b | `envs/rna-seq.yml` | 2.7.10b |
| | Picard MarkDuplicates (batch 1, flag only) | 3.0.0 | `envs/rna-seq.yml` | 3.0.0 |
| | StringTie | 2.2.1 | `envs/rna-seq.yml` | 2.2.1 |
| | featureCounts (Subread) | 2.0.6 | `envs/rna-seq.yml` | 2.0.6 |
| | GENCODE annotation | v44 basic | `config.sh` (`GTF`) | v44 |
| Differential expression / GSEA (`05_differential_expression/`) | R | 4.2.2 | `envs/R_sessionInfo_local.txt` | 4.2.2 |
| | DESeq2 | 1.36.0 | `envs/R_sessionInfo_local.txt` | 1.36.0 |
| | fgsea | 1.22.0 | `envs/R_sessionInfo_local.txt` | 1.22.0 |
| | msigdbr / MSigDB | 26.1.0 / 2026.1.Hs | `envs/R_sessionInfo_local.txt` (`02_gsea_and_RS_genes.R` also writes `gsea_session_versions.txt`) | 26.1.0 / 2026.1.Hs |
| | org.Hs.eg.db | 3.15.0 | `envs/R_sessionInfo_local.txt` | stated |
| Figures (`06_figures/*.py`) | Python / matplotlib / pandas | not pinned (patched `fig2C_volcano.py` prints versions at run time) | — | — |

Current-export versions that differ from run time (`envs/ampsuite.yml`, exported 2026-10-26): bioconda `ampliconsuite` 1.4.0,
CNVkit 0.9.12. The git clone `~/AmpliconSuite-pipeline` is at `v1.4.0-2-gb409006` (Sep 2025); the runs were launched
2024-02-09 (batch 2, KAN) and 2024-02-19 (batch 1, CHP134) with the clone at 1.2.2; both batches report CNVkit 0.9.10,
AmpliconSuite-pipeline 1.2.2 and AmpliconClassifier 1.1.1 (`run_records/`).

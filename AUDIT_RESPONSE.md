# Response to the code audit (Neuroblastoma_Code_Fixes.docx)

Repository `Neuroblastoma_ecDNA-main`, manuscript v11, figure file `Figure_rerun.pptx/.pdf`. Status as of 2026-10-27.

Audit principle followed throughout: results were not altered to match manuscript values; original 2024 outputs are preserved
in `data/original_2024_readcounts/`, every re-analysis is documented in `REANALYSIS.md`, and the one change that altered
numbers (point 4) is reported with its before/after comparison below.

## Point-by-point status

| # | Audit item | Status | What was done | Where |
|---|---|---|---|---|
| 1 | Correct DELLY VCF filtering | Closed | awk '\s' regex replaced by explicit tab-delimited CHROM / FILTER column tests; stops if a sample keeps no record. Applied to the 21 DELLY VCFs, the fixed filter reproduces byte for byte the 21 SV files used for ShatterSeek, so variant counts, chromosome distribution and the chromothripsis calls (17 lines, 30 events) are unchanged. | 03_chromothripsis/02_filter_delly_vcf.sh; REANALYSIS.md 'Audit points 1-4' |
| 2 | Reconcile DESeq2 count-matrix input | Closed | DESeq2 reads the version-stripped counts_21lines.txt from 08_merge_counts.sh --join (not merged_counts.txt), strips residual version suffixes, maps symbols explicitly and stops unless > 50 % map (35,377 / 62,700; 99 duplicate symbols dropped). Re-run on the original read-level matrix reproduces 498 / 374 exactly; fragment-level matrix gives 486 / 404. | 05_differential_expression/01_deseq2_ecDNA.R; data/deseq2_gsea/, data/original_2024_readcounts/ |
| 3 | Make count merging gene-ID safe | Closed | 08_merge_counts.sh runs 08b automatically; every sample's gene-ID column must equal the reference (content and order) or the merge aborts; --join requires identical ID columns across batches. Demonstrated with a shuffled and a truncated synthetic input (exit 1, nothing merged). | 04_rnaseq/08_merge_counts.sh, 08b_featurecounts_sanity_check.sh |
| 4 | Verify paired-end featureCounts behavior | Closed - results changed | Subread 2.0.6 with -p alone counted reads (Assigned / aligned pairs = 1.46-1.67). Re-counted with -p --countReadPairs -B -C; DE signature 498 / 374 -> 486 / 404 (log2FC r = 0.958), Hallmark 18 -> 22 significant sets (all 18 retained), RS-panel stars unchanged. Manuscript and Figures 2C-E updated; original outputs preserved. | 04_rnaseq/07_featurecounts.sh; REANALYSIS.md (first section and validation table) |
| 5 | Clarify duplicate handling | Closed | BAM provenance per branch documented: StringTie batch 1 on the Picard-flagged BAM (flag ignored; verified identical FPKM for CHP212 with/without flagging), batch 2 on the STAR BAM in -e mode; featureCounts on unflagged BAMs, no duplicate filtering. Methods rewritten to the executed commands. | 04_rnaseq/05_picard_markduplicates.sh, 06_stringtie_batch1.sh, 06_stringtie_batch2.sh; manuscript Methods; REANALYSIS.md 'Audit point 5' |
| 6 | Pin computational versions | Closed | Run records (<sample>_run_metadata.json, AA log) give CNVkit 0.9.10 (CBS), AmpliconSuite 1.2.2, AmpliconArchitect 1.3.r7, AmpliconClassifier 1.1.1; manuscript corrected 0.9.12 -> 0.9.10 (v10). All tool versions tabulated with their evidence; conda exports, R sessionInfo and the Python package set archived. | README.md 'Software versions'; envs/versions.md, envs/*.yml, envs/run_records/, envs/figures_python.txt |
| 7 | Pin MSigDB Hallmark gene sets | Closed | msigdbr(collection = 'H') with stopifnot(packageVersion('msigdbr') >= '26.1.0'); the 50-set membership used (7,331 gene rows) is archived and re-read; NES / FDR reproduced byte for byte. Pinning changed no value. | 05_differential_expression/02_gsea_and_RS_genes.R; data/deseq2_gsea/hallmark_genesets_msigdbr_26.1.0.tsv, gsea_session_versions.txt |
| 8 | Make group definitions explicit | Closed | Groups come from samples.tsv (MYCN_amp and ecDNA_group columns). The DESeq2 script asserts valid labels, exactly 6 no_amp / 15 amp, and that MYCN_amp == ecDNA_group for all 21 lines (they coincide in this cohort); results byte-identical. | samples.tsv; 05_differential_expression/01_deseq2_ecDNA.R |
| 9 | Standardize replication-stress panel size | Closed | Single 29-gene panel defined once, written to rs_genes.txt and re-read by the FPKM script; both scripts stopifnot 29 unique symbols and that all 29 are present in the DESeq2 results. 'n = 30' string and '30 29' text corrected; Methods and Fig 2E legend say 29. | 02_gsea_and_RS_genes.R, 03_export_FPKM_tables.R; data/deseq2_gsea/rs_genes.txt |
| 10 | Clarify Figure 2E FDR encoding | Closed | Dual encoding kept by author decision: colour marks a deliberately lenient exploratory threshold (FDR < 0.2) for direction across the panel, asterisks carry FDR < 0.05 / 0.01 / 0.001. Methods and legend now say so explicitly and define grey as FDR >= 0.2 (manuscript v11); no colour or star changed. | manuscript v11 Methods and Fig 2E legend; 02_gsea_and_RS_genes.R comment; REANALYSIS.md 'Audit point 10' |
| 11 | Provide frozen figure inputs | Closed | data/ holds master.csv, FPKM_CN_MYCN.txt, MYCN_FPKM.txt, MYCN_CN.tsv, the count matrix, DESeq2 / GSEA results, AmpliconArchitect / AmpliconClassifier output, chr2 copy-number tables (heatmap_condition.csv is superseded by samples.tsv). bash 06_figures/make_figures.sh regenerates every Fig 1-2 panel from these files; verified locally, panels match the figure file. | data/; 06_figures/make_figures.sh; REANALYSIS.md 'Audit point 11' |
| 12 | Eliminate manual data handoffs | Closed | master.csv MYCN / CHEK1 / copy-number columns, FPKM_CN_MYCN.txt and MYCN_FPKM.txt are built by scripts from AmpliconClassifier, CNVkit and StringTie output; ecDNA columns are asserted against AmpliconClassifier (21/21); sample naming centralised in samples.tsv (sample_id, wgs_id, dna_id). 01_check_figure_tables.py cross-checks IDs, MYCN status, copy number and FPKM across all figure tables (10/10 pass). | 06_figures/00_build_figure_tables.py, 01_check_figure_tables.py; 02_ecDNA/03_extract_MYCN_CN.py, 04_check_ecDNA_columns.py |
| 13 | Add executable QC and failure checks | Closed | All scripts stop on missing/empty inputs (set -euo pipefail, require_file), missing sample lists, missing groups, duplicate gene IDs, failed joins, missing figure inputs and unedited config; logs kept under figures/logs/. Non-zero exits demonstrated for corrupted master.csv, wrong group labels, shuffled counts, missing lists, empty ShatterSeek folder, bad SV records. | config.example.sh helpers; REANALYSIS.md 'Code review before release' and 'Audit points 1-4' |
| 14 | Add Figure 3-4 analysis code | Open (out of scope of this response) | Functional-assay analyses were not performed by the maintainer of this repository; no scripts are included. Stated in README and REANALYSIS. | REANALYSIS.md 'Audit point 14' |
| 15 | Provide a single reproducible workflow | Closed | No placeholder or hard-coded path in any script (cluster scripts read config.sh, which refuses to run unedited; local scripts take arguments and stop on placeholders). README gives a per-stage Workflow table (inputs, commands, environments, outputs), a Data availability section and all versions. Stages 5-6 run cleanly from the repository alone (make_figures.sh); raw data access documented. | README.md 'Workflow', 'Usage', 'Data availability'; envs/ |

## Required final validation: before / after

"Before" = original 2024 analysis (read-level featureCounts, hand-entered tables). "After" = the current repository
(fragment-level counts, scripted tables). Every difference is explained in the Note column and in `REANALYSIS.md`.

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

## Where to look

- `REANALYSIS.md`: every re-run, its reason, its numbers, and the code reviews.
- `README.md`: figure-to-code map, workflow table, usage, versions, data availability.
- `data/`: frozen inputs and results; `envs/`: environments, run records, version evidence.
- `bash 06_figures/make_figures.sh`: regenerates Figures 1-2 from `data/` with all checks.


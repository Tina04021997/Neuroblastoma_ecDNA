#!/bin/bash
# config.sh — the single place for site-specific paths, scheduler settings and helpers.
#
#   cp config.example.sh config.sh      # config.sh is git-ignored; edit the values below
#
# Every script sources this file with
#   source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
# and must not contain any absolute path of its own. Job-generator scripts write the
# resolved values into the per-sample SLURM/PBS scripts they emit.
#
# Values recorded here for the published analysis: TSCC (UCSD), Feb 2024 (RNA-seq and
# WGS processing) and Oct 2026 (featureCounts fragment re-count).

# ── Repository and sample sheet ──────────────────────────────────────────────
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLES_TSV="${REPO_DIR}/samples.tsv"          # 21 lines; column 1 = sample_id (RNA/WGS file prefix)

# ── Scheduler ────────────────────────────────────────────────────────────────
SLURM_PARTITION=<partition>
SLURM_ACCOUNT=<account>
PBS_QUEUE=<queue>
PBS_ACCOUNT=<account>

# ── Software ─────────────────────────────────────────────────────────────────
CONDA_BASE=/path/to/conda                       # e.g. $HOME/miniforge3; envs: evc_main, rna-seq, ampsuite, ASCAT
AMPSUITE_DIR=/path/to/AmpliconSuite-pipeline    # git clone; v1.2.2 at run time (envs/versions.md)
AMPSUITE_ENV=${CONDA_BASE}/envs/ampsuite
CYCLEVIZ_DIR=/path/to/CycleViz                  # v0.2.2 (git clone)
DELLY_BIN=/path/to/delly_v1.7.2_linux_x86_64bit

# ── References ───────────────────────────────────────────────────────────────
REF_DIR=/path/to/references
WGS_REF_FASTA=${REF_DIR}/GRCh38.d1.vd1/GRCh38.d1.vd1.fa
STAR_REF=${REF_DIR}/STAR/GRCh38                 # STAR index built from the two files below
RNA_GENOME_FA=${STAR_REF}/GRCh38.primary_assembly.genome.fa
GTF=${STAR_REF}/gencode.v44.basic.annotation.gtf
RRNA_DB=${REF_DIR}/SortMeRNA/rRNA_databases_v4.3.6
DELLY_EXCLUDE=${REF_DIR}/Delly/hg38.cen

# ── WGS ──────────────────────────────────────────────────────────────────────
WGS_FASTQ_DIR=/path/to/wgs/fastq                # <sample>_1.fastq.gz, <sample>_2.fastq.gz
WGS_ALIGN_DIR=/path/to/wgs/alignment            # <sample>/<sample>_tumor_mkdp.bam
ECDNA_DIR=/path/to/ecDNA_analysis               # AmpliconSuite output, one folder per sample
DELLY_DIR=/path/to/delly                        # <sample>/<sample>.vcf
TMP_DIR=/path/to/tmp

# ── RNA-seq ──────────────────────────────────────────────────────────────────
RNA_DIR=/path/to/rnaseq                         # one folder per batch; both batches use the same layout
RNA_FASTQ_DIR=${RNA_DIR}/fastq                  # <sample>/<sample>_1.fq.gz, <sample>_2.fq.gz
FASTQC_DIR=${RNA_DIR}/fastqc
TRIM_DIR=${RNA_DIR}/trim_galore
SORTRNA_DIR=${RNA_DIR}/sort_rna
STAR_DIR=${RNA_DIR}/STAR
MARKDP_DIR=${RNA_DIR}/markDP
STRINGTIE_DIR=${RNA_DIR}/stringtie
FC_DIR=${RNA_DIR}/featureCounts
RNA_LIST=${RNA_DIR}/list.txt                    # sample IDs of THIS batch, one per line, in matrix column order

# ── Count matrix / DESeq2 (05_differential_expression) ───────────────────────
BATCH1_COUNTS=/path/to/batch1/featureCounts/final_counts.txt   # 08_merge_counts.sh, batch 1 (10 lines)
BATCH2_COUNTS=/path/to/batch2/featureCounts/final_counts.txt   # 08_merge_counts.sh, batch 2 (11 lines)
COUNTS_21=/path/to/deseq2/counts_21lines.txt                  # 08_merge_counts.sh --join -> 01_deseq2_ecDNA.R
DESEQ_DIR=/path/to/deseq2                                     # DESeq2 / GSEA outputs

# ── Helpers ──────────────────────────────────────────────────────────────────
# conda's activation hooks are not `set -u`-clean on every version, so toggle -u around them.
activate_env() {                                 # usage: activate_env rna-seq
    set +u
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "$1"
    set -u
}

require_file() {                                 # usage: require_file path [path ...]
    local f
    for f in "$@"; do
        [ -s "$f" ] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
    done
}

sample_ids() {                                   # usage: for s in $(sample_ids [list.txt]); do ...
    local list="${1:-$RNA_LIST}"                 # default: this batch's list.txt; blank lines skipped
    require_file "$list"
    grep -v '^[[:space:]]*$' "$list"
}

log_version() {                                  # usage: log_version "featureCounts -v" out.txt
    { echo "# $(date -u +%Y-%m-%dT%H:%MZ)  host=$(hostname)"; eval "$1" 2>&1; } > "$2" || true
}

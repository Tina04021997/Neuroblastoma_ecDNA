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
SAMPLES_TSV="${REPO_DIR}/samples.tsv"          # 21 lines; sample_id = RNA file prefix, dna_id = WGS file prefix (3 differ)

# ── Scheduler ────────────────────────────────────────────────────────────────
SLURM_PARTITION="<partition>"
SLURM_ACCOUNT="<account>"
PBS_QUEUE="<queue>"
PBS_ACCOUNT="<account>"

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

wgs_ids() {                                      # usage: ids=$(wgs_ids) || exit 1   -- the dna_id column of samples.tsv
    require_file "$SAMPLES_TSV"                  # (WGS file prefix; differs from sample_id for NB10, SKNBE2, SKNrI)
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="dna_id") c=i; if(!c){print "samples.tsv: no dna_id column" > "/dev/stderr"; exit 1}; next} $c!=""{print $c}' "$SAMPLES_TSV"
}

sample_ids() {                                   # usage: ids=$(sample_ids [list.txt]) || exit 1; for s in $ids; do ...
                                                 # (assign first: a failing $(...) directly in a "for" list does not stop a set -e script)
    local list="${1:-$RNA_LIST}"                 # default: this batch's list.txt; blank lines skipped
    require_file "$list"
    grep -v '^[[:space:]]*$' "$list"
}

log_version() {                                  # usage: log_version "featureCounts -v" out.txt
    { echo "# $(date -u +%Y-%m-%dT%H:%MZ)  host=$(hostname)"; eval "$1" 2>&1; } > "$2" || true
}

# Refuse to run with an unedited template: any value still in <angle brackets> or at /path/to/... stops the sourcing
# script with a clear message (audit point 15). Derived values (${REF_DIR}/..., ${RNA_DIR}/...) are covered by their base.
# Fill in every value listed here even if you only run one part of the pipeline (put a dummy value in the unused ones).
for _v in SLURM_PARTITION SLURM_ACCOUNT PBS_QUEUE PBS_ACCOUNT CONDA_BASE AMPSUITE_DIR CYCLEVIZ_DIR DELLY_BIN REF_DIR \
          WGS_FASTQ_DIR WGS_ALIGN_DIR ECDNA_DIR DELLY_DIR TMP_DIR RNA_DIR BATCH1_COUNTS BATCH2_COUNTS COUNTS_21 DESEQ_DIR; do
    case "${!_v}" in
        *"<"*|/path/to/*) echo "ERROR: ${_v} is still a placeholder in config.sh (edit config.sh first)" >&2; return 1 2>/dev/null || exit 1;;
    esac
done
unset _v

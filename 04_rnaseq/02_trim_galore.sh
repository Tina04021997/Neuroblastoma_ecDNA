#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -t 100:00:00
#SBATCH --mem=100G
#SBATCH --job-name=trim_galore
#SBATCH --output=logs.o/trim.o
#SBATCH --error=logs.e/trim.e
# Partition/account: pass at submission, e.g.  sbatch -p $SLURM_PARTITION -A $SLURM_ACCOUNT 02_trim_galore.sh

# Author: Tina Yang — Trim Galore 0.6.10 (cutadapt) on the raw paired FASTQs of one batch; original run Feb 6, 2024.
# Output: $TRIM_DIR/<sample>_1_val_1.fq.gz, <sample>_2_val_2.fq.gz (+ FastQC of the trimmed reads).

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
activate_env rna-seq

mkdir -p "${TRIM_DIR}/logs.o" "${TRIM_DIR}/logs.e"
cd "${TRIM_DIR}"
log_version "trim_galore --version; cutadapt --version" trim_galore.versions.txt

ids=$(sample_ids) || exit 1          # a failed substitution inside "for" would not stop the script

for file in $ids; do
    require_file "${RNA_FASTQ_DIR}/${file}/${file}_1.fq.gz" "${RNA_FASTQ_DIR}/${file}/${file}_2.fq.gz"
    trim_galore --paired "${RNA_FASTQ_DIR}/${file}/${file}_1.fq.gz" "${RNA_FASTQ_DIR}/${file}/${file}_2.fq.gz" \
        --fastqc --path_to_cutadapt "$(which cutadapt)" --output_dir "${TRIM_DIR}"
done

#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -t 100:00:00
#SBATCH --mem=100G
#SBATCH --job-name=FastQC
#SBATCH --output=logs.o/FASTQC.o
#SBATCH --error=logs.e/FASTQC.e
# Partition/account: pass at submission, e.g.  sbatch -p $SLURM_PARTITION -A $SLURM_ACCOUNT 01_fastqc.sh
# (the SBATCH header cannot read config.sh).

# Author: Tina Yang — FastQC 0.12.1 on the raw paired FASTQs of one sequencing batch; original run Feb 6, 2024.
# Run once per batch from $RNA_DIR (config.sh): samples are read from $RNA_LIST.

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
activate_env rna-seq

mkdir -p "${FASTQC_DIR}/logs.o" "${FASTQC_DIR}/logs.e"
cd "${FASTQC_DIR}"
log_version "fastqc --version" fastqc.versions.txt

ids=$(sample_ids) || exit 1          # a failed substitution inside "for" would not stop the script

for file in $ids; do
    require_file "${RNA_FASTQ_DIR}/${file}/${file}_1.fq.gz" "${RNA_FASTQ_DIR}/${file}/${file}_2.fq.gz"
    fastqc -o "${FASTQC_DIR}" "${RNA_FASTQ_DIR}/${file}/${file}_1.fq.gz" "${RNA_FASTQ_DIR}/${file}/${file}_2.fq.gz"
done

echo "!!!! FastQC done !!!!"

#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -t 100:00:00
#SBATCH --mem=100G
#SBATCH -p <partition>
#SBATCH -A <account>
#SBATCH --job-name=FastQC
#SBATCH --output=logs.o/FASTQC.o
#SBATCH --error=logs.e/FASTQC.e

# Author: Tina Yang
# Date: Feb 6, 2024
# Run FastQC

# Set up envs
mamba activate rna-seq

# ── Paths: edit before running ───────────────────────────────────────────────
WORK_DIR=/path/to/rnaseq                             # contains list.txt (one sample ID per line)
OUT_DIR=${WORK_DIR}/fastqc
FASTQ_DIR=/path/to/rnaseq/fastq                      # ${FASTQ_DIR}/<sample>/<sample>_1.fq.gz, _2.fq.gz
# ─────────────────────────────────────────────────────────────────────────────

## Create output folder
cd $WORK_DIR
mkdir fastqc
cd fastqc
mkdir logs.o
mkdir logs.e

## Run FASTQC
for file in $(cat $WORK_DIR/list.txt); do fastqc -o $OUT_DIR $FASTQ_DIR/${file}/${file}_1.fq.gz $FASTQ_DIR/${file}/${file}_2.fq.gz; done

echo "!!!! FastQC done !!!!"

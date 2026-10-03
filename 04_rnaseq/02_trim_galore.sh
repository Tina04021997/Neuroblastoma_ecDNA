#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -t 100:00:00
#SBATCH --mem=100G
#SBATCH -p <partition>
#SBATCH -A <account>
#SBATCH --job-name=trim_galore
#SBATCH --output=logs.o/trim.o
#SBATCH --error=logs.e/trim.e

# Author: Tina Yang
# Date: Feb 6, 2024
# Run trim_galore

# Set up envs
mamba activate rna-seq

# ── Paths: edit before running ───────────────────────────────────────────────
WORK_DIR=/path/to/rnaseq                             # contains list.txt (one sample ID per line)
OUT_DIR=${WORK_DIR}/trim_galore
FASTQ_DIR=/path/to/rnaseq/fastq                      # ${FASTQ_DIR}/<sample>/<sample>_1.fq.gz, _2.fq.gz
# ─────────────────────────────────────────────────────────────────────────────

## Create output folder
cd $WORK_DIR
mkdir trim_galore
cd trim_galore
mkdir logs.o
mkdir logs.e

## Run Trim Galore
for file in $(cat $WORK_DIR/list.txt); do trim_galore --paired $FASTQ_DIR/${file}/${file}_1.fq.gz $FASTQ_DIR/${file}/${file}_2.fq.gz --fastqc --path_to_cutadapt $(which cutadapt) --output_dir $OUT_DIR; done

#!/bin/bash
# Job generator: writes one SLURM script per sample (<sample>_STAR.sh).
# Run them one at a time: step 3 rebuilds a shared SJ_index directory.

# ── Paths: edit before running ───────────────────────────────────────────────
WORK_DIR=/path/to/rnaseq                             # contains list.txt (one sample ID per line)
STAR_REF=/path/to/references/STAR/GRCh38                    # STAR index; contains GRCh38.primary_assembly.genome.fa + gencode.v44.basic.annotation.gtf
# ─────────────────────────────────────────────────────────────────────────────

cd ${WORK_DIR}
mkdir STAR
cd STAR
mkdir logs.o
mkdir logs.e

for file in $(cat ${WORK_DIR}/list.txt); do echo -e "#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -t 100:00:00
#SBATCH --mem=100G
#SBATCH -p <partition>
#SBATCH -A <account>
#SBATCH --job-name=${file}_STAR
#SBATCH --output=${WORK_DIR}/STAR/logs.o/${file}.o
#SBATCH --error=${WORK_DIR}/STAR/logs.e/${file}.e

# Author: Tina Yang
# Date: Feb 7, 2024
# STAR 2-pass mapping script for PE
# STAR v2.7.10b

# Reference: GRCh38
# Steps:
#1. Indexing genome with annotations
#2. 1-pass mapping with indexed genome
#3. Indexing genome with annotations and SJ.out.tab files
#4. 2-pass mapping with new indexed genome with annotations and SJ.out.tab files

# Set up envs
mamba activate rna-seq

# Dataset
OUT_DIR='${WORK_DIR}/STAR'
FQ_PATH='${WORK_DIR}/sort_rna'
ANNOTATION_PATH='${STAR_REF}'



#  # Step1 Indexing (run once, not per sample)
#  STAR \
#  --runMode genomeGenerate \
#  --genomeDir \$ANNOTATION_PATH \
#  --genomeFastaFiles \$ANNOTATION_PATH/GRCh38.primary_assembly.genome.fa \
#  --sjdbGTFfile \$ANNOTATION_PATH/gencode.v44.basic.annotation.gtf \
#  --runThreadN 32

## Create output folder
mkdir \"\$OUT_DIR/${file}\"
cd \"\$OUT_DIR/${file}\"


# Step2 Mapping-1-pass
STAR \
--genomeDir \"\$ANNOTATION_PATH\" \
--readFilesIn \"\$FQ_PATH/${file}/out/other_fwd.fq\" \"\$FQ_PATH/${file}/out/other_rev.fq\" \
--outFileNamePrefix ${file}_PASS1_result \
--runThreadN 64

# Step3 Indexing genome with annotations and SJ.out.tab file
STAR \
--runMode genomeGenerate \
--genomeDir \"\$ANNOTATION_PATH/SJ_index\" \
--genomeFastaFiles \"\$ANNOTATION_PATH/GRCh38.primary_assembly.genome.fa\" \
--sjdbGTFfile \"\$ANNOTATION_PATH/gencode.v44.basic.annotation.gtf\" \
--runThreadN 64 \
--sjdbFileChrStartEnd \"\$OUT_DIR/${file}/${file}_PASS1_resultSJ.out.tab\"


# Step 4 2-pass mapping output alignments in BAM format and sort BAM by coordinate
STAR \
--genomeDir \"\$ANNOTATION_PATH/SJ_index\" \
--readFilesIn \"\$FQ_PATH/${file}/out/other_fwd.fq\" \"\$FQ_PATH/${file}/out/other_rev.fq\" \
--outFileNamePrefix ${file}_PASS2_result \
--runThreadN 64 \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif
" > "${file}_STAR.sh"; done

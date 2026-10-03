#!/bin/bash
# Job generator: writes one SLURM script per sample (<sample>_sortRNA.sh), then submit with sbatch.

# ── Paths: edit before running ───────────────────────────────────────────────
WORK_DIR=/path/to/rnaseq                             # contains list.txt (one sample ID per line)
RRNA_DATABASE=/path/to/references/SortMeRNA/rRNA_databases_v4.3.6
# ─────────────────────────────────────────────────────────────────────────────

cd ${WORK_DIR}
mkdir sort_rna
cd sort_rna
mkdir logs.o
mkdir logs.e

for file in $(cat ${WORK_DIR}/list.txt); do echo "#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -t 100:00:00
#SBATCH --mem=100G
#SBATCH -p <partition>
#SBATCH -A <account>
#SBATCH --job-name=${file}_SortMeRNA
#SBATCH --output=${WORK_DIR}/sort_rna/logs.o/${file}.o
#SBATCH --error=${WORK_DIR}/sort_rna/logs.e/${file}.e

# Author: Tina Yang
# Date: Feb 6, 2024
# Run SortMeRNA

# Set up envs
mamba activate rna-seq

OUT_DIR='${WORK_DIR}/sort_rna'
FASTQ_DIR='${WORK_DIR}/trim_galore'
rRNA_DATABASE='${RRNA_DATABASE}'

## Create output folder
mkdir \"\$OUT_DIR/${file}\"
cd \"\$OUT_DIR/${file}\"

sortmerna --ref \"\$rRNA_DATABASE/smr_v4.3_default_db.fasta\" --reads \"\$FASTQ_DIR/${file}_1_val_1.fq.gz\" --reads \"\$FASTQ_DIR/${file}_2_val_2.fq.gz\" -workdir \"\$OUT_DIR/${file}\" -fastx --other --paired_in --out2
" > "${file}_sortRNA.sh"; done

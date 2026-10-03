#!/bin/bash
# Job generator: writes one SLURM script per sample (<sample>_featureCounts.sh), then submit with sbatch.

# ── Paths: edit before running ───────────────────────────────────────────────
WORK_DIR=/path/to/rnaseq                             # contains list.txt (one sample ID per line)
ANNOTATION_PATH=/path/to/references/STAR/GRCh38             # contains gencode.v44.basic.annotation.gtf
# ─────────────────────────────────────────────────────────────────────────────

cd ${WORK_DIR}
mkdir featureCounts
cd featureCounts
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
#SBATCH --job-name=${file}_featureCounts
#SBATCH --output=${WORK_DIR}/featureCounts/logs.o/${file}.o
#SBATCH --error=${WORK_DIR}/featureCounts/logs.e/${file}.e

# Author: Tina Yang
# Date: Feb 7, 2024
# FeatureCount

# Set up envs
mamba activate rna-seq

# Dataset
OUT_DIR='${WORK_DIR}/featureCounts'
STAR_PATH='${WORK_DIR}/STAR'
ANNOTATION_PATH='${ANNOTATION_PATH}'

cd \"\$OUT_DIR\"

featureCounts \\
\"\$STAR_PATH/${file}/${file}_PASS2_resultAligned.sortedByCoord.out.bam\" \\
-a \"\$ANNOTATION_PATH/gencode.v44.basic.annotation.gtf\" \\
-p \\
-o ${file}.featureCounts.txt \\
-T 64
" > "${file}_featureCounts.sh"; done

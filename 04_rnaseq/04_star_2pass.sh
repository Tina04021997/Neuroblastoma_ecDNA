#!/bin/bash
# Writes one SLURM script per sample (<sample>_STAR.sh) into $STAR_DIR, then submit with sbatch.
# Run the generated jobs ONE AT A TIME: step 3 rebuilds the shared $STAR_REF/SJ_index directory.
# STAR 2.7.10b, 2-pass mapping of the non-rRNA reads from 03_sortmerna.sh against the GRCh38 primary assembly
# with the GENCODE v44 basic GTF ($STAR_REF, built once with the commented step 1). Original run Feb 7, 2024.
# Usage: bash 04_star_2pass.sh   (samples from $RNA_LIST of the current batch, config.sh)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

require_file "${RNA_GENOME_FA}" "${GTF}"
mkdir -p "${STAR_DIR}/logs.o" "${STAR_DIR}/logs.e"
cd "${STAR_DIR}"

for file in $(sample_ids); do
    require_file "${SORTRNA_DIR}/${file}/out/other_fwd.fq" "${SORTRNA_DIR}/${file}/out/other_rev.fq"
    cat > "${file}_STAR.sh" <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -t 100:00:00
#SBATCH --mem=100G
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -A ${SLURM_ACCOUNT}
#SBATCH --job-name=${file}_STAR
#SBATCH --output=${STAR_DIR}/logs.o/${file}.o
#SBATCH --error=${STAR_DIR}/logs.e/${file}.e

# Author: Tina Yang — STAR 2-pass mapping (PE) for ${file}; original run Feb 7, 2024
# Steps: 1. index genome with annotation (once) 2. 1-pass mapping 3. re-index with SJ.out.tab 4. 2-pass mapping
set -euo pipefail
source "${REPO_DIR}/config.sh"
activate_env rna-seq

#  # Step 1 — indexing (run once, not per sample)
#  STAR --runMode genomeGenerate --genomeDir "${STAR_REF}" \\
#       --genomeFastaFiles "${RNA_GENOME_FA}" --sjdbGTFfile "${GTF}" --runThreadN 32

mkdir -p "${STAR_DIR}/${file}"
cd "${STAR_DIR}/${file}"
log_version "STAR --version" ${file}.star.versions.txt

# Step 2 — 1-pass mapping
STAR \\
--genomeDir "${STAR_REF}" \\
--readFilesIn "${SORTRNA_DIR}/${file}/out/other_fwd.fq" "${SORTRNA_DIR}/${file}/out/other_rev.fq" \\
--outFileNamePrefix ${file}_PASS1_result \\
--runThreadN 64

# Step 3 — re-index with this sample's splice junctions (shared SJ_index: run samples one at a time)
STAR \\
--runMode genomeGenerate \\
--genomeDir "${STAR_REF}/SJ_index" \\
--genomeFastaFiles "${RNA_GENOME_FA}" \\
--sjdbGTFfile "${GTF}" \\
--runThreadN 64 \\
--sjdbFileChrStartEnd "${STAR_DIR}/${file}/${file}_PASS1_resultSJ.out.tab"

# Step 4 — 2-pass mapping, coordinate-sorted BAM
STAR \\
--genomeDir "${STAR_REF}/SJ_index" \\
--readFilesIn "${SORTRNA_DIR}/${file}/out/other_fwd.fq" "${SORTRNA_DIR}/${file}/out/other_rev.fq" \\
--outFileNamePrefix ${file}_PASS2_result \\
--runThreadN 64 \\
--outSAMtype BAM SortedByCoordinate \\
--outSAMstrandField intronMotif
EOF
    echo "wrote ${file}_STAR.sh"
done

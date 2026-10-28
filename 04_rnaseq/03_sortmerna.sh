#!/bin/bash
# Writes one SLURM script per sample (<sample>_sortRNA.sh) into $SORTRNA_DIR, then submit with sbatch.
# SortMeRNA 4.3.6 (smr_v4.3_default_db) removes rRNA reads from the trimmed FASTQs; the non-rRNA reads
# ($SORTRNA_DIR/<sample>/out/other_fwd.fq, other_rev.fq) go to STAR (04_star_2pass.sh). Original run Feb 6, 2024.
# Usage: bash 03_sortmerna.sh   (samples from $RNA_LIST of the current batch, config.sh)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

require_file "${RRNA_DB}/smr_v4.3_default_db.fasta"
mkdir -p "${SORTRNA_DIR}/logs.o" "${SORTRNA_DIR}/logs.e"
cd "${SORTRNA_DIR}"

ids=$(sample_ids) || exit 1          # a failed substitution inside "for" would not stop the script

for file in $ids; do
    require_file "${TRIM_DIR}/${file}_1_val_1.fq.gz" "${TRIM_DIR}/${file}_2_val_2.fq.gz"
    cat > "${file}_sortRNA.sh" <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -t 100:00:00
#SBATCH --mem=100G
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -A ${SLURM_ACCOUNT}
#SBATCH --job-name=${file}_SortMeRNA
#SBATCH --output=${SORTRNA_DIR}/logs.o/${file}.o
#SBATCH --error=${SORTRNA_DIR}/logs.e/${file}.e

# Author: Tina Yang — SortMeRNA for ${file}; original run Feb 6, 2024
set -euo pipefail
source "${REPO_DIR}/config.sh"
activate_env rna-seq

mkdir -p "${SORTRNA_DIR}/${file}"
cd "${SORTRNA_DIR}/${file}"
log_version "sortmerna --version" ${file}.sortmerna.versions.txt

sortmerna --ref "${RRNA_DB}/smr_v4.3_default_db.fasta" \\
    --reads "${TRIM_DIR}/${file}_1_val_1.fq.gz" --reads "${TRIM_DIR}/${file}_2_val_2.fq.gz" \\
    -workdir "${SORTRNA_DIR}/${file}" -fastx --other --paired_in --out2
EOF
    echo "wrote ${file}_sortRNA.sh"
done

#!/bin/bash
# Writes one SLURM script per sample (<sample>_featureCounts.sh) into $FC_DIR, then submit with sbatch.
# featureCounts (Subread 2.0.6) on the unflagged STAR 2-pass BAM, GENCODE v44 basic GTF, gene level.
#   -p --countReadPairs : count FRAGMENTS (read pairs), not reads — required with Subread >= 2.0.2 (audit point 4;
#                         the original Feb 2024 run used -p alone, which counts reads in that version; re-counted Oct 2026)
#   -B                  : both ends of a pair must align      -C : discard chimeric pairs
#   no duplicate filtering (duplicates retained; see REANALYSIS.md)
# Output: $FC_DIR/<sample>.featureCounts.txt (column 7 = counts) -> 08b_featurecounts_sanity_check.sh, 08_merge_counts.sh
# Usage : bash 07_featurecounts.sh   (samples from $RNA_LIST of the current batch, config.sh)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

require_file "${GTF}"
mkdir -p "${FC_DIR}/logs.o" "${FC_DIR}/logs.e"
cd "${FC_DIR}"

for file in $(sample_ids); do
    bam="${STAR_DIR}/${file}/${file}_PASS2_resultAligned.sortedByCoord.out.bam"
    require_file "$bam"
    cat > "${file}_featureCounts.sh" <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH --mem=30G
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -A ${SLURM_ACCOUNT}
#SBATCH --job-name=${file}_featureCounts
#SBATCH --output=${FC_DIR}/logs.o/${file}.o
#SBATCH --error=${FC_DIR}/logs.e/${file}.e

# Author: Tina Yang — featureCounts for ${file}; original run Feb 7, 2024, fragment re-count Oct 2026
set -euo pipefail
source "${REPO_DIR}/config.sh"
activate_env rna-seq

cd "${FC_DIR}"
log_version "featureCounts -v" ${file}.featureCounts.versions.txt

featureCounts \\
    "${bam}" \\
    -a "${GTF}" \\
    -p --countReadPairs -B -C \\
    -o ${file}.featureCounts.txt \\
    -T 2
EOF
    echo "wrote ${file}_featureCounts.sh"
done

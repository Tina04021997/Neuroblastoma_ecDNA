#!/bin/bash
# Writes one SLURM script per sample (<sample>_Delly.sh) into $DELLY_DIR, then submit with sbatch.
# DELLY v1.7.2 static binary ($DELLY_BIN), run on the duplicate-marked WGS BAM (01_wgs_alignment), -q 20 -s 15,
# centromere exclusion list $DELLY_EXCLUDE. Original run: Feb 2026 (VCF header fileDate=20260213).
# Usage: bash 01_delly_call.sh [sample_list.txt]   (default: all 21 lines, dna_id column of samples.tsv = WGS file prefix)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

LIST="${1:-}"
if [ -z "$LIST" ]; then LIST=$(mktemp); wgs_ids > "$LIST" || exit 1; fi     # default: dna_id column (WGS file prefix)
require_file "${DELLY_BIN}" "${DELLY_EXCLUDE}" "${WGS_REF_FASTA}"
mkdir -p "${DELLY_DIR}/logs.o" "${DELLY_DIR}/logs.e"
cd "${DELLY_DIR}"

ids=$(sample_ids "$LIST") || exit 1          # a failed substitution inside "for" would not stop the script

for sample in $ids; do
    require_file "${WGS_ALIGN_DIR}/${sample}/${sample}_tumor_mkdp.bam"
    cat > "${sample}_Delly.sh" <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 20:00:00
#SBATCH --mem=60G
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -A ${SLURM_ACCOUNT}
#SBATCH --job-name=${sample}_Delly
#SBATCH --output=${DELLY_DIR}/logs.o/${sample}.o
#SBATCH --error=${DELLY_DIR}/logs.e/${sample}.e

set -euo pipefail
source "${REPO_DIR}/config.sh"

mkdir -p "${DELLY_DIR}/${sample}"
cd "${DELLY_DIR}/${sample}"
log_version "${DELLY_BIN} --version" "${sample}.delly.versions.txt"

"${DELLY_BIN}" call \\
    -x "${DELLY_EXCLUDE}" \\
    -g "${WGS_REF_FASTA}" \\
    -q 20 \\
    -s 15 \\
    "${WGS_ALIGN_DIR}/${sample}/${sample}_tumor_mkdp.bam" > "${sample}.vcf"
EOF
    echo "wrote ${sample}_Delly.sh"
done

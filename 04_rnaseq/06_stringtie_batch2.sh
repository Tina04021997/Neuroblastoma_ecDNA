#!/bin/bash
# StringTie 2.2.1 — sequencing batch 2 (11 lines: KAN, KANR, KCNR, KPNYN, LAN1, LANS, NB10, NB1643, NB69, SKNBE2, SKNrI).
# Writes one SLURM script per sample (<sample>_stringtie.sh) into $STRINGTIE_DIR, then submit with sbatch.
# Input : $STAR_DIR/<sample>/<sample>_PASS2_resultAligned.sortedByCoord.out.bam (unflagged STAR 2-pass BAM; batch 2 was
#         not run through Picard MarkDuplicates).
# Mode  : reference-ONLY quantification (-e): expression is estimated for the GTF transcripts only, no novel assembly.
#         Batch 1 was run in assembly mode on the Picard-flagged BAM — see 06_stringtie_batch1.sh. Gene-level FPKM for
#         annotated genes is unaffected by the mode (REANALYSIS.md).
# Output: $STRINGTIE_DIR/<sample>/<sample>.gene.abundance.txt (-A), .transcripts.gtf, .coverage.gtf, .ballgown/.
#         Original run Feb 2024 (SLURM, TSCC).
# Usage : bash 06_stringtie_batch2.sh   (samples from $RNA_LIST of batch 2, config.sh)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

require_file "${GTF}"
mkdir -p "${STRINGTIE_DIR}/logs.o" "${STRINGTIE_DIR}/logs.e"
cd "${STRINGTIE_DIR}"

for file in $(sample_ids); do
    bam="${STAR_DIR}/${file}/${file}_PASS2_resultAligned.sortedByCoord.out.bam"
    require_file "$bam"
    cat > "${file}_stringtie.sh" <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -t 30:00:00
#SBATCH --mem=200G
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -A ${SLURM_ACCOUNT}
#SBATCH --job-name=${file}_stringtie
#SBATCH --output=${STRINGTIE_DIR}/logs.o/${file}.o
#SBATCH --error=${STRINGTIE_DIR}/logs.e/${file}.e

# Author: Tina Yang — StringTie (reference-only mode, -e) for ${file}, batch 2
set -euo pipefail
source "${REPO_DIR}/config.sh"
activate_env rna-seq

mkdir -p "${STRINGTIE_DIR}/${file}"
cd "${STRINGTIE_DIR}/${file}"
log_version "stringtie --version" ${file}.stringtie.versions.txt

stringtie \\
    "${bam}" \\
    -G "${GTF}" \\
    -o ${file}.transcripts.gtf \\
    -A ${file}.gene.abundance.txt \\
    -C ${file}.coverage.gtf \\
    -b ${file}.ballgown \\
    -e
EOF
    echo "wrote ${file}_stringtie.sh"
done

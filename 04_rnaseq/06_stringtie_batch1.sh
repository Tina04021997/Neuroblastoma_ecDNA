#!/bin/bash
# StringTie 2.2.1 — sequencing batch 1 (10 lines: CHP134, CHP212, IMR32, Kelly, LAN6, NBLS, NGP, SKNAS, SKNSH, TR14).
# Writes one PBS script per sample (<sample>_stringtie.sh) into $STRINGTIE_DIR, then submit with qsub.
# Input : $MARKDP_DIR/<sample>.bam (Picard-flagged STAR BAM, 05_picard_markduplicates.sh); StringTie ignores the
#         duplicate flag (verified for CHP212: identical gene FPKM with and without flagging, REANALYSIS.md).
# Mode  : reference-GUIDED ASSEMBLY (no -e): novel transcripts (STRG.*) are assembled in addition to the GTF genes.
#         Batch 2 was run in reference-only mode (-e) on the unflagged STAR BAM — see 06_stringtie_batch2.sh.
# Output: <sample>.gene.abundance.txt (-A; FPKM per gene, read by 05_differential_expression/03_export_FPKM_tables.R),
#         <sample>.transcripts.gtf, <sample>.coverage.gtf, <sample>.ballgown/. Original run Feb 2024 (PBS, old TSCC).
# Usage : bash 06_stringtie_batch1.sh   (samples from $RNA_LIST of batch 1, config.sh)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

require_file "${GTF}"
mkdir -p "${STRINGTIE_DIR}/logs_o" "${STRINGTIE_DIR}/logs_e"
cd "${STRINGTIE_DIR}"

for file in $(sample_ids); do
    require_file "${MARKDP_DIR}/${file}.bam"
    cat > "${file}_stringtie.sh" <<EOF
#!/bin/bash
#PBS -q ${PBS_QUEUE}
#PBS -A ${PBS_ACCOUNT}
#PBS -l nodes=1:ppn=32
#PBS -l walltime=300:00:00
#PBS -o ${STRINGTIE_DIR}/logs_o
#PBS -e ${STRINGTIE_DIR}/logs_e
#PBS -N ${file}_stringtie
#PBS -V

# Author: Tina Yang — StringTie (assembly mode) for ${file}, batch 1
set -euo pipefail
source "${REPO_DIR}/config.sh"
activate_env rna-seq

cd "${STRINGTIE_DIR}"
log_version "stringtie --version" ${file}.stringtie.versions.txt

stringtie \\
    "${MARKDP_DIR}/${file}.bam" \\
    -G "${GTF}" \\
    -o ${file}.transcripts.gtf \\
    -A ${file}.gene.abundance.txt \\
    -C ${file}.coverage.gtf \\
    -b ${file}.ballgown
EOF
    echo "wrote ${file}_stringtie.sh"
done

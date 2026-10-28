#!/bin/bash
# Writes one SLURM script per sample (<sample>_AA.sh) into $ECDNA_DIR, then submit with sbatch.
# PrepareAA.py from a git clone of AmpliconSuite-pipeline ($AMPSUITE_DIR; v1.2.2 at run time, AmpliconClassifier 1.1.1),
# with python / samtools / CNVkit (0.9.10 at run time) / AmpliconArchitect (ampliconarchitectlib) taken from the conda env `ampsuite`
# (envs/ampsuite.yml is a later export of that env; see envs/README.md).
# GRCh38, CNVkit seeds, --cngain 4.5 --cnsize_min 50000, --downsample -1 (no downsampling), 32 threads.
# Per-run versions are written by the pipeline to <sample>/<sample>_run_metadata.json.
# Usage: bash 01_run_ampliconsuite.sh [sample_list.txt]   (default: all 21 lines, dna_id column of samples.tsv = WGS file prefix)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

LIST="${1:-}"
if [ -z "$LIST" ]; then LIST=$(mktemp); wgs_ids > "$LIST" || exit 1; fi     # default: dna_id column (WGS file prefix)
mkdir -p "${ECDNA_DIR}/logs.o" "${ECDNA_DIR}/logs.e"
cd "${ECDNA_DIR}"

ids=$(sample_ids "$LIST") || exit 1          # a failed substitution inside "for" would not stop the script

for file in $ids; do
    require_file "${WGS_ALIGN_DIR}/${file}/${file}_tumor_mkdp.bam"
    cat > "${file}_AA.sh" <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -t 200:00:00
#SBATCH --mem=1000G
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -A ${SLURM_ACCOUNT}
#SBATCH --job-name=${file}_ecDNA
#SBATCH --output=${ECDNA_DIR}/logs.o/${file}_AA.o
#SBATCH --error=${ECDNA_DIR}/logs.e/${file}_AA.e

# Author: Tina Yang — AmpliconSuite (PrepareAA -> AA -> AC) for ${file}; original run Feb 2024
set -euo pipefail
source "${REPO_DIR}/config.sh"
activate_env ampsuite

mkdir -p "${ECDNA_DIR}/${file}"
cd "${ECDNA_DIR}/${file}"
log_version "python ${AMPSUITE_DIR}/PrepareAA.py --version; ${AMPSUITE_ENV}/bin/cnvkit.py version" ${file}.versions.txt

${AMPSUITE_DIR}/PrepareAA.py -s ${file} -t 32 --ref GRCh38 -o "${ECDNA_DIR}/${file}" \\
    --cnvkit_dir ${AMPSUITE_ENV}/bin/cnvkit.py \\
    --bam "${WGS_ALIGN_DIR}/${file}/${file}_tumor_mkdp.bam" \\
    --aa_python_interpreter ${AMPSUITE_ENV}/bin/python3 \\
    --samtools_path ${AMPSUITE_ENV}/bin/samtools \\
    --cngain 4.5 --cnsize_min 50000 --downsample -1 --run_AA --run_AC
EOF
    echo "wrote ${file}_AA.sh"
done

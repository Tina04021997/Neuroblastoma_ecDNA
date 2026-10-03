#!/bin/bash
# Job generator: writes one SLURM script per sample (<sample>_AA.sh), then submit with sbatch.

# ── Paths: edit before running ───────────────────────────────────────────────
WORK_PATH=/path/to/ecDNA_analysis                           # AmpliconSuite output; contains list.txt (one sample ID per line)
BAM_DIR=/path/to/wgs/alignment                              # ${BAM_DIR}/<sample>/<sample>_tumor_mkdp.bam (01_wgs_alignment)
AMPSUITE_DIR=/path/to/AmpliconSuite-pipeline                # contains PrepareAA.py
AMPSUITE_ENV=/path/to/conda/envs/ampsuite                   # conda env providing cnvkit.py, python3, samtools
# ─────────────────────────────────────────────────────────────────────────────

cd $WORK_PATH

for file in $(cat $WORK_PATH/list.txt); do echo -e "#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -t 200:00:00
#SBATCH --mem=1000G
#SBATCH -p <partition>
#SBATCH -A <account>
#SBATCH --job-name=${file}_ecDNA
#SBATCH --output=${WORK_PATH}/logs.o/${file}_AA.o
#SBATCH --error=${WORK_PATH}/logs.e/${file}_AA.e

# Author: Tina Yang
# Date: Feb 9, 2024
# Run AA + AC


OUT_DIR='${WORK_PATH}'
BAM_DIR='${BAM_DIR}'


mkdir \"\$OUT_DIR/${file}\"
cd \"\$OUT_DIR/${file}\"

mamba activate ampsuite

${AMPSUITE_DIR}/PrepareAA.py -s ${file} -t 64 --ref GRCh38 -o \"\$OUT_DIR/${file}\" --cnvkit_dir ${AMPSUITE_ENV}/bin/cnvkit.py --bam \"\$BAM_DIR/${file}/${file}_tumor_mkdp.bam\" --aa_python_interpreter ${AMPSUITE_ENV}/bin/python3 --samtools_path ${AMPSUITE_ENV}/bin/samtools --cngain 4.5 --cnsize_min 50000 --downsample -1 --run_AA --run_AC
" > "${file}_AA.sh"; done

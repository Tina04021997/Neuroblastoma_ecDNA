#!/bin/bash
# Job generator: writes one SLURM script per sample (<sample>_Delly.sh), then submit with sbatch.

# ── Paths: edit before running ───────────────────────────────────────────────
sample_file=/path/to/sample_list.txt                        # one sample ID per line
OUT_DIR=/path/to/delly                                      # output: ${OUT_DIR}/<sample>/<sample>.vcf
BAM_DIR=/path/to/wgs/alignment                              # ${BAM_DIR}/<sample>/<sample>_tumor_mkdp.bam
DELLY=/path/to/delly_v1.7.2_linux_x86_64bit                 # DELLY v1.7.2 static binary
EXCLUDE=/path/to/references/Delly/hg38.cen                  # centromere exclusion list
REF_FASTA=/path/to/references/GRCh38.d1.vd1/GRCh38.d1.vd1.fa
CONDA_BASE=/path/to/conda                                   # output of `conda info --base`
# ─────────────────────────────────────────────────────────────────────────────

# Loop through each sample in the sample file
while IFS= read -r sample; do
    echo "#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 20:00:00
#SBATCH --mem=60G
#SBATCH -p <partition>
#SBATCH -A <account>
#SBATCH --job-name=${sample}_Delly
#SBATCH --output=./logs.o/${sample}.o
#SBATCH --error=./logs.e/${sample}.e

cd ${OUT_DIR}
mkdir -p ${sample}
cd ${sample}

source ${CONDA_BASE}/etc/profile.d/conda.sh
source ~/.bashrc
conda activate ASCAT

${DELLY} call \\
    -x ${EXCLUDE} \\
    -g ${REF_FASTA} \\
    -q 20 \\
    -s 15 \\
    ${BAM_DIR}/${sample}/${sample}_tumor_mkdp.bam > ${sample}.vcf

" > "${sample}_Delly.sh"
done < "$sample_file"

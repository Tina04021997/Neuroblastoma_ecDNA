#!/bin/bash
# Job generator: writes one PBS script per sample (<sample>_stringtie.sh), then submit with qsub.

# ── Paths: edit before running ───────────────────────────────────────────────
SAMPLE_LIST=/path/to/sample_list.txt                        # one sample ID per line
MARKDP_DIR=/path/to/rnaseq/markDP                           # output of 05_picard_markduplicates.sh
GTF=/path/to/references/STAR/GRCh38/gencode.v44.basic.annotation.gtf
OUT_DIR=/path/to/rnaseq/stringtie
CONDA_BASE=/path/to/conda                                   # output of `conda info --base`
# ─────────────────────────────────────────────────────────────────────────────

for file in $(cat ${SAMPLE_LIST}); do echo -e "#!/bin/bash
#PBS -q <queue>
#PBS -A <account>
#PBS -l nodes=1:ppn=32
#PBS -l walltime=300:00:00
#PBS -o ${OUT_DIR}/logs_o
#PBS -e ${OUT_DIR}/logs_e
#PBS -N ${file}_stringtie
#PBS -V

# Set up conda envs
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate rna-seq

cd ${OUT_DIR}/

stringtie \\
    ${MARKDP_DIR}/${file}.bam \\
    -G ${GTF} \\
    -o ${file}.transcripts.gtf \\
    -A ${file}.gene.abundance.txt \\
    -C ${file}.coverage.gtf \\
    -b ${file}.ballgown \\
" > "${file}_stringtie.sh"; done

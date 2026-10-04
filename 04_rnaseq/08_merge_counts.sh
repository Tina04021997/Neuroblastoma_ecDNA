#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -t 100:00:00
#SBATCH --mem=100G
#SBATCH -p <partition>
#SBATCH -A <account>
#SBATCH --job-name=FC_process
#SBATCH --output=logs.o/process.o
#SBATCH --error=logs.e/process.e

# Author: Tina Yang
# Date: Feb 7, 2024
# Process FeatureCounts
# Run order: sections 1-3 below -> 08b_featurecounts_sanity_check.sh (creates reference_column.txt) -> section 4.

# Set up envs
mamba activate rna-seq

# ── Paths: edit before running ───────────────────────────────────────────────
WORK_DIR=/path/to/rnaseq                      
OUT_DIR=${WORK_DIR}/featureCounts
STAR_DIR=${WORK_DIR}/STAR                   
# ─────────────────────────────────────────────────────────────────────────────

cd $OUT_DIR

# 1. Generate individual counts (geneID:counts)
for file in $(cat ../list.txt); do cut -f1,7 ${file}.featureCounts.txt > ${file}_counts.txt; done

# 2. Change the column name of the counts.txt
for file in $(cat $WORK_DIR/list.txt); do sed -i "s,$STAR_DIR/${file}/${file}_PASS2_resultAligned.sortedByCoord.out.bam,${file},g" ${file}_counts.txt; done

# 3. Remove the first row
for file in $(cat $WORK_DIR/list.txt); do sed -i '1d' ${file}_counts.txt; done


# 4. Concatenate all the individual files into one counts.txt
cp reference_column.txt merged_counts.txt

# Iterate over the samples
for file in $(cat ../list.txt); do
# Extract the second column from the current file and append it to the output file
awk '{print $2}' ${file}_counts.txt | paste merged_counts.txt - > temp.txt
mv temp.txt merged_counts.txt
done

# Remove the character after the decimal in geneid
awk -F'\t' 'BEGIN {OFS = FS} NR > 1 {gsub(/\..*/, "", $1)} 1' merged_counts.txt > final_counts.txt

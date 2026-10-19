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


# 5. Join the two sequencing batches on Geneid -> 21-sample matrix for DESeq2.
#    Both inputs are version-stripped final_counts.txt files from this script, one per batch.
BATCH1=/path/to/batch1/final_counts.txt   # CHP134 ... TR14   (10 lines)
BATCH2=/path/to/batch2/final_counts.txt   # KAN ... SKNrI     (11 lines)
OUT21=counts_21lines.txt

cmp <(cut -f1 "$BATCH1") <(cut -f1 "$BATCH2") || { echo "ERROR: gene ID columns differ between batches" >&2; exit 1; }
paste "$BATCH1" <(cut -f2- "$BATCH2") > "$OUT21"
[ "$(awk -F'\t' 'NR>1 && $1 ~ /\./' "$OUT21" | wc -l)" -eq 0 ] || { echo "ERROR: versioned IDs in $OUT21" >&2; exit 1; }
echo "$OUT21: $(($(wc -l < "$OUT21")-1)) genes x $(($(head -1 "$OUT21" | awk -F'\t' '{print NF}')-1)) samples"

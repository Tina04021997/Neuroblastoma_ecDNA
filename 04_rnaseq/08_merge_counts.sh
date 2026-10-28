#!/bin/bash
# Build the gene x sample count matrix from the per-sample featureCounts outputs of ONE batch (steps 1-4), then join
# the two batches on Ensembl gene ID into the 21-sample matrix for DESeq2 (step 5).
# Steps 1-4 run from $FC_DIR of the current batch (config.sh); step 5 needs both batches' final_counts.txt (set below).
# The paste in step 4 is guarded: 08b_featurecounts_sanity_check.sh is run automatically, writes reference_column.txt and
# compares every sample's gene-ID column with it (content AND order); the run aborts on any mismatch, so a shuffled input
# cannot change counts. Usage: bash 08_merge_counts.sh            # steps 1-4 for this batch
#        bash 08_merge_counts.sh --join     # step 5 only: BATCH1_COUNTS + BATCH2_COUNTS -> counts_21lines.txt

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"      # absolute: step 4 runs after cd "$FC_DIR"
source "${HERE}/../config.sh"

if [ "${1:-}" != "--join" ]; then
    cd "${FC_DIR}"
    # 1. gene ID + count column per sample
    ids=$(sample_ids) || exit 1          # a failed substitution inside "for" would not stop the script
    for file in $ids; do
        require_file "${file}.featureCounts.txt"
        cut -f1,7 "${file}.featureCounts.txt" > "${file}_counts.txt"
    done
    # 2. header: replace the BAM path by the sample ID;  3. drop the featureCounts comment line
    ids=$(sample_ids) || exit 1          # a failed substitution inside "for" would not stop the script
    for file in $ids; do
        sed -i "s,${STAR_DIR}/${file}/${file}_PASS2_resultAligned.sortedByCoord.out.bam,${file},g" "${file}_counts.txt"
        sed -i '1d' "${file}_counts.txt"
    done
    # 4. paste columns, guarded by the gene-ID check (08b writes reference_column.txt and exits 1 on any mismatch)
    bash "${HERE}/08b_featurecounts_sanity_check.sh" || { echo "ERROR: gene-ID columns differ between samples; nothing merged" >&2; exit 1; }
    [ -s reference_column.txt ] || { echo "ERROR: reference_column.txt missing" >&2; exit 1; }
    cp reference_column.txt merged_counts.txt
    ids=$(sample_ids) || exit 1          # a failed substitution inside "for" would not stop the script
    for file in $ids; do
        cmp -s <(cut -f1 "${file}_counts.txt") reference_column.txt \
            || { echo "ERROR: gene ID column of ${file}_counts.txt differs from reference_column.txt (order or content)" >&2; exit 1; }
        paste merged_counts.txt <(cut -f2 "${file}_counts.txt") > temp.txt
        mv temp.txt merged_counts.txt
    done
    # strip the Ensembl version suffix (ENSG00000xxxxxx.12 -> ENSG...); DESeq2 maps symbols on unversioned IDs
    awk -F'\t' 'BEGIN {OFS = FS} NR > 1 {gsub(/\..*/, "", $1)} 1' merged_counts.txt > final_counts.txt
    echo "final_counts.txt: $(($(wc -l < final_counts.txt)-1)) genes x $(($(head -1 final_counts.txt | awk -F'\t' '{print NF}')-1)) samples"
    exit 0
fi

# 5. Join the two sequencing batches on Geneid -> 21-sample matrix for 05_differential_expression/01_deseq2_ecDNA.R.
#    Both inputs are version-stripped final_counts.txt files produced above, one per batch (paths in config.sh).
require_file "${BATCH1_COUNTS}" "${BATCH2_COUNTS}"
OUT21="${COUNTS_21}"
cmp <(cut -f1 "${BATCH1_COUNTS}") <(cut -f1 "${BATCH2_COUNTS}") || { echo "ERROR: gene ID columns differ between batches" >&2; exit 1; }
paste "${BATCH1_COUNTS}" <(cut -f2- "${BATCH2_COUNTS}") > "$OUT21"
[ "$(awk -F'\t' 'NR>1 && $1 ~ /\./' "$OUT21" | wc -l)" -eq 0 ] || { echo "ERROR: versioned IDs in $OUT21" >&2; exit 1; }
echo "$OUT21: $(($(wc -l < "$OUT21")-1)) genes x $(($(head -1 "$OUT21" | awk -F'\t' '{print NF}')-1)) samples"

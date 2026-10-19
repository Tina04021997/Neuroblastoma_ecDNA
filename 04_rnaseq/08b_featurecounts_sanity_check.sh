#!/bin/bash

# Make sure the first column of featureCounts output is all the same,
# and write reference_column.txt (gene IDs) used by 08_merge_counts.sh.

# ── Paths: edit before running ───────────────────────────────────────────────
SAMPLE="KAN"                                              
WORK_DIR=/path/to/rnaseq/featureCounts
# ─────────────────────────────────────────────────────────────────────────────

cd $WORK_DIR

cut -d $'\t' --complement -f 2 ${SAMPLE}_counts.txt > reference_column.txt

status=0
for file in *_counts.txt; do
    if ! cmp -s <(cut -f1 "$file") reference_column.txt; then
        echo "MISMATCH: $file gene IDs differ from reference_column.txt" >&2
        status=1
    fi
done
[ "$status" -eq 0 ] && echo "OK: all $(ls *_counts.txt | wc -l) samples share the reference gene-ID column"
exit "$status"

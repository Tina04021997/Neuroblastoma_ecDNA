#!/bin/bash

# Make sure the first column of featureCounts output is all the same,
# and write reference_column.txt (gene IDs) used by 08_merge_counts.sh.

# ── Paths: edit before running ───────────────────────────────────────────────
SAMPLE="KAN"                                              
WORK_DIR=/path/to/rnaseq/featureCounts
# ─────────────────────────────────────────────────────────────────────────────

cd $WORK_DIR

cut -d $'\t' --complement -f 2 ${SAMPLE}_counts.txt > reference_column.txt

for file in *_counts.txt; do cut -d $'\t' --complement -f 2 $file > ${file%_counts.txt}_column.txt; done

for file in *_column.txt; do diff -q "$file" reference_column.txt; done

rm *_column.txt

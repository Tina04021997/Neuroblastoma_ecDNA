#!/bin/bash
# Check that every <sample>_counts.txt in $FC_DIR has the same gene-ID column (content AND order) and write
# reference_column.txt, which 08_merge_counts.sh compares against before pasting columns together.
# Run after step 1-3 of 08_merge_counts.sh have produced the <sample>_counts.txt files (or re-run 08 after this).
# Usage: bash 08b_featurecounts_sanity_check.sh [reference_sample]   (default: first sample in $RNA_LIST)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

cd "${FC_DIR}"
REF_SAMPLE="${1:-$(sample_ids | head -1)}"
require_file "${REF_SAMPLE}_counts.txt"
cut -f1 "${REF_SAMPLE}_counts.txt" > reference_column.txt

status=0
for file in *_counts.txt; do
    if ! cmp -s <(cut -f1 "$file") reference_column.txt; then
        echo "MISMATCH: $file gene IDs differ from reference_column.txt" >&2
        status=1
    fi
done
[ "$status" -eq 0 ] && echo "OK: all $(ls *_counts.txt | wc -l) samples share the reference gene-ID column ($(wc -l < reference_column.txt) rows)"
exit "$status"

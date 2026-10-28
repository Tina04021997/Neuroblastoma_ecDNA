#!/bin/bash
# Keep PASS calls on chr1-22 and chrX (header removed) as input for 03_format_sv_for_shatterseek.py.
# Column-based filter: $1 = CHROM, $7 = FILTER (tab-delimited VCF). Requires GNU or BSD awk.
# Reads $DELLY_DIR/<sample>/<sample>.vcf (01_delly_call.sh), writes $DELLY_DIR/<sample>_sv.txt.
# Usage: bash 02_filter_delly_vcf.sh [sample_list.txt]   (default: all 21 lines, dna_id column of samples.tsv = WGS file prefix)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

LIST="${1:-}"
if [ -z "$LIST" ]; then LIST=$(mktemp); wgs_ids > "$LIST" || exit 1; fi     # default: dna_id column (WGS file prefix)
cd "${DELLY_DIR}"

ids=$(sample_ids "$LIST") || exit 1          # a failed substitution inside "for" would not stop the script

for sample in $ids; do
    in="${sample}/${sample}.vcf"
    out="${sample}_sv.txt"
    require_file "$in"
    awk -F'\t' '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X)$/ && $7=="PASS"' "$in" > "$out"
    [ -s "$out" ] || { echo "ERROR: no PASS SVs kept for $sample" >&2; exit 1; }
    echo "$sample: $(wc -l < "$out") PASS SVs"
done

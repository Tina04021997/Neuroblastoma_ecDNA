#!/bin/bash
set -euo pipefail
# Keep PASS calls on chr1-22 and chrX (header removed) as input for 03_format_sv_for_shatterseek.py.
# Column-based filter: $1 = CHROM, $7 = FILTER (tab-delimited VCF). Requires GNU or BSD awk.
# Usage: bash 02_filter_delly_vcf.sh sample_list.txt

while read -r sample; do
    in="${sample}/${sample}.vcf"
    out="${sample}_sv.txt"
    [ -s "$in" ] || { echo "ERROR: missing or empty $in" >&2; exit 1; }
    awk -F'\t' '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X)$/ && $7=="PASS"' "$in" > "$out"
    [ -s "$out" ] || { echo "ERROR: no PASS SVs kept for $sample" >&2; exit 1; }
    echo "$sample: $(wc -l < "$out") PASS SVs"
done < "$1"

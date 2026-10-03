#!/bin/bash
# Keep PASS calls on chr1-22 and chrX, header removed, as input for 03_format_sv_for_shatterseek.py.
# Taken from the original pipeline notes.
#
# NOTE: 01_delly_call.sh writes ${sample}.vcf directly (stdout), so the old
# "bcftools view sv.bcf" conversion step is not needed. Confirm which form was used.
#
# Usage: bash 02_filter_delly_vcf.sh sample_list.txt

while read -r sample; do
    awk '/^chr([1-9]|1[0-9]|2[0-2]|X)\s.*\sPASS\s/' ${sample}/${sample}.vcf > ${sample}_sv.txt
done < "$1"

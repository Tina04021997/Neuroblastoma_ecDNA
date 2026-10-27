#!/bin/bash
# Picard MarkDuplicates 3.0.0 on the STAR 2-pass BAMs — sequencing batch 1 (10 lines) only.
# Duplicates are FLAGGED, not removed (no --REMOVE_DUPLICATES); the metrics files give the per-sample duplication
# rate (QC, see REANALYSIS.md). StringTie (06_stringtie_batch1.sh) reads these BAMs but ignores the flag;
# featureCounts (07_featurecounts.sh) reads the unflagged STAR BAMs. Batch 2 skipped this step.
# Usage: bash 05_picard_markduplicates.sh   (samples from $RNA_LIST of batch 1, config.sh)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
activate_env rna-seq

mkdir -p "${MARKDP_DIR}"
log_version "picard MarkDuplicates --version" "${MARKDP_DIR}/picard.versions.txt"

for file in $(sample_ids); do
    bam="${STAR_DIR}/${file}/${file}_PASS2_resultAligned.sortedByCoord.out.bam"
    require_file "$bam"
    echo "start working on ${file} mark duplicate ..."
    picard MarkDuplicates -Xmx500g \
        --INPUT "$bam" \
        --OUTPUT "${MARKDP_DIR}/${file}.bam" \
        --METRICS_FILE "${MARKDP_DIR}/${file}.MarkDuplicates.metrics.txt"
    echo "done ${file} mark duplicate ..."
done

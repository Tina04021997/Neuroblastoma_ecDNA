#!/bin/bash

# ── Paths: edit before running ───────────────────────────────────────────────
SAMPLE_LIST=/path/to/sample_list.txt                    
STAR_DIR=/path/to/rnaseq/STAR                            
OUT_DIR=/path/to/rnaseq/markDP
# ─────────────────────────────────────────────────────────────────────────────

for file in $(cat ${SAMPLE_LIST}); do 
    (
        echo "start working on ${file} mark duplicate ..."
        picard MarkDuplicates -Xmx500g --INPUT ${STAR_DIR}/${file}/${file}_PASS2_resultAligned.sortedByCoord.out.bam --OUTPUT ${OUT_DIR}/${file}.bam --METRICS_FILE ${OUT_DIR}/${file}.MarkDuplicates.metrics.txt 
        echo "done ${file} mark duplicate ..."
    ) 
done

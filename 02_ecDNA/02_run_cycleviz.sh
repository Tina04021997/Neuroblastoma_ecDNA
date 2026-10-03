#!/bin/bash
# Fig 1B: circular visualization of the NGP MYCN ecDNA (AmpliconArchitect amplicon 1, cycle 1).
# Run inside the AmpliconArchitect output directory for NGP.

CYCLEVIZ=/path/to/CycleViz

python ${CYCLEVIZ}/CycleViz.py \
  --ref GRCh38 \
  --cycles_file NGP_amplicon1_cycles.txt \
  --cycle 1 \
  -g NGP_amplicon1_graph.txt \
  --rotate_to_min \
  --figure_size_style small \
  --gene_subset_list MYCN \
  --gene_fontsize 5 \
  --tick_fontsize 5 \
  --gene_spacing 3.2

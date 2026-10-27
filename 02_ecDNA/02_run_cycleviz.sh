#!/bin/bash
# Fig 1B: circular visualization of the NGP MYCN ecDNA (AmpliconArchitect amplicon 1, cycle 1).
# Reads the AA cycles/graph files of NGP from the final AmpliconSuite run ($ECDNA_DIR, 01_run_ampliconsuite.sh)
# and writes the figure into $ECDNA_DIR/NGP/cycleviz/. CycleViz 0.2.2 (git clone at $CYCLEVIZ_DIR).
# Usage: bash 02_run_cycleviz.sh [sample] [amplicon] [cycle]   (defaults: NGP 1 1)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"
activate_env ampsuite

SAMPLE="${1:-NGP}"; AMP="${2:-1}"; CYCLE="${3:-1}"
AA_DIR="${ECDNA_DIR}/${SAMPLE}/${SAMPLE}_AA_results"
CYCLES="${AA_DIR}/${SAMPLE}_amplicon${AMP}_cycles.txt"
GRAPH="${AA_DIR}/${SAMPLE}_amplicon${AMP}_graph.txt"
require_file "${CYCLES}" "${GRAPH}" "${CYCLEVIZ_DIR}/CycleViz.py"

OUT="${ECDNA_DIR}/${SAMPLE}/cycleviz"
mkdir -p "${OUT}"; cd "${OUT}"
log_version "python ${CYCLEVIZ_DIR}/CycleViz.py --version" "${SAMPLE}_cycleviz.versions.txt"

python "${CYCLEVIZ_DIR}/CycleViz.py" \
  --ref GRCh38 \
  --cycles_file "${CYCLES}" \
  --cycle "${CYCLE}" \
  -g "${GRAPH}" \
  --rotate_to_min \
  --figure_size_style small \
  --gene_subset_list MYCN \
  --gene_fontsize 5 \
  --tick_fontsize 5 \
  --gene_spacing 3.2

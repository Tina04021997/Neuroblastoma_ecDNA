#!/bin/bash
# Regenerate every panel of Figures 1 and 2 from the frozen inputs in data/ (audit point 11).
#
#   bash 06_figures/make_figures.sh [out_dir]          (default: figures/ in the repository root, git-ignored)
#
# Inputs (all under data/, written by the pipeline scripts named in README.md):
#   master.csv + ecDNA_summary/           -> Fig 1A  (fig1A_oncoprint.ipynb)
#   (00_build_figure_tables.py / 01_check_figure_tables.py: how the tables are built and cross-checked)
#   chr2_cn/NGP_cn.txt, SKNYrI_cn.txt     -> Fig 1B  chromosome-2 copy-number panels (fig1B_chr2_copy_number.py);
#                                            the NGP cycle diagram is drawn by 02_ecDNA/02_run_cycleviz.sh on the cluster
#   FPKM_CN_MYCN.txt                      -> Fig 2A  (fig2A_MYCN_CN_vs_FPKM.py)
#   MYCN_FPKM.txt                         -> Fig 2B  (fig2B_MYCN_FPKM_boxplot.py)
#   deseq2_gsea/deseq2_results_ecDNA.csv  -> Fig 2C  (fig2C_volcano.py) and Fig 2D/2E (02_gsea_and_RS_genes.R)
#
# Interpreters: PYTHON, RSCRIPT, JUPYTER (defaults python3, Rscript, jupyter) and JUPYTER_KERNEL (default python3; the
# kernel, not $PYTHON, runs the notebook). Python needs pandas, numpy, matplotlib, scipy, adjustText; R needs tidyverse,
# fgsea, msigdbr >= 26.1.0, ggplot2. This script does not source config.sh. Every step stops the script on failure and
# its full log is kept under <out_dir>/logs/.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(dirname "$HERE")"
OUT="$(mkdir -p "${1:-$REPO/figures}" && cd "${1:-$REPO/figures}" && pwd)"
LOG="$OUT/logs"; mkdir -p "$LOG"
PYTHON="${PYTHON:-python3}"; RSCRIPT="${RSCRIPT:-Rscript}"; JUPYTER="${JUPYTER:-jupyter}"; JUPYTER_KERNEL="${JUPYTER_KERNEL:-python3}"
DATA="$REPO/data"
for f in master.csv FPKM_CN_MYCN.txt MYCN_FPKM.txt deseq2_gsea/deseq2_results_ecDNA.csv chr2_cn/NGP_cn.txt chr2_cn/SKNYrI_cn.txt; do
  [ -s "$DATA/$f" ] || { echo "missing input: $DATA/$f" >&2; exit 1; }
done
[ -s "$REPO/samples.tsv" ] || { echo "missing input: $REPO/samples.tsv" >&2; exit 1; }
# Cell-line display names (figure labels / titles) = samples.tsv cell_line, the manuscript's cell-line list
label_of() { awk -F'\t' -v id="$1" 'NR==1{for(i=1;i<=NF;i++) if($i=="cell_line") c=i} NR>1 && $1==id {print $c}' "$REPO/samples.tsv"; }
NGP_LABEL=$(label_of NGP); SKNFI_LABEL=$(label_of SKNrI)
[ -n "$NGP_LABEL" ] && [ -n "$SKNFI_LABEL" ] || { echo "cell_line lookup failed in $REPO/samples.tsv (NGP, SKNrI)" >&2; exit 1; }
for tool in "$PYTHON" "$RSCRIPT" "$JUPYTER"; do command -v "$tool" >/dev/null || { echo "not found: $tool" >&2; exit 1; }; done
echo "repo: $REPO"; echo "out:  $OUT"; "$PYTHON" --version; "$RSCRIPT" --version 2>&1 | head -1

fail() { echo "FAILED: $1 (see $2)" >&2; tail -20 "$2" >&2; exit 1; }

echo "== consistency check: master.csv ecDNA columns vs AmpliconClassifier"
"$PYTHON" "$REPO/02_ecDNA/04_check_ecDNA_columns.py" --master "$DATA/master.csv" --samples "$REPO/samples.tsv" --dir "$DATA/ampliconclassifier" > "$LOG/check_ecDNA.log" 2>&1 || fail "ecDNA column check" "$LOG/check_ecDNA.log"
tail -1 "$LOG/check_ecDNA.log"
echo "== consistency check: IDs / MYCN status / copy number / FPKM across the figure tables"
"$PYTHON" "$HERE/01_check_figure_tables.py" --samples "$REPO/samples.tsv" --data "$DATA" > "$LOG/check_figure_tables.log" 2>&1 || fail "figure-table cross-check" "$LOG/check_figure_tables.log"
tail -1 "$LOG/check_figure_tables.log"

echo "== Fig 1A oncoprint (notebook, kernel $JUPYTER_KERNEL)"
rm -f "$HERE/neuroblastoma_heatmap.pdf" "$HERE/neuroblastoma_heatmap.png"
( cd "$HERE" && "$JUPYTER" nbconvert --to notebook --execute fig1A_oncoprint.ipynb \
    --output "$LOG/fig1A_oncoprint_executed.ipynb" --ExecutePreprocessor.kernel_name="$JUPYTER_KERNEL" \
    --ExecutePreprocessor.timeout=600 ) > "$LOG/fig1A_nbconvert.log" 2>&1 || fail "Fig 1A notebook" "$LOG/fig1A_nbconvert.log"
[ -s "$HERE/neuroblastoma_heatmap.pdf" ] || fail "Fig 1A notebook produced no neuroblastoma_heatmap.pdf" "$LOG/fig1A_nbconvert.log"
mv "$HERE/neuroblastoma_heatmap.pdf" "$OUT/fig1A_oncoprint.pdf"; mv "$HERE/neuroblastoma_heatmap.png" "$OUT/fig1A_oncoprint.png"

echo "== Fig 1B chr2 copy number (NGP ecDNA+, SKNFI ecDNA-)"
"$PYTHON" "$HERE/fig1B_chr2_copy_number.py" --cn "$DATA/chr2_cn/NGP_cn.txt"    --sample NGP   --label "$NGP_LABEL"   --out_dir "$OUT" > "$LOG/fig1B_NGP.log" 2>&1 || fail "Fig 1B NGP" "$LOG/fig1B_NGP.log"
"$PYTHON" "$HERE/fig1B_chr2_copy_number.py" --cn "$DATA/chr2_cn/SKNYrI_cn.txt" --sample SKNFI --label "$SKNFI_LABEL" --ecdna_negative --out_dir "$OUT" > "$LOG/fig1B_SKNFI.log" 2>&1 || fail "Fig 1B SKNFI" "$LOG/fig1B_SKNFI.log"
mv "$OUT/NGP_chr2_cn.png" "$OUT/fig1B_NGP_chr2_cn.png"; mv "$OUT/SKNFI_chr2_cn.png" "$OUT/fig1B_SKNFI_chr2_cn.png"

echo "== Fig 2A / 2B"
"$PYTHON" "$HERE/fig2A_MYCN_CN_vs_FPKM.py" "$DATA/FPKM_CN_MYCN.txt" "$OUT" "$REPO/samples.tsv" > "$LOG/fig2A.log" 2>&1 || fail "Fig 2A" "$LOG/fig2A.log"
mv "$OUT/MYCN_CN_vs_GE.pdf" "$OUT/fig2A_MYCN_CN_vs_FPKM.pdf"
"$PYTHON" "$HERE/fig2B_MYCN_FPKM_boxplot.py" "$DATA/MYCN_FPKM.txt" "$OUT" > "$LOG/fig2B.log" 2>&1 || fail "Fig 2B" "$LOG/fig2B.log"
mv "$OUT/MYCN_FPKM_boxplot.pdf" "$OUT/fig2B_MYCN_FPKM_boxplot.pdf"
grep -E "Mann|p ?=|p-value" "$LOG/fig2B.log" || true

echo "== Fig 2C volcano"
"$PYTHON" "$HERE/fig2C_volcano.py" "$DATA/deseq2_gsea/deseq2_results_ecDNA.csv" "$OUT/fig2C_volcano.pdf" > "$LOG/fig2C.log" 2>&1 || fail "Fig 2C" "$LOG/fig2C.log"
grep -E "^Up:" "$LOG/fig2C.log" || fail "Fig 2C printed no Up/Down counts" "$LOG/fig2C.log"

echo "== Fig 2D / 2E (fgsea + RS panel)"
GSEA_OUT="$OUT/gsea_rerun"; rm -rf "$GSEA_OUT"; mkdir -p "$GSEA_OUT"       # never reuse a previous run's files
"$RSCRIPT" "$REPO/05_differential_expression/02_gsea_and_RS_genes.R" "$DATA/deseq2_gsea/deseq2_results_ecDNA.csv" "$GSEA_OUT" > "$LOG/gsea.log" 2>&1 || fail "02_gsea_and_RS_genes.R" "$LOG/gsea.log"
grep -E "Hallmark|RS genes" "$LOG/gsea.log" || true
[ -s "$GSEA_OUT/GSEA_hallmark_dotplot.pdf" ] && [ -s "$GSEA_OUT/RS_genes_log2FC.pdf" ] || fail "GSEA script produced no plots" "$LOG/gsea.log"
cp "$GSEA_OUT/GSEA_hallmark_dotplot.pdf" "$OUT/fig2D_hallmark_dotplot.pdf"; cp "$GSEA_OUT/RS_genes_log2FC.pdf" "$OUT/fig2E_RS_genes_log2FC.pdf"
if cmp -s "$GSEA_OUT/hallmark_GSEA_results.csv" "$DATA/deseq2_gsea/hallmark_GSEA_results.csv"; then
  echo "hallmark_GSEA_results.csv: identical to frozen data/deseq2_gsea copy"
else
  echo "FAILED: re-computed hallmark_GSEA_results.csv differs from data/deseq2_gsea/hallmark_GSEA_results.csv (check fgsea / msigdbr versions in $GSEA_OUT/gsea_session_versions.txt)" >&2; exit 1
fi

echo "== done: $(ls "$OUT" | grep -c '^fig.*\.\(pdf\|png\)$') panel files in $OUT (logs and executed notebook in $LOG)"; ls -1 "$OUT" | grep '^fig'

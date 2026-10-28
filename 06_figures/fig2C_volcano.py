"""Figure 2C — volcano plot of DESeq2 results (ecDNA+ vs ecDNA−).

Input : deseq2_results_ecDNA.csv from 05_differential_expression/01_deseq2_ecDNA.R
Output: volcano_ecDNA.pdf

Genes with NA padj (removed by DESeq2 independent filtering) are not plotted, so the
legend counts (n up / down / not significant) sum to the number of genes tested, not
to the number of rows in the CSV. Fragment-level counts (Oct 2026 re-count) give
486 up / 404 down / 16,245 not significant.
"""
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import adjustText
from adjustText import adjust_text

# ── Paths: edit here or pass on the command line ─────────────────────────────
# Usage: python fig2C_volcano.py [deseq2_results_ecDNA.csv] [volcano_ecDNA.pdf]
DESEQ2_CSV = sys.argv[1] if len(sys.argv) > 1 else '/path/to/deseq2/deseq2_results_ecDNA.csv'
OUT_PDF    = sys.argv[2] if len(sys.argv) > 2 else '/path/to/figures/volcano_ecDNA.pdf'
if DESEQ2_CSV.startswith('/path/to/') or OUT_PDF.startswith('/path/to/'):
    sys.exit('usage: python fig2C_volcano.py <deseq2_results_ecDNA.csv> <volcano_ecDNA.pdf>  (or edit the defaults)')
# ─────────────────────────────────────────────────────────────────────────────

print(f"python {sys.version.split()[0]} | pandas {pd.__version__} | matplotlib {matplotlib.__version__} "
      f"| adjustText {getattr(adjustText, '__version__', 'unknown')}")

# ── Data ──────────────────────────────────────────────────────────────────────
df = pd.read_csv(DESEQ2_CSV)
df = df.rename(columns={'Unnamed: 0': 'ensembl'})
n_total = len(df)
df = df.dropna(subset=['padj', 'log2FoldChange']).reset_index(drop=True)
print(f"Genes in CSV: {n_total}; with padj (plotted): {len(df)}")

# ── Thresholds ────────────────────────────────────────────────────────────────
FC_THR  = 1.0
FDR_THR = 0.05
N_LABEL = 15                      # label the N most significant DE genes (by padj)
LEFT_LABELS = {'TBC1D17'}         # cosmetic: labels placed manually to the left to avoid
                                  # the dense centre; applied only if the gene is in the top N

# ── Significance category ─────────────────────────────────────────────────────
def categorise(row):
    sig  = row['padj'] < FDR_THR
    up   = row['log2FoldChange'] >  FC_THR
    down = row['log2FoldChange'] < -FC_THR
    if sig and up:   return 'up'
    elif sig and down: return 'down'
    else:             return 'ns'

df['category']       = df.apply(categorise, axis=1)
df['neg_log10_padj'] = -np.log10(df['padj'])

# ── Gene selections ───────────────────────────────────────────────────────────
sig_df = df[df['category'] != 'ns'].copy()
top_n  = sig_df.nsmallest(N_LABEL, 'padj')['Symbol'].tolist()

# ── Style ─────────────────────────────────────────────────────────────────────
BG       = 'white'
C_UP     = '#C0392B'
C_DOWN   = '#2471A3'
C_NS     = '#D5D8DC'
C_SPINE  = '#BDC3C7'
FONT_COL = '#2C3E50'
FONT_SCALE = 1.4  # uniformly enlarge every figure font by 40%

# ── Figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 8))
fig.patch.set_facecolor(BG)
ax.set_facecolor(BG)

# Y axis
y_cap = df['neg_log10_padj'].max() * 1.08
ax.set_ylim(-0.3, y_cap)

# ── Scatter layers ────────────────────────────────────────────────────────────
for cat, color, z, alpha in [
    ('ns',   C_NS,   1, 0.40),
    ('down', C_DOWN, 2, 0.80),
    ('up',   C_UP,   3, 0.80),
]:
    sub = df[df['category'] == cat]
    ax.scatter(sub['log2FoldChange'], sub['neg_log10_padj'],
               c=color, s=6, alpha=alpha, linewidths=0, zorder=z)

# ── Threshold lines ───────────────────────────────────────────────────────────
ax.axhline(-np.log10(FDR_THR), color=C_SPINE, lw=0.8, ls='--', zorder=0)
ax.axvline( FC_THR,             color=C_SPINE, lw=0.8, ls='--', zorder=0)
ax.axvline(-FC_THR,             color=C_SPINE, lw=0.8, ls='--', zorder=0)

# ── Top-N labels using adjustText with tight bounds ──────────────────────────
label_df = df[df['Symbol'].isin(top_n)].copy()

# Points to avoid (all plotted points) — pass to adjustText
all_x = df['log2FoldChange'].values
all_y = df['neg_log10_padj'].values

texts = []
point_coords = []

for _, row in label_df.iterrows():
    x = row['log2FoldChange']
    y = row['neg_log10_padj']
    c = C_DOWN if row['log2FoldChange'] < 0 else C_UP

    if row['Symbol'] in LEFT_LABELS:
        ax.annotate(row['Symbol'],
            xy     = (x, y),
            xytext = (x - 5.0, y + 0.5),
            fontsize=10 * FONT_SCALE, color=c,
            ha='right', va='center', zorder=9,
            arrowprops=dict(arrowstyle='-', color='#AAAAAA',
                            lw=0.5, shrinkA=2, shrinkB=2))
        continue

    t = ax.text(x, y, row['Symbol'],
                fontsize=10 * FONT_SCALE, color=c,
                fontweight='normal', zorder=7,
                ha='center', va='center')
    texts.append(t)
    point_coords.append((x, y))

# adjustText — keep labels close to points, avoid all scatter points
adjust_text(
    texts,
    x          = all_x,
    y          = all_y,
    ax         = ax,
    expand_text   = (1.15, 1.3),
    expand_points = (1.2,  1.4),
    force_text    = (0.15, 0.25),
    force_points  = (0.10, 0.15),
    lim           = 300,
    only_move     = {'points': 'xy', 'texts': 'xy'},
    arrowprops    = dict(arrowstyle='-', color='#AAAAAA',
                         lw=0.5, shrinkA=2, shrinkB=2),
)

# ── Axes ──────────────────────────────────────────────────────────────────────
ax.set_xlabel('log₂ fold change  (ecDNA+ / ecDNA−)',
              fontsize=11 * FONT_SCALE, color=FONT_COL, labelpad=10)
ax.set_ylabel('−log₁₀ (adjusted p-value)',
              fontsize=11 * FONT_SCALE, color=FONT_COL, labelpad=10)
ax.spines[['top','right']].set_visible(False)
ax.spines[['left','bottom']].set_color(C_SPINE)
ax.tick_params(axis='both', colors='#555', labelsize=12 * FONT_SCALE, length=4)
ax.yaxis.grid(False)
ax.xaxis.grid(False)

# ── Legend ────────────────────────────────────────────────────────────────────
n_up   = len(df[df['category'] == 'up'])
n_down = len(df[df['category'] == 'down'])
n_ns   = len(df[df['category'] == 'ns'])

patches = [
    mpatches.Patch(facecolor=C_UP,   edgecolor='none', label=f'Higher in ecDNA+  (n={n_up})'),
    mpatches.Patch(facecolor=C_DOWN, edgecolor='none', label=f'Higher in ecDNA−  (n={n_down})'),
    mpatches.Patch(facecolor=C_NS,   edgecolor='none', label=f'Not significant  (n={n_ns})'),
]
ax.legend(handles=patches, fontsize=10 * FONT_SCALE, frameon=False,
          loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.06))

plt.tight_layout(pad=1.5)
fig.savefig(OUT_PDF, dpi=900, bbox_inches='tight')
plt.close(fig)
print(f"Saved: {OUT_PDF}")
print(f"Up: {n_up}  |  Down: {n_down}  |  NS: {n_ns}")
print(f"Top {N_LABEL} labelled: {top_n}")

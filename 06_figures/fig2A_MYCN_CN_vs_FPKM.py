import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── Data ──────────────────────────────────────────────────────────────────────
df = pd.read_csv('./FPKM_CN_MYCN.txt', sep='\t')   # manually assembled table, see data/README.md

# ── Style tokens ──────────────────────────────────────────────────────────────
BG       = 'white'
C_TRUE   = '#C0392B'
C_FALSE  = '#AEB6BF'
C_SPINE  = '#BDC3C7'
FONT_COL = '#2C3E50'

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 6))
fig.patch.set_facecolor(BG)
ax.set_facecolor(BG)

for amp, color, label, zorder, size in [
    (False, C_FALSE, 'ecDNA− (MYCN non-amp)', 3, 13),
    (True,  C_TRUE,  'ecDNA+ (MYCN amp)',     4, 35),
]:
    sub = df[df['MYCNamp'] == amp]
    ax.scatter(sub['GE'], sub['CN'],
               color=color, s=size, zorder=zorder,
               edgecolors='white', linewidths=0.6,
               label=label)

# ── Sample labels — amplified only ───────────────────────────────────────────
offsets = {
    'CHP134': (  8,   6),
    'IMR32':  (  8,  -8),
    'CHP212': (  8,  -6),
    'KAN':    (  8,   6),
    'KANR':   (  8,  -6),
    'KPNYN':  (  8,  -6),
    'Kelly':  (  8,   6),
    'NB1643': (  8,  -8),
    'NB10':   (  8,   6),
    'NGP':    (  8,   6),
    'TR14':   (  8,   6),
    'KCNR':   (  8,   6),
    'LAN1':   (  8,   6),
    'LANS':   (  8,   6),
    'SKNBE2': (  8,   6),
}

for _, row in df.iterrows():
    if not row['MYCNamp']:   # skip non-amplified
        continue
    dx, dy = offsets.get(row['Sample'], (8, 6))
    ax.annotate(
        row['Sample'],
        xy=(row['GE'], row['CN']),
        xytext=(row['GE'] + dx, row['CN'] + dy),
        fontsize=8.5, color=FONT_COL, va='center',
        arrowprops=dict(arrowstyle='-', color=C_SPINE, lw=0.5)
        if abs(dx) > 10 or abs(dy) > 8 else None,
    )

# ── Axes styling ──────────────────────────────────────────────────────────────
ax.set_xlabel('MYCN gene expression (FPKM)', fontsize=11, color=FONT_COL, labelpad=10)
ax.set_ylabel('MYCN copy number',            fontsize=11, color=FONT_COL, labelpad=10)
ax.yaxis.grid(False)
ax.xaxis.grid(False)
ax.spines[['top', 'right']].set_visible(False)
ax.spines[['left', 'bottom']].set_color(C_SPINE)
ax.tick_params(axis='both', colors='#555', labelsize=9)
ax.tick_params(axis='x', length=4)
ax.tick_params(axis='y', length=4)

# ── Legend: unboxed, horizontal, top ─────────────────────────────────────────
legend_patches = [
    mpatches.Patch(facecolor=C_TRUE,  edgecolor='none', label='ecDNA+ (MYCN amp)'),
    mpatches.Patch(facecolor=C_FALSE, edgecolor='none', label='ecDNA− (MYCN non-amp)'),
]
ax.legend(handles=legend_patches, fontsize=10.5, frameon=False,
          loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.06))

plt.tight_layout(pad=1.5)
out_dir = '.'
fig.savefig(f'{out_dir}/MYCN_CN_vs_GE.pdf', dpi=900, bbox_inches='tight')
plt.close(fig)
print("Saved: MYCN_CN_vs_GE.pdf")
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy import stats

# ── Data ──────────────────────────────────────────────────────────────────────
df = pd.read_csv('./MYCN_FPKM.txt', sep='\t')   # manually assembled table, see data/README.md

# ── Style tokens ──────────────────────────────────────────────────────────────
BG       = 'white'
C_TRUE   = '#C0392B'
C_FALSE  = '#AEB6BF'
C_SPINE  = '#BDC3C7'
FONT_COL = '#2C3E50'

# ── Groups ────────────────────────────────────────────────────────────────────
amp     = df[df['MYCN_amp'] == True]['MYCN'].values
non_amp = df[df['MYCN_amp'] == False]['MYCN'].values

# ── Mann-Whitney U test ───────────────────────────────────────────────────────
stat, pval = stats.mannwhitneyu(amp, non_amp, alternative='greater')
if pval < 0.001:
    pval_str = '***'
elif pval < 0.01:
    pval_str = '**'
elif pval < 0.05:
    pval_str = '*'
else:
    pval_str = f'p = {pval:.3f}'

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5, 6))
fig.patch.set_facecolor(BG)
ax.set_facecolor(BG)

positions = [0, 1]
data      = [non_amp, amp]
colors    = [C_FALSE, C_TRUE]
labels    = ['ecDNA−\n(MYCN non-amp)', 'ecDNA+\n(MYCN amp)']

# Box plots
bp = ax.boxplot(data, positions=positions, widths=0.45,
                patch_artist=True, showfliers=False,
                medianprops=dict(color='white', linewidth=2),
                whiskerprops=dict(color=C_SPINE, linewidth=1.2),
                capprops=dict(color=C_SPINE, linewidth=1.2),
                boxprops=dict(linewidth=0))

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.5)

# Jittered dots overlaid
np.random.seed(42)
for xi, vals, color in zip(positions, data, colors):
    jitter = np.random.uniform(-0.08, 0.08, len(vals))
    ax.scatter(xi + jitter, vals, color=color, s=30, zorder=3,
               edgecolors='white', linewidths=0.5, alpha=0.9)

# ── p-value bracket ───────────────────────────────────────────────────────────
y_max = max(np.max(amp), np.max(non_amp))
y_br  = y_max * 1.08
y_top = y_max * 1.22
ax.set_ylim(-10, y_top)
ax.plot([0, 0, 1, 1], [y_br * 0.97, y_br, y_br, y_br * 0.97],
        color=FONT_COL, linewidth=1.2, clip_on=False)
ax.text(0.5, y_br * 1.01, pval_str, ha='center', va='bottom',
        fontsize=10, color=FONT_COL, fontweight='normal', clip_on=False)

# ── Axes styling ──────────────────────────────────────────────────────────────
ax.set_xticks(positions)
ax.set_xticklabels(labels, fontsize=11, fontweight='normal', color=FONT_COL)
ax.set_ylabel('MYCN expression (FPKM)', fontsize=11, color=FONT_COL, labelpad=10)
ax.yaxis.grid(False)
ax.xaxis.grid(False)
ax.spines[['top', 'right']].set_visible(False)
ax.spines[['left', 'bottom']].set_color(C_SPINE)
ax.tick_params(axis='y', colors='#555', labelsize=9)
ax.tick_params(axis='x', length=0)

# ── Legend: unboxed, horizontal, top ─────────────────────────────────────────
legend_patches = [
    mpatches.Patch(facecolor=C_TRUE,  edgecolor='none', label='ecDNA+ (MYCN amp)'),
    mpatches.Patch(facecolor=C_FALSE, edgecolor='none', label='ecDNA− (MYCN non-amp)'),
]
ax.legend(handles=legend_patches, fontsize=10.5, frameon=False,
          loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.06))

plt.tight_layout(pad=1.5)
out_dir = '.'
fig.savefig(f'{out_dir}/MYCN_FPKM_boxplot.pdf', dpi=900, bbox_inches='tight')
plt.close(fig)
print("Saved: MYCN_FPKM_boxplot.pdf")
#!/usr/bin/env python3
"""
Chr2 Copy Number Plot with MYCN highlight
==========================================
Plots CN segments for chromosome 2 only, in the step-line style.
MYCN locus (hg38: chr2:15,940,587-15,946,097) is highlighted with a
shaded band and label.

Usage:
    python plot_chr2_cn.py --cn CHP134_cn.txt --sample CHP134
    python CN.py --all --data_dir /path/to/folder
"""

import argparse
import glob
import os
import traceback

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter

# ── MYCN locus (hg19) ────────────────────────────────────────────────────────
MYCN_START = 15_940_587
MYCN_END   = 15_946_097

# ── Colors (from CT_ecDNA_plot.py) ───────────────────────────────────────────
COL = {
    "cn_seg":     "#202121",   # all CN segments — same as original cn_normal/cn_amp
    "diploid":    "#AAAAAA",   # dashed CN=2 reference line
    "mycn_band":  "#E64B35",   # MYCN highlight band (ecDNA red)
    "mycn_label": "#C0392B",
    "chrom_bg":   "#F7F7F7",
    "amp_label":  "#202121",
}

CN_CAP           = 15
DEFAULT_DATA_DIR = "."   # folder with <sample>_cn.txt files (03_chromothripsis/04_format_cn_for_shatterseek.py)


# ══════════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ══════════════════════════════════════════════════════════════════════════════

def load_cn(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={"chromosome": str})
    df["chromosome"] = df["chromosome"].str.lstrip("chr")
    df["CN"] = pd.to_numeric(df["CN"], errors="coerce").fillna(2)
    return df


# ══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def mb_fmt(x, pos=None):
    return f"{x/1e6:.0f} Mb"


def merge_amp_labels(segs, min_gap):
    """Merge nearby amplified segments so CN labels don't overlap."""
    if not segs:
        return []
    segs = sorted(segs, key=lambda x: x[0])
    merged, (cs, ce, ccn) = [], segs[0]
    for s, e, cn in segs[1:]:
        if s - ce < min_gap:
            ce  = max(ce, e)
            ccn = max(ccn, cn)
        else:
            merged.append((cs, ce, ccn))
            cs, ce, ccn = s, e, cn
    merged.append((cs, ce, ccn))
    return merged


# ══════════════════════════════════════════════════════════════════════════════
# PLOT
# ══════════════════════════════════════════════════════════════════════════════

def plot_chr2(cn_df: pd.DataFrame,
              sample: str,
              cn_cap: int          = CN_CAP,
              out_path: str        = None,
              ecdna_positive: bool = True):

    seg = cn_df[cn_df["chromosome"] == "2"].sort_values("start").copy()
    if seg.empty:
        print(f"  SKIP {sample}: no chr2 segments found.")
        return

    win_start = 0
    win_end   = int(seg["end"].max())

    # ── Figure layout ────────────────────────────────────────────────────────
    #   Row 0 (thin):  ideogram / MYCN marker strip
    #   Row 1 (tall):  CN track
    FIG_W, FIG_H = 9.3, 3.2
    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor="white")

    gs = fig.add_gridspec(
        2, 1,
        height_ratios=[0.12, 0.88],
        hspace=0.04,
        left=0.06, right=0.97,
        top=0.82, bottom=0.18,
    )
    ax_ideo = fig.add_subplot(gs[0])
    ax_cn   = fig.add_subplot(gs[1], sharex=ax_ideo)

    # ── Ideogram strip ───────────────────────────────────────────────────────
    ax_ideo.set_xlim(win_start, win_end)
    ax_ideo.set_ylim(0, 1)
    ax_ideo.add_patch(mpatches.FancyBboxPatch(
        (win_start, 0.1), win_end - win_start, 0.8,
        boxstyle="round,pad=0.0", linewidth=0.6,
        edgecolor="#888888", facecolor="#D9D9D9"
    ))
    # MYCN band on ideogram
    ax_ideo.add_patch(mpatches.Rectangle(
        (MYCN_START, 0.05), MYCN_END - MYCN_START, 0.9,
        linewidth=0, facecolor=COL["mycn_band"], alpha=0.9, zorder=3
    ))
    ax_ideo.text(
        (MYCN_START + MYCN_END) / 2, 1.15, "MYCN",
        ha="center", va="bottom", fontsize=7.5,
        fontweight="bold", color=COL["mycn_label"],
        transform=ax_ideo.transData, clip_on=False
    )
    ax_ideo.set_yticks([])
    ax_ideo.set_xticks([])
    for sp in ax_ideo.spines.values():
        sp.set_visible(False)
    ax_ideo.text(
        -0.005, 0.5, "Chr 2",
        transform=ax_ideo.transAxes,
        ha="right", va="center",
        fontsize=8, fontweight="bold", color="#333333"
    )

    # ── CN track ─────────────────────────────────────────────────────────────
    ax_cn.set_xlim(win_start, win_end)
    ax_cn.set_facecolor(COL["chrom_bg"])

    # MYCN shaded region
    ax_cn.axvspan(MYCN_START, MYCN_END,
                  color=COL["mycn_band"], alpha=0.12, zorder=0, linewidth=0)
    ax_cn.axvline(MYCN_START, color=COL["mycn_band"],
                  linewidth=0.8, linestyle="--", alpha=0.5, zorder=1)
    ax_cn.axvline(MYCN_END,   color=COL["mycn_band"],
                  linewidth=0.8, linestyle="--", alpha=0.5, zorder=1)

    # diploid reference
    ax_cn.axhline(2, color=COL["diploid"], linewidth=0.7,
                  linestyle="--", alpha=0.8, zorder=2)

    # CN segments
    amp_segs = []
    for _, row in seg.iterrows():
        raw_cn = row["CN"]
        y      = min(raw_cn, cn_cap)
        lw     = 3.0 if raw_cn >= cn_cap else 2.5
        ax_cn.plot([row["start"], row["end"]], [y, y],
                   color=COL["cn_seg"], linewidth=lw,
                   solid_capstyle="butt", alpha=0.95, zorder=4)
        if raw_cn >= cn_cap:
            amp_segs.append((row["start"], row["end"], raw_cn))

    # CN labels for capped (amplified) segments
    min_gap = (win_end - win_start) * 0.02
    for ms, me, mcn in merge_amp_labels(amp_segs, min_gap):
        ax_cn.text(
            (ms + me) / 2, cn_cap + 0.35, f"CN≈{int(mcn)}",
            ha="center", va="bottom", fontsize=6,
            color=COL["amp_label"], fontweight="bold", clip_on=True
        )

    # cap dashed line
    ax_cn.axhline(cn_cap, color="#BBBBBB", linewidth=0.4,
                  linestyle=":", alpha=0.5, zorder=1)

    # y-axis
    y_max   = max(min(seg["CN"].max(), cn_cap) + 2, 8)
    yticks  = [t for t in [2, 5, 10, cn_cap] if t <= y_max]
    ax_cn.set_ylim(0, y_max)
    ax_cn.set_yticks(yticks)
    ax_cn.set_yticklabels(
        [str(t) if t < cn_cap else f"≥{cn_cap}" for t in yticks],
        fontsize=7
    )
    ax_cn.set_ylabel("Copy Number", fontsize=8, labelpad=4, color="#444444")

    # x-axis
    ax_cn.xaxis.set_major_formatter(FuncFormatter(mb_fmt))
    ax_cn.tick_params(axis="x", labelsize=7, length=3)
    ax_cn.set_xlabel("Chr2 position", fontsize=8, labelpad=3)

    for sp in ["top", "right"]:
        ax_cn.spines[sp].set_visible(False)
    ax_cn.spines["left"].set_color("#AAAAAA")
    ax_cn.spines["left"].set_linewidth(0.6)
    ax_cn.spines["bottom"].set_color("#AAAAAA")
    ax_cn.spines["bottom"].set_linewidth(0.6)
    ax_cn.tick_params(axis="y", length=2, labelcolor="#444444")

    # MYCN x-axis label
    mycn_mid = (MYCN_START + MYCN_END) / 2
    ax_cn.annotate(
        "", xy=(mycn_mid, -0.08), xycoords=("data", "axes fraction"),
        xytext=(mycn_mid, -0.02), textcoords=("data", "axes fraction"),
        arrowprops=dict(arrowstyle="-", color=COL["mycn_band"],
                        lw=1.2, alpha=0.7),
        annotation_clip=False
    )

    # ── Legend ───────────────────────────────────────────────────────────────
    handles = [
        plt.Line2D([0], [0], color=COL["cn_seg"], lw=2.5, label="CN segment"),
        plt.Line2D([0], [0], color=COL["diploid"], lw=1,
                   linestyle="--", label="CN = 2 (diploid)"),
        mpatches.Patch(facecolor=COL["mycn_band"], alpha=0.5,
                       label="MYCN locus"),
    ]
    fig.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(0.97, 0.98),
        bbox_transform=fig.transFigure,
        fontsize=7, frameon=True, framealpha=0.95,
        edgecolor="lightgrey", borderpad=0.6,
        handlelength=1.6, labelspacing=0.3,
    )

    # ── Title ────────────────────────────────────────────────────────────────
    tag = "  [ecDNA+]" if ecdna_positive else "  [ecDNA−]"
    fig.suptitle(
        f"{sample}",
        fontsize=10, fontweight="bold", y=0.97,
        color="#1A1A1A"
    )

    # ── Save ─────────────────────────────────────────────────────────────────
    if out_path is None:
        out_path = f"{sample}_chr2_cn.png"
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=900, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved → {out_path}")
    return out_path


# ══════════════════════════════════════════════════════════════════════════════
# BATCH HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def discover_cn_files(data_dir: str):
    files = sorted(glob.glob(os.path.join(data_dir, "*_cn.txt")))
    results = []
    for f in files:
        base   = os.path.splitext(os.path.basename(f))[0]
        sample = base[:-3] if base.endswith("_cn") else base
        results.append((sample, f))
    return results


def process_sample(sample, cn_path, out_dir, cn_cap, ecdna_positive):
    print(f"\n{'='*55}\nProcessing: {sample}\n{'='*55}")
    try:
        cn_df    = load_cn(cn_path)
        chr2_seg = cn_df[cn_df["chromosome"] == "2"]
        print(f"  Chr2 segments: {len(chr2_seg)}  |  "
              f"CN range: {chr2_seg['CN'].min():.0f}–{chr2_seg['CN'].max():.0f}")
        out_path = os.path.join(out_dir, f"{sample}_chr2_cn.png")
        plot_chr2(cn_df, sample=sample, cn_cap=cn_cap,
                  out_path=out_path, ecdna_positive=ecdna_positive)
    except Exception as exc:
        print(f"  ERROR: {exc}")
        traceback.print_exc()


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Chr2 CN step-line plot with MYCN highlight"
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--all", action="store_true",
                      help="Batch-process ALL *_cn.txt in --data_dir")
    mode.add_argument("--cn",  default=None,
                      help="Path to a single <sample>_cn.txt")

    parser.add_argument("--data_dir", default=DEFAULT_DATA_DIR,
                        help=f"Folder to scan (default: {DEFAULT_DATA_DIR})")
    parser.add_argument("--out_dir",  default=None,
                        help="Output folder (default: same as --data_dir)")
    parser.add_argument("--sample",   default=None,
                        help="Sample name override (single-file mode only)")
    parser.add_argument("--cn_cap",   type=int, default=CN_CAP,
                        help=f"Cap display CN (default: {CN_CAP})")
    parser.add_argument("--ecdna_positive", action="store_true", default=True)
    parser.add_argument("--ecdna_negative", dest="ecdna_positive",
                        action="store_false")

    args = parser.parse_args()

    if args.all:
        samples = discover_cn_files(args.data_dir)
        if not samples:
            print(f"No *_cn.txt files found in: {args.data_dir}")
            return
        out_dir = args.out_dir or args.data_dir
        os.makedirs(out_dir, exist_ok=True)
        print(f"Found {len(samples)} sample(s) in {args.data_dir}")
        ok = fail = 0
        for sample, cn_path in samples:
            try:
                process_sample(sample, cn_path, out_dir,
                               args.cn_cap, args.ecdna_positive)
                ok += 1
            except Exception:
                fail += 1
        print(f"\nDone.  {ok} succeeded, {fail} failed.")

    elif args.cn:
        sample = args.sample or os.path.splitext(
            os.path.basename(args.cn))[0]
        if sample.endswith("_cn"):
            sample = sample[:-3]
        out_dir = args.out_dir or os.path.dirname(os.path.abspath(args.cn))
        os.makedirs(out_dir, exist_ok=True)
        process_sample(sample, args.cn, out_dir,
                       args.cn_cap, args.ecdna_positive)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
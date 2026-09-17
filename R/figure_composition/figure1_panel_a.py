#!/usr/bin/env python3
"""
Programmatic recreation of Figure 1 Panel A — Study Design Overview.
Matches the colour palette and style of the original Figure1_Diagram.png.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import os

# ── Colours (matched to original Figure1_Diagram.png) ────────────────────
COL_ASTRO    = "#F4A460"
COL_MICRO    = "#5B9BD5"
COL_OLIGO    = "#70AD47"
COL_GT       = "#C5E0B4"
COL_VAL      = "#FFD966"
COL_PSEUDO   = "#E8D5F0"   # lavender — PSEUDOBULK DE
COL_COMPARE  = "#D5C8E8"   # purple — COMPARE
COL_ACCUM    = "#FFF2CC"
COL_MODEL    = "#D6E4F0"
COL_EDGE     = "#888888"
WHITE        = "#FFFFFF"

# Font scale factor — bump everything up
FS = 1.35


def rbox(ax, x, y, w, h, bg=WHITE, edge=COL_EDGE, lw=1.2, rounding=0.015):
    p = FancyBboxPatch((x, y), w, h,
                       boxstyle=f"round,pad=0,rounding_size={rounding}",
                       facecolor=bg, edgecolor=edge, linewidth=lw,
                       transform=ax.transAxes, zorder=2,
                       mutation_aspect=0.8)
    ax.add_patch(p)


def txt(ax, x, y, s, fs=8.5, bold=False, italic=False, ha="center", va="center", zorder=3):
    ax.text(x, y, s, ha=ha, va=va, fontsize=fs * FS,
            fontweight="bold" if bold else "normal",
            fontstyle="italic" if italic else "normal",
            transform=ax.transAxes, zorder=zorder)


def arrow(ax, x0, y0, x1, y1, dashed=False, lw=1.3):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="->", color="#555555", lw=lw,
                                linestyle="dashed" if dashed else "solid"),
                xycoords="axes fraction", textcoords="axes fraction", zorder=1)


def hline(ax, x0, x1, y, lw=1.3):
    ax.plot([x0, x1], [y, y], color="#555555", lw=lw,
            transform=ax.transAxes, zorder=1)


def draw_panel_a(ax, n_runs=10, run_label="runs"):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    cx = 0.50          # centre x for most boxes
    lcx = 0.22         # left branch centre
    rcx = 0.62         # right branch centre

    # ─── 1. PSYCHAD COHORT ───────────────────────────────────────────────
    bw, bh = 0.58, 0.095
    bx = cx - bw/2
    by = 0.870
    rbox(ax, bx, by, bw, bh)
    txt(ax, cx, by + bh - 0.015, "PSYCHAD COHORT", fs=10, bold=True)
    # dots
    dot_y = by + 0.052
    for i, dx in enumerate(np.linspace(bx + 0.05, bx + bw - 0.05, 30)):
        ax.plot(dx, dot_y, "o", color=[COL_ASTRO, COL_MICRO, COL_OLIGO][i % 3],
                markersize=3.5, transform=ax.transAxes, zorder=4)
    # legend
    for i, (c, lab) in enumerate([(COL_ASTRO, "Astro"), (COL_MICRO, "Micro"), (COL_OLIGO, "Oligo")]):
        lx = bx + 0.10 + i * 0.15
        ax.plot(lx, by + 0.028, "s", color=c, markersize=5, transform=ax.transAxes, zorder=4)
        txt(ax, lx + 0.020, by + 0.028, lab, fs=6.5, ha="left")

    # arrow
    arrow(ax, cx, by, cx, 0.825)

    # ─── 2. 50/50 DONOR SPLIT ───────────────────────────────────────────
    s_y = 0.755; s_h = 0.072
    rbox(ax, bx, s_y, bw, s_h)
    txt(ax, cx, s_y + s_h - 0.014, "50 / 50 DONOR SPLIT", fs=10, bold=True)
    # sub-boxes
    sw, sh = 0.21, 0.020
    rbox(ax, bx + 0.03, s_y + 0.028, sw, sh, bg=COL_GT, edge=COL_GT, lw=0.5, rounding=0.007)
    txt(ax, bx + 0.03 + sw/2, s_y + 0.028 + sh/2, "Ground truth", fs=7)
    rbox(ax, bx + bw - sw - 0.03, s_y + 0.028, sw, sh, bg=COL_VAL, edge=COL_VAL, lw=0.5, rounding=0.007)
    txt(ax, bx + bw - sw/2 - 0.03, s_y + 0.028 + sh/2, "Validation", fs=7)
    # runs annotation INSIDE the box
    txt(ax, cx, s_y + 0.010, f"\u00d7 {n_runs} {run_label}, sex-stratified", fs=7, italic=True)

    # ─── Fork ────────────────────────────────────────────────────────────
    fork_y = 0.715
    arrow(ax, cx, s_y, cx, fork_y + 0.005)
    hline(ax, lcx, rcx, fork_y)
    arrow(ax, lcx, fork_y, lcx, fork_y - 0.030)
    arrow(ax, rcx, fork_y, rcx, fork_y - 0.030)

    # ─── 3. FULL DE ANALYSIS (left) ─────────────────────────────────────
    lw_ = 0.32; lx = lcx - lw_/2
    ly = 0.610; lh = 0.070
    rbox(ax, lx, ly, lw_, lh, bg=COL_PSEUDO, edge="#9B72B0")
    txt(ax, lcx, ly + lh - 0.015, "FULL DE ANALYSIS", fs=9, bold=True)
    txt(ax, lcx, ly + lh - 0.035, "DESeq2, edgeR & limma-voom", fs=7.5)
    txt(ax, lcx, ly + lh - 0.058, "\u2192 \u22650.2 |LFC| threshold", fs=7.5, italic=True)

    # ─── 4. GRID DOWNSAMPLING (right) ───────────────────────────────────
    rw = 0.44; rx = rcx - rw/2
    ry = 0.525; rh = 0.160
    rbox(ax, rx, ry, rw, rh)
    txt(ax, rcx, ry + rh - 0.015, "GRID DOWNSAMPLING", fs=9, bold=True)
    # mode sub-boxes — evenly distributed down the box
    mw, mh = 0.38, 0.019
    mx = rcx - mw/2
    my1 = ry + rh - 0.052
    rbox(ax, mx, my1, mw, mh, bg="#E2EFDA", edge="#A9D18E", lw=0.8, rounding=0.007)
    txt(ax, rcx, my1 + mh/2, "Cells mode", fs=7.5, bold=True)
    txt(ax, rcx, my1 - 0.011, "12 donor counts \u00d7 14 cell counts", fs=6.5)

    my2 = my1 - 0.046
    rbox(ax, mx, my2, mw, mh, bg="#D6E4F0", edge="#8DB4E2", lw=0.8, rounding=0.007)
    txt(ax, rcx, my2 + mh/2, "UMI mode", fs=7.5, bold=True)
    txt(ax, rcx, my2 - 0.011, "12 donor counts \u00d7 13 UMI budgets", fs=6.5)

    txt(ax, rcx, ry + 0.020, "\u00d7 20 iterations per completed grid point", fs=6.5, italic=True)

    # arrow — longer gap before PSEUDOBULK DE
    arrow(ax, rcx, ry, rcx, ry - 0.045)

    # ─── 5. PSEUDOBULK DE ───────────────────────────────────────────────
    py_ = 0.400; ph = 0.070; pw = 0.38; px = rcx - pw/2
    rbox(ax, px, py_, pw, ph, bg=COL_PSEUDO, edge="#9B72B0")
    txt(ax, rcx, py_ + ph - 0.015, "PSEUDOBULK DE", fs=9, bold=True)
    txt(ax, rcx, py_ + ph - 0.035, "DESeq2, edgeR & limma-voom", fs=7)
    txt(ax, rcx, py_ + ph - 0.052, "at each completed grid point", fs=7)

    # arrow — longer gap before COMPARE
    arrow(ax, rcx, py_, rcx, py_ - 0.040)

    # ─── 6. COMPARE ─────────────────────────────────────────────────────
    cy_ = 0.325; ch = 0.032; cw = 0.36; cxb = rcx - cw/2
    rbox(ax, cxb, cy_, cw, ch, bg=COL_ACCUM, edge="#D4A849")
    txt(ax, rcx, cy_ + ch/2, "Per-gene: detected (1) or missed (0)", fs=7.5)

    # dashed arrow from Full DE down + across to Compare
    arrow(ax, lcx, ly, lcx, cy_ + ch/2, dashed=True)
    arrow(ax, lcx, cy_ + ch/2, cxb, cy_ + ch/2, dashed=True)
    txt(ax, lcx - 0.025, 0.50, "reference DEGs", fs=7, italic=True, va="center")
    # rotate
    ax.texts[-1].set_rotation(90)

    # arrows down from compare to accumulate
    arrow(ax, cx, cy_, cx, cy_ - 0.035)

    # ─── 7. ACCUMULATE ~22M OBSERVATIONS ────────────────────────────────
    aw = 0.68; ah = 0.145; ax_ = cx - aw/2; ay = 0.150
    rbox(ax, ax_, ay, aw, ah, bg=COL_ACCUM, edge="#D4A849")
    txt(ax, cx, ay + ah - 0.015, "ACCUMULATE OBSERVATIONS", fs=9, bold=True)
    # sub-boxes — taller with more spread
    asw, ash = 0.28, 0.080
    asy = ay + 0.030
    rbox(ax, ax_ + 0.03, asy, asw, ash, bg="#FFF9E6", edge="#E6C94F", lw=0.8)
    txt(ax, ax_ + 0.03 + asw/2, asy + ash - 0.015, "Experiment design", fs=7.5, bold=True)
    txt(ax, ax_ + 0.03 + asw/2, asy + ash - 0.044, r"$\mathit{donors},\ \mathit{mean\ cells},\ \mathit{mean\ UMIs}$", fs=6.5)
    txt(ax, ax_ + 0.03 + asw/2, asy + ash - 0.067, "(3 parameters)", fs=6.5)

    rbox(ax, ax_ + aw - asw - 0.03, asy, asw, ash, bg="#FFF9E6", edge="#E6C94F", lw=0.8)
    txt(ax, ax_ + aw - asw/2 - 0.03, asy + ash - 0.015, "Gene characteristics", fs=7.5, bold=True)
    txt(ax, ax_ + aw - asw/2 - 0.03, asy + ash - 0.044, r"$|\mathit{LFC}|,\ \mathit{mean\ expression},\ \mathit{dispersion}$", fs=6.5)
    txt(ax, ax_ + aw - asw/2 - 0.03, asy + ash - 0.067, "(3 parameters)", fs=6.5)

    txt(ax, cx, ay + 0.011, "\u2192 detected (0 / 1)", fs=7.5, italic=True)

    # arrow
    arrow(ax, cx, ay, cx, ay - 0.030)

    # ─── 8. XGBOOST POWER MODEL ─────────────────────────────────────────
    mw_ = 0.72; mh_ = 0.11; mx_ = cx - mw_/2; my_ = 0.010
    rbox(ax, mx_, my_, mw_, mh_, bg=COL_MODEL, edge="#6C8EBF")
    txt(ax, cx, my_ + mh_ - 0.018, "XGBOOST POWER MODEL", fs=9, bold=True)
    txt(ax, cx, my_ + 0.062, "gradient-boosted trees, monotone-constrained", fs=7.5, italic=True)
    txt(ax, cx, my_ + 0.036, r"$\hat{P} = f(\mathit{donors},\ \mathit{mean\ cells},\ \mathit{mean\ UMIs},\ |\mathit{LFC}|,\ \mathit{mean\ expression},\ \mathit{dispersion}$", fs=7.5)
    txt(ax, cx, my_ + 0.011, r"$\mathit{cell\ type},\ \mathit{DE\ method})$", fs=7.5)


def generate_figures(output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Widened Panel A only (no "A" label — standalone figure)
    fig, ax = plt.subplots(figsize=(9, 7.2))
    draw_panel_a(ax, n_runs=10, run_label="runs")
    fig.tight_layout()
    p = os.path.join(output_dir, "Figure1_PanelA_wide.png")
    fig.savefig(p, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"Saved: {p}")

    return p


if __name__ == "__main__":
    out = os.path.join(os.path.dirname(__file__), "output", "figure1")
    generate_figures(out)

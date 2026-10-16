#!/usr/bin/env python3
# ==============================================================================
# plot_oligo_expr_disp_panels.py -- Figure 6: predicted detectability of
# oligodendrocyte genes by functional class, against mean expression (A) and
# DESeq2 dispersion (B), at |LFC| = 2.0, 1.0 and 0.2.
# ------------------------------------------------------------------------------
# Input : Datalink/Gene_Classification/oligo_gene_power_multiLFC.csv
#         (R/gene_classification/classify_oligo_multi_lfc.R)
# Output: Datalink/Gene_Classification/oligo_paper_expr_disp_panels.{png,pdf}
#         Datalink/Gene_Classification/oligo_gene_group_power_summary.csv
#         (per column and |LFC|: genes, expressed genes, and % of expressed genes
#         in each class; the numbers in the column headers, read by the text)
#
# Columns (x axis), ordered by the median mean expression of their expressed genes:
#   * seven Lambert TF families (DNA-binding domain) plus "Other TF" (all other
#     TFs), purple; TFs not expressed in oligodendrocytes are grey points on a
#     floor below the dotted line and are counted as "n.e." in the label;
#   * twelve PANTHER protein classes of expressed non-TF genes;
#   * LoF-intolerant genes (gnomAD pLI >= 0.9, expressed), blue, always last.
# Points are each gene's position coloured by its DESeq2 detectability class at
# the panel's |LFC| (50 donors, 1,000 cells per donor). At most POINT_CAP genes
# per column are drawn: one fixed random subsample (seed SEED) is used for all
# six panels so every gene sits at the same place throughout. Boxes summarise
# all expressed genes of the column (1.5 IQR whiskers, no outliers). Column
# headers give the % Undetectable | % Detectable among expressed genes.
#
# Reconstructed 2026-10-01 from the lost /tmp/oligo_annot/plot_paper_panels.py
# (figure_expr_disp; diary 2026-09-02/03): same layout, groups, colours and
# numbers. Point positions differ from the original because its subsampling
# seed was not recorded.
# ==============================================================================
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import FixedLocator  # noqa: E402

ROOT = os.environ.get("POWER_ANALYSIS_ROOT") or os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", ".."))
GC_DIR = os.path.join(ROOT, "Datalink", "Gene_Classification")
OUT_STEM = os.path.join(GC_DIR, "oligo_paper_expr_disp_panels")

SEED = 1234
POINT_CAP = 260
LFC_ROWS = [("200", "2.0"), ("100", "1.0"), ("020", "0.2")]

TF_FAMILIES = ["C2H2 ZF", "Homeodomain", "bHLH", "bZIP", "Forkhead",
               "Nuclear receptor", "HMG/Sox"]
PANTHER_CLASSES = {
    "intermediate filament": "Intermediate filament",
    "cytokine": "Cytokine",
    "G-protein coupled receptor": "GPCR",
    "growth factor": "Growth factor",
    "intercellular signal molecule": "Signal molecule",
    "transmembrane signal receptor": "TM signal receptor",
    "ion channel": "Ion channel",
    "RNA splicing factor": "RNA splicing",
    "ribosomal protein": "Ribosomal",
    "RNA metabolism protein": "RNA metabolism",
    "chromatin/chromatin-binding, or -regulatory protein": "Chromatin regulator",
    "ubiquitin-protein ligase": "Ubiquitin ligase",
}

C_UNDET, C_UNSURE, C_DET, C_NE = "#d62728", "#ff7f0e", "#2ca02c", "#9e9e9e"
C_TF, C_LOF = "#6a3d9a", "#1f77b4"
CLASS_COLOUR = {"Undetectable": C_UNDET, "Unsure": C_UNSURE, "Detectable": C_DET}

plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42})


def as_bool(s):
    return s.astype(str).str.lower().isin(["true", "1"])


def build_groups(df):
    """List of dicts: name, kind (tf|class|lof), genes (DataFrame incl. non-expressed TFs)."""
    tf = df[df.is_tf]
    groups = []
    for fam in TF_FAMILIES:
        groups.append({"name": fam, "kind": "tf", "genes": tf[tf.lambert_dbd == fam]})
    groups.append({"name": "Other TF", "kind": "tf", "genes": tf[~tf.lambert_dbd.isin(TF_FAMILIES)]})
    non_tf = df[~df.is_tf & df.expressed]
    for pc, label in PANTHER_CLASSES.items():
        groups.append({"name": label, "kind": "class", "genes": non_tf[non_tf.panther_protein_class == pc]})
    for g in groups:
        e = g["genes"][g["genes"].expressed]
        g["median_expr"] = float(np.median(e.basemean)) if len(e) else -np.inf
    groups.sort(key=lambda g: g["median_expr"])
    lof = df[df.lof_intol & df.expressed]
    groups.append({"name": "LoF-intolerant", "kind": "lof", "genes": lof,
                   "median_expr": float(np.median(lof.basemean))})
    rng = np.random.default_rng(SEED)
    for g in groups:
        genes = g["genes"]
        expr = genes[genes.expressed]
        if len(expr) > POINT_CAP:
            expr = expr.iloc[np.sort(rng.choice(len(expr), POINT_CAP, replace=False))]
        g["shown"] = expr
        g["jitter"] = rng.uniform(-0.25, 0.25, len(expr))
        ne = genes[~genes.expressed]
        g["ne"] = ne
        g["ne_jitter"] = rng.uniform(-0.25, 0.25, len(ne))
    return groups


def group_summary(groups):
    """One row per column and |LFC|: class percentages among expressed genes."""
    rows = []
    for g in groups:
        expr = g["genes"][g["genes"].expressed]
        n = len(expr)
        for tag, label in LFC_ROWS:
            cls = expr[f"class_lfc{tag}"]
            row = {"group": g["name"], "kind": g["kind"], "abs_lfc": float(label),
                   "n_genes": len(g["genes"]), "n_expressed": n}
            for c in ("Undetectable", "Unsure", "Detectable"):
                row[f"pct_{c.lower()}"] = 100 * (cls == c).sum() / n if n else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def draw_panel(ax, groups, value_col, lfc_tag, lfc_label, floor_y, line_y, ylim, clip=None, ticks=None):
    cls_col = f"class_lfc{lfc_tag}"
    for i, g in enumerate(groups):
        if g["kind"] == "tf":
            ax.axvspan(i - 0.5, i + 0.5, color=C_TF, alpha=0.06, lw=0, zorder=0)
        elif g["kind"] == "lof":
            ax.axvspan(i - 0.5, i + 0.5, color=C_LOF, alpha=0.10, lw=0, zorder=0)
        expr = g["genes"][g["genes"].expressed]
        vals = expr[value_col].to_numpy(float)
        vals = vals[np.isfinite(vals) & (vals > 0)]
        if clip is not None:
            vals = np.clip(vals, *clip)
        if len(vals):
            edge = C_LOF if g["kind"] == "lof" else "#555555"
            lw = 2.0 if g["kind"] == "lof" else 1.0
            ax.boxplot([vals], positions=[i], widths=0.6, whis=1.5, showfliers=False,
                       patch_artist=True, zorder=2,
                       boxprops=dict(facecolor="white", edgecolor=edge, linewidth=lw),
                       whiskerprops=dict(color="black", linewidth=1.0),
                       capprops=dict(color="black", linewidth=1.0),
                       medianprops=dict(color="black", linewidth=2.2))
        shown = g["shown"]
        y = shown[value_col].to_numpy(float)
        ok = np.isfinite(y) & (y > 0)
        if clip is not None:
            y = np.clip(y, *clip)
        colours = shown[cls_col].map(CLASS_COLOUR).fillna(C_UNDET).to_numpy()
        ax.scatter(i + g["jitter"][ok], y[ok], s=14, c=colours[ok], alpha=0.65,
                   linewidths=0, zorder=3)
        if len(g["ne"]):
            ax.scatter(i + g["ne_jitter"], np.full(len(g["ne"]), floor_y), s=14, c=C_NE,
                       alpha=0.65, linewidths=0, zorder=3)
        # header: % Undetectable | % Detectable among expressed genes
        n = len(expr)
        if n:
            pu = round(100 * (expr[cls_col] == "Undetectable").sum() / n)
            pd_ = round(100 * (expr[cls_col] == "Detectable").sum() / n)
            tr = ax.get_xaxis_transform()
            ax.text(i - 0.06, 1.02, f"{pu}", color=C_UNDET, ha="right", va="bottom",
                    fontsize=8.5, fontweight="bold", transform=tr)
            ax.text(i, 1.02, "|", color="#777777", ha="center", va="bottom", fontsize=8.5, transform=tr)
            ax.text(i + 0.06, 1.02, f"{pd_}", color=C_DET, ha="left", va="bottom",
                    fontsize=8.5, fontweight="bold", transform=tr)
    ax.axhline(line_y, color="#999999", ls=":", lw=1, zorder=1)
    ax.set_yscale("log")
    ax.set_ylim(*ylim)
    if ticks is not None:
        ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.set_xlim(-0.5, len(groups) - 0.5)
    ax.grid(axis="y", which="major", color="#cccccc", ls=":", lw=0.8, zorder=0)
    ax.tick_params(axis="y", labelsize=12)
    ax.text(0.008, 0.95, f"|LFC| = {lfc_label}", transform=ax.transAxes, ha="left", va="top",
            fontsize=17, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#999999"))


def main():
    df = pd.read_csv(os.path.join(GC_DIR, "oligo_gene_power_multiLFC.csv"))
    df["expressed"] = as_bool(df.expressed_in_oligo)
    df["is_tf"] = as_bool(df.is_tf_lambert)
    df["lof_intol"] = as_bool(df.lof_intol)
    groups = build_groups(df)
    group_summary(groups).to_csv(os.path.join(GC_DIR, "oligo_gene_group_power_summary.csv"),
                                 index=False)

    disp = df.loc[df.expressed, "dispersion"].to_numpy(float)
    disp = disp[np.isfinite(disp) & (disp > 0)]
    d_lo, d_hi = np.percentile(disp, [0.5, 99.7])

    fig = plt.figure(figsize=(16.38, 19.0))
    gs = fig.add_gridspec(7, 1, height_ratios=[1, 1, 1, 0.12, 1, 1, 1])
    fig.subplots_adjust(left=0.105, right=0.99, top=0.955, bottom=0.175, hspace=0.16)
    axes = []
    for block, rows in (("A", [0, 1, 2]), ("B", [4, 5, 6])):
        for k, r in enumerate(rows):
            ax = fig.add_subplot(gs[r])
            axes.append(ax)
            tag, label = LFC_ROWS[k]
            if block == "A":
                draw_panel(ax, groups, "basemean", tag, label, floor_y=3.3e-3, line_y=5e-3,
                           ylim=(1.4e-3, 1.3e5), ticks=[1e-1, 1e1, 1e3, 1e5])
                ax.axhline(15, color="#888888", ls="--", lw=1.2, zorder=1)
                if k == 0:
                    ax.text(len(groups) - 0.55, 15 * 1.15, "mean expression 15 (approx. power floor)",
                            ha="right", va="bottom", fontsize=9, color="#777777")
                ax.set_ylabel("Gene counts\n(log)", fontsize=15)
            else:
                draw_panel(ax, groups, "dispersion", tag, label, floor_y=0.5 * d_lo,
                           line_y=0.8 * d_lo, ylim=(0.42 * d_lo, 1.3 * d_hi), clip=(d_lo, d_hi))
                ax.set_ylabel("DESeq2 dispersion\n(log)", fontsize=15)
            ax.set_xticks(range(len(groups)))
            if k == 0:
                title = "Gene counts" if block == "A" else "Dispersion"
                ax.set_title(title, fontsize=19, fontweight="bold", y=1.075)
                ax.text(-0.075, 1.12, block, transform=ax.transAxes, fontsize=24,
                        fontweight="bold", ha="left", va="bottom")
            ax.set_xticklabels([])

    # x labels on the bottom panel only
    labels = []
    for g in groups:
        n_all = len(g["genes"])
        n_ne = int((~g["genes"].expressed).sum())
        labels.append(f"{g['name']}\n(n={n_all}" + (f", {n_ne} n.e.)" if n_ne else ")"))
    bottom = axes[-1]
    bottom.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor", fontsize=14)
    for tick, g in zip(bottom.get_xticklabels(), groups):
        if g["kind"] == "tf":
            tick.set_color(C_TF); tick.set_fontweight("bold")
        elif g["kind"] == "lof":
            tick.set_color(C_LOF); tick.set_fontweight("bold")

    handles = [
        Line2D([], [], marker="o", ls="", color=C_UNDET, markersize=9, label="Undetectable (power ≤ 0.4)"),
        Line2D([], [], marker="o", ls="", color=C_UNSURE, markersize=9, label="Unsure (0.4–0.8)"),
        Line2D([], [], marker="o", ls="", color=C_DET, markersize=9, label="Detectable (power > 0.8)"),
        Line2D([], [], marker="o", ls="", color=C_NE, markersize=9, label="Not expressed in Oligo"),
        Patch(facecolor=C_TF, alpha=0.3, label="TF family (Lambert DBD)"),
        Patch(facecolor="white", edgecolor="#555555", label="Non-TF class (PANTHER)"),
        Patch(facecolor=C_LOF, alpha=0.3, label="LoF-intolerant (pLI ≥ 0.9)"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.045), ncol=4,
               fontsize=12, frameon=True)

    fig.savefig(OUT_STEM + ".png", dpi=200)
    fig.savefig(OUT_STEM + ".pdf", metadata={"CreationDate": None})
    print("wrote", OUT_STEM + ".{png,pdf}")
    for g in groups:
        print(f"  {g['name']:<22} n={len(g['genes']):>5}  median expr={g['median_expr']:.3g}")


if __name__ == "__main__":
    main()

"""
Supplementary figure composition — cross-method / cross-cell-type comparisons.

Purpose: Tile existing per-(cell_type, method) PNGs from the paper's Datalink
output into supplementary comparison panels. Non-destructive: reads source PNGs
only, never modifies a source figure or any paper draft.

This is a sibling of compose_figures.py (which builds the MAIN-text composites).
This script builds the SUPPLEMENTARY cross-method / cross-cell-type grids.

Two source kinds are tiled (paths derived from the cartesian pipeline):
  surface : <method>/Continuous_Viz_Log/Mean_<Metric>_Total_Continuous_SemiLog.png   (script 07)
  figTS   : <method>/Plots/Figure_TS/Figure_<bin>_<Metric>_<mode>.png                (script 09)
where <method> in {deseq2, edger, limma, all_methods} and <bin> in {T (Tiny 0.2-1.0), S (Small 1.0-2.0)}.

Figures built (grid of panels, one panel per method OR per cell type):
  REQUESTED (Chris):
    supp_oligo_recall_surface_bymethod      Oligo recall  surface, 4 methods (2x2)
    supp_oligo_precision_surface_bymethod    Oligo precision surface, 4 methods (2x2)
    supp_oligo_recall_tiny_bymethod          Oligo Tiny-bin recall vs expression, 4 methods (2x2)
    supp_oligo_precision_tiny_bymethod       Oligo Tiny-bin precision vs expression, 4 methods (2x2)
  SUGGESTED (companions / extra angles):
    supp_oligo_recall_small_bymethod         Oligo Small-bin recall vs expression, 4 methods (2x2)
    supp_oligo_precision_small_bymethod      Oligo Small-bin precision vs expression, 4 methods (2x2)
    supp_allmethods_recall_surface_bycelltype   all_methods recall surface, Astro|Micro|Oligo (1x3)
    supp_allmethods_precision_surface_bycelltype all_methods precision surface, Astro|Micro|Oligo (1x3)

Run:
    python3 R/figure_composition/compose_supplementary.py                 # all figures, cells view
    python3 R/figure_composition/compose_supplementary.py --only tiny      # substring filter on name
    python3 R/figure_composition/compose_supplementary.py --mode reads     # reads view instead of cells
    python3 R/figure_composition/compose_supplementary.py --datalink /path/to/Datalink

Outputs: R/figure_composition/output_supplementary/<name>.png
Publish: copies to Datalink/Composed_Figures/Supplementary/<name>.png
Log:     R/figure_composition/output_supplementary/build_log.txt

Any figure whose source PNGs are not all present is SKIPPED with a warning
(so a partial pipeline run still produces whatever is available).
"""

from __future__ import annotations

import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent.parent  # .../Empirical_Power_Analysis
DEFAULT_DATALINK = PROJECT_ROOT / "Datalink"
OUT_ROOT = HERE / "output_supplementary"

# ---------------------------------------------------------------------------
# Pretty labels
# ---------------------------------------------------------------------------

METHOD_LABEL = {
    "deseq2": "DESeq2",
    "edger": "edgeR",
    "limma": "limma-voom",
    "all_methods": "All methods (pooled)",
}
CELL_LABEL = {
    "Astro": "Astrocytes",
    "Micro": "Microglia",
    "Oligo": "Oligodendrocytes",
}
BIN_LABEL = {"T": "Tiny |LFC| 0.2-1.0", "S": "Small |LFC| 1.0-2.0"}

METHODS_4 = ["deseq2", "edger", "limma", "all_methods"]
CELLTYPES_3 = ["Astro", "Micro", "Oligo"]

# ---------------------------------------------------------------------------
# Layout style (single clean supplementary style)
# ---------------------------------------------------------------------------

PANEL_W = 1600            # each panel scaled to this width (px)
GAP = 28                  # gap between panels
PAD = 40                  # outer padding
CAPTION_H = 66            # per-panel caption strip height
SUPER_H = 96              # super-title strip height
BG = (255, 255, 255)
FG = (0, 0, 0)
CAPTION_BG = (238, 238, 238)


def load_font(size: int, bold: bool = True) -> ImageFont.ImageFont:
    suffix = "Bold" if bold else ""
    candidates = [
        f"/usr/share/fonts/truetype/dejavu/DejaVuSans{'-Bold' if bold else ''}.ttf",
        f"/usr/share/fonts/truetype/liberation/LiberationSans{'-' + suffix if bold else '-Regular'}.ttf",
        "/usr/share/fonts/truetype/freefont/FreeSansBold.ttf",
    ]
    for c in candidates:
        if Path(c).exists():
            return ImageFont.truetype(c, size=size)
    return ImageFont.load_default()


CAPTION_FONT = load_font(40, bold=True)
SUPER_FONT = load_font(52, bold=True)


# ---------------------------------------------------------------------------
# Image helpers
# ---------------------------------------------------------------------------

def load_rgb(path: Path) -> Image.Image:
    img = Image.open(path)
    if img.mode != "RGB":
        img = img.convert("RGB")
    return img


def fit_to_width(img: Image.Image, target_w: int) -> Image.Image:
    if img.width == target_w:
        return img
    scale = target_w / img.width
    new_h = max(1, int(round(img.height * scale)))
    return img.resize((target_w, new_h), Image.LANCZOS)


def _text_size(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.ImageFont) -> tuple[int, int]:
    l, t, r, b = draw.textbbox((0, 0), text, font=font)
    return r - l, b - t


def make_panel(src: Path, caption: str) -> Image.Image:
    """A single panel: source image scaled to PANEL_W with a caption bar on top."""
    img = fit_to_width(load_rgb(src), PANEL_W)
    panel = Image.new("RGB", (PANEL_W, CAPTION_H + img.height), BG)
    draw = ImageDraw.Draw(panel)
    draw.rectangle((0, 0, PANEL_W, CAPTION_H), fill=CAPTION_BG)
    tw, th = _text_size(draw, caption, CAPTION_FONT)
    draw.text(((PANEL_W - tw) // 2, (CAPTION_H - th) // 2 - 2), caption, fill=FG, font=CAPTION_FONT)
    panel.paste(img, (0, CAPTION_H))
    return panel


def grid_compose(panels: list[Image.Image], ncol: int, super_title: str) -> Image.Image:
    """Lay panels out in a grid (row-major), letterboxing rows to a common height."""
    nrow = (len(panels) + ncol - 1) // ncol
    # Per-row height = max panel height in that row (panels share PANEL_W already).
    rows = [panels[r * ncol:(r + 1) * ncol] for r in range(nrow)]
    row_h = [max(p.height for p in row) for row in rows]

    inner_w = ncol * PANEL_W + (ncol - 1) * GAP
    total_w = inner_w + 2 * PAD
    total_h = SUPER_H + PAD + sum(row_h) + GAP * (nrow - 1) + PAD

    canvas = Image.new("RGB", (total_w, total_h), BG)
    draw = ImageDraw.Draw(canvas)

    # Super-title, centred.
    tw, th = _text_size(draw, super_title, SUPER_FONT)
    draw.text(((total_w - tw) // 2, (SUPER_H - th) // 2), super_title, fill=FG, font=SUPER_FONT)

    y = SUPER_H + PAD
    for r, row in enumerate(rows):
        # Centre a short final row (e.g. 3 panels in a 2-col grid) — but for our
        # specs ncol divides the panel count, so this just left-aligns normally.
        used_w = len(row) * PANEL_W + (len(row) - 1) * GAP
        x0 = PAD + (inner_w - used_w) // 2
        x = x0
        for p in row:
            # Top-align each panel within the row band.
            canvas.paste(p, (x, y))
            x += PANEL_W + GAP
        y += row_h[r] + GAP
    return canvas


# ---------------------------------------------------------------------------
# Source path builders
# ---------------------------------------------------------------------------

def surface_path(dl: Path, cell: str, method: str, metric: str, mode: str) -> Path:
    return (dl / "Cartesian_Downsampling" / cell / "Roussos" / mode / method /
            "Continuous_Viz_Log" / f"Mean_{metric}_Total_Continuous_SemiLog.png")


def figts_path(dl: Path, cell: str, method: str, metric: str, mode: str, bin_tag: str) -> Path:
    return (dl / "Cartesian_Downsampling" / cell / "Roussos" / mode / method /
            "Plots" / "Figure_TS" / f"Figure_{bin_tag}_{metric}_{mode}.png")


# ---------------------------------------------------------------------------
# Figure specs
# ---------------------------------------------------------------------------

@dataclass
class SuppFig:
    name: str
    grid_by: str            # "method" | "cell_type"
    metric: str             # "Recall" | "Precision"
    kind: str               # "surface" | "figTS"
    fixed_cell: str = "Oligo"          # used when grid_by == "method"
    fixed_method: str = "all_methods"  # used when grid_by == "cell_type"
    bin_tag: str = ""       # "T" | "S" for kind == figTS
    requested: bool = False
    ncol: int = 2

    def panel_specs(self) -> list[tuple[str, str, str]]:
        """Return list of (cell, method, caption) for each panel."""
        if self.grid_by == "method":
            return [(self.fixed_cell, m, METHOD_LABEL[m]) for m in METHODS_4]
        return [(c, self.fixed_method, CELL_LABEL[c]) for c in CELLTYPES_3]

    def source_for(self, cell: str, method: str, mode: str) -> Path:
        if self.kind == "surface":
            return surface_path(DATALINK, cell, method, self.metric, mode)
        return figts_path(DATALINK, cell, method, self.metric, mode, self.bin_tag)

    def super_title(self, mode: str) -> str:
        axis = "donors x reads" if mode == "reads" else "donors x cells"
        if self.kind == "surface":
            what = f"{self.metric} across the {axis} grid (semi-log surface)"
        else:
            what = f"{self.metric} across mean-expression, {BIN_LABEL[self.bin_tag]} ({axis})"
        if self.grid_by == "method":
            return f"{CELL_LABEL[self.fixed_cell]}: {what}, by DE method"
        return f"{METHOD_LABEL[self.fixed_method]}: {what}, by cell type"


def _bymethod_figs() -> list[SuppFig]:
    """Per-method (4-panel) comparison figures for every cell type.

    For each cell type we build the requested surface + Tiny-bin recall/precision
    panels plus the Small-bin companions, so the downsampled grid is shown split
    by DE method (DESeq2 | edgeR | limma | all_methods) for Astro, Micro and Oligo.
    """
    figs: list[SuppFig] = []
    for cell in CELLTYPES_3:
        c = cell.lower()
        figs += [
            SuppFig(f"supp_{c}_recall_surface_bymethod", "method", "Recall", "surface", fixed_cell=cell, requested=True),
            SuppFig(f"supp_{c}_precision_surface_bymethod", "method", "Precision", "surface", fixed_cell=cell, requested=True),
            SuppFig(f"supp_{c}_recall_tiny_bymethod", "method", "Recall", "figTS", bin_tag="T", fixed_cell=cell, requested=True),
            SuppFig(f"supp_{c}_precision_tiny_bymethod", "method", "Precision", "figTS", bin_tag="T", fixed_cell=cell, requested=True),
            SuppFig(f"supp_{c}_recall_small_bymethod", "method", "Recall", "figTS", bin_tag="S", fixed_cell=cell),
            SuppFig(f"supp_{c}_precision_small_bymethod", "method", "Precision", "figTS", bin_tag="S", fixed_cell=cell),
        ]
    return figs


SUPP_FIGURES: list[SuppFig] = _bymethod_figs() + [
    # ---- cross-cell-type hierarchy at all_methods ----
    SuppFig("supp_allmethods_recall_surface_bycelltype", "cell_type", "Recall", "surface", ncol=3),
    SuppFig("supp_allmethods_precision_surface_bycelltype", "cell_type", "Precision", "surface", ncol=3),
]


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

DATALINK = DEFAULT_DATALINK  # rebound in main()


def build_one(fig: SuppFig, mode: str) -> tuple[Image.Image | None, list[str]]:
    """Return (image | None, missing_source_strings)."""
    panels: list[Image.Image] = []
    missing: list[str] = []
    for cell, method, caption in fig.panel_specs():
        src = fig.source_for(cell, method, mode)
        if not src.exists():
            missing.append(str(src))
            continue
        panels.append(make_panel(src, caption))
    if missing:
        return None, missing
    return grid_compose(panels, fig.ncol, fig.super_title(mode)), []


def main() -> None:
    global DATALINK
    mode = "cells"
    only_filter: list[str] = []
    argv = sys.argv[1:]
    i = 0
    while i < len(argv):
        a = argv[i]
        if a.startswith("--only="):
            only_filter = [s.strip() for s in a.split("=", 1)[1].split(",")]
        elif a == "--only" and i + 1 < len(argv):
            only_filter = [s.strip() for s in argv[i + 1].split(",")]; i += 1
        elif a.startswith("--mode="):
            mode = a.split("=", 1)[1].strip()
        elif a == "--mode" and i + 1 < len(argv):
            mode = argv[i + 1].strip(); i += 1
        elif a.startswith("--datalink="):
            DATALINK = Path(a.split("=", 1)[1]).expanduser().resolve()
        elif a == "--datalink" and i + 1 < len(argv):
            DATALINK = Path(argv[i + 1]).expanduser().resolve(); i += 1
        i += 1

    if mode not in ("cells", "reads"):
        print(f"ERROR: --mode must be 'cells' or 'reads', got '{mode}'", file=sys.stderr)
        sys.exit(1)

    figures = SUPP_FIGURES
    if only_filter:
        figures = [f for f in SUPP_FIGURES if any(sub in f.name for sub in only_filter)]
        if not figures:
            print(f"ERROR: --only {only_filter} matched no figures", file=sys.stderr)
            sys.exit(1)

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    composed_dir = DATALINK / "Composed_Figures" / "Supplementary"
    composed_dir.mkdir(parents=True, exist_ok=True)

    log_lines: list[str] = []
    built = 0
    skipped = 0
    for fig in figures:
        img, missing = build_one(fig, mode)
        if img is None:
            skipped += 1
            line = f"[SKIP] {fig.name} ({mode}): {len(missing)} missing source(s); first: {missing[0]}"
            print(line, file=sys.stderr)
            log_lines.append(line)
            continue
        out_path = OUT_ROOT / f"{fig.name}_{mode}.png"
        img.save(out_path, optimize=True)
        shutil.copy2(out_path, composed_dir / out_path.name)
        built += 1
        tag = "REQ" if fig.requested else "sug"
        size_kb = out_path.stat().st_size // 1024
        line = f"[{tag}] {fig.name} ({mode}): {img.width}x{img.height} -> {out_path} ({size_kb} KB)"
        print(line)
        log_lines.append(line)

    log_lines.append(f"built={built} skipped={skipped} mode={mode} datalink={DATALINK}")
    (OUT_ROOT / "build_log.txt").write_text("\n".join(log_lines) + "\n")
    print(f"\nDone: {built} built, {skipped} skipped ({mode} view).")


if __name__ == "__main__":
    main()

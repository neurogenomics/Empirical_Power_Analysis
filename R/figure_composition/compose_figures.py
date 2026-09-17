"""
Figure composition pipeline — non-destructive.

Purpose: Take existing PNGs from the paper's Datalink output and compose them into
Nature Comms-style faceted panels (A/B/C/D). Does NOT modify any source figure or
any paper draft.

Outputs: R/figure_composition/output/<variant>/<figure_name>.png
Log:     R/figure_composition/output/build_log.txt

Two figures are produced, each in three style variants so Chris can pick.

Figure 1 — Grid Oligo Composite:
    A (top-left, half): Oligo recall continuous surface (semi-log)
    B (top-right, half): Oligo precision continuous surface (semi-log)
    C (full-width, middle): Oligo recall discrete grid
    D (full-width, bottom): Oligo precision discrete grid

Figure 2 — Bulk Validation Composite (Oligo only, NoXIST correlation):
    A (full-width, top): Oligo bulk recovery stratified bar chart
    B (full-width, bottom): Oligo bulk-anchored LFC correlation (XIST excluded)

Variants:
    v1_tight     — minimal padding, inline bold A/B labels in panel corners
    v2_spaced    — medium padding, bold Nature-style labels top-left outside panels
    v3_bordered  — hairline panel borders, labels in a dedicated header strip

Run:
    python3 R/figure_composition/compose_figures.py
"""

from __future__ import annotations

import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from PIL import Image, ImageDraw, ImageFont

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent.parent  # .../Power_Analysis_Paper
DATALINK = PROJECT_ROOT / "Datalink"
OUT_ROOT = HERE / "output"

# Source images (see FIGURE_TRACKING.md for provenance and supplement mapping)
FIG1_SOURCES = {
    "A_recall_surface": DATALINK / "Cartesian_Downsampling/Oligo/Roussos/cells/all_methods/Continuous_Viz_Log/Mean_Recall_Total_Continuous_SemiLog.png",
    "B_precision_surface": DATALINK / "Cartesian_Downsampling/Oligo/Roussos/cells/all_methods/Continuous_Viz_Log/Mean_Precision_Total_Continuous_SemiLog.png",
    "C_recall_grid": DATALINK / "Cartesian_Downsampling/Oligo/Roussos/cells/all_methods/Plots/Mean_Recall_Total.png",
    "D_precision_grid": DATALINK / "Cartesian_Downsampling/Oligo/Roussos/cells/all_methods/Plots/Mean_Precision_Total.png",
}

# Figure 2 now focuses on oligodendrocytes only, using the NoXIST correlation.
FIG2_SOURCE_BAR = DATALINK / "Bulk_Analysis/Oligo/Recovery_Stats_Stratified.png"
FIG2_SOURCE_CORR = DATALINK / "Bulk_Analysis/Oligo/Compare_Correlation_BulkAnchored_NoXIST.png"

# Figure 3 (main text): cross-dataset reproducibility — TWO stacked panels.
# The raw IoU panel has been retired to supplementary at Chris's request, so
# the main composite is now: A = universe-controlled IoU, B = FDP heatmap.
# Both panels share the same x-axis labels (the universe-controlled IoU now
# renders rotated column labels via the heatmap script change in 02b).
FIG3_SOURCES = {
    "A_iou_controlled": DATALINK / "Full_Ground_Truth/Oligo/Diagnostic_Results/Expanded_IoU_Matrix_Universe_Controlled.png",
    "B_fdp": DATALINK / "Full_Ground_Truth/FDP_Heatmap.png",
}

# Figure 6 (FigPower) sources: A = power curve overlay (no legend), B = simple
# sensitivity curves (carries the shared legend at the bottom).
FIG6_SOURCES = {
    "A_overlay": DATALINK / "Power_Curve_Fitting/Plots/with_shallow/full_per_iteration/Power_Curve_Overlay_NoLegend.png",
    "B_sensitivity": DATALINK / "Power_Curve_Fitting/Plots/with_shallow/full_per_iteration/CellType_Sensitivity_Simple.png",
}

# ---------------------------------------------------------------------------
# Style variants
# ---------------------------------------------------------------------------

@dataclass
class Style:
    name: str
    outer_pad: int          # whitespace outside the whole figure
    gap: int                # gap between panels
    label_size: int         # panel label font size (pt-ish)
    label_offset: tuple[int, int]  # (x, y) relative to panel top-left
    label_outside: bool     # draw label in dedicated strip above panel
    label_strip: int        # height of header strip if label_outside
    panel_border: int       # border width around each panel (0 = none)
    bg: tuple[int, int, int] = (255, 255, 255)
    fg: tuple[int, int, int] = (0, 0, 0)


STYLES = [
    Style(
        name="v1_tight",
        outer_pad=12,
        gap=8,
        label_size=48,
        label_offset=(14, 8),
        label_outside=False,
        label_strip=0,
        panel_border=0,
    ),
    Style(
        name="v2_spaced",
        outer_pad=40,
        gap=36,
        label_size=60,
        label_offset=(-6, -12),
        label_outside=True,
        label_strip=80,
        panel_border=0,
    ),
    Style(
        name="v3_bordered",
        outer_pad=28,
        gap=22,
        label_size=54,
        label_offset=(8, 6),
        label_outside=False,
        label_strip=0,
        panel_border=0,
    ),
]


def load_font(size: int) -> ImageFont.ImageFont:
    # Try to find a bold sans-serif font; fall back to default.
    candidates = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
        "/usr/share/fonts/truetype/freefont/FreeSansBold.ttf",
    ]
    for c in candidates:
        if Path(c).exists():
            return ImageFont.truetype(c, size=size)
    return ImageFont.load_default()


# ---------------------------------------------------------------------------
# Composition helpers
# ---------------------------------------------------------------------------

def load_rgb(path: Path) -> Image.Image:
    img = Image.open(path)
    if img.mode != "RGB":
        img = img.convert("RGB")
    return img


def fit_to_width(img: Image.Image, target_w: int) -> Image.Image:
    """Resize preserving aspect ratio to a target width."""
    if img.width == target_w:
        return img
    scale = target_w / img.width
    new_h = max(1, int(round(img.height * scale)))
    return img.resize((target_w, new_h), Image.LANCZOS)


def fit_to_height(img: Image.Image, target_h: int) -> Image.Image:
    if img.height == target_h:
        return img
    scale = target_h / img.height
    new_w = max(1, int(round(img.width * scale)))
    return img.resize((new_w, target_h), Image.LANCZOS)


def hstack(imgs: list[Image.Image], gap: int, bg: tuple[int, int, int]) -> Image.Image:
    """Horizontally stack images scaled to the same height."""
    h = max(i.height for i in imgs)
    scaled = [fit_to_height(i, h) for i in imgs]
    total_w = sum(i.width for i in scaled) + gap * (len(scaled) - 1)
    canvas = Image.new("RGB", (total_w, h), bg)
    x = 0
    for im in scaled:
        canvas.paste(im, (x, 0))
        x += im.width + gap
    return canvas


def draw_panel(
    canvas: Image.Image,
    panel: Image.Image,
    box: tuple[int, int, int, int],  # (x, y, w, h)
    label: str,
    style: Style,
    font: ImageFont.ImageFont,
) -> None:
    """Paste panel into box (x, y, w, h) with letterbox and draw label."""
    x, y, w, h = box
    # Scale panel to fit inside box preserving aspect ratio (letterboxed).
    ar_panel = panel.width / panel.height
    ar_box = w / h
    if ar_panel >= ar_box:
        # Fit by width
        new_w = w
        new_h = max(1, int(round(w / ar_panel)))
    else:
        new_h = h
        new_w = max(1, int(round(h * ar_panel)))
    px = x + (w - new_w) // 2
    py = y + (h - new_h) // 2
    resized = panel.resize((new_w, new_h), Image.LANCZOS)
    canvas.paste(resized, (px, py))

    draw = ImageDraw.Draw(canvas)
    if style.panel_border > 0:
        draw.rectangle(
            (x - style.panel_border // 2, y - style.panel_border // 2,
             x + w + style.panel_border // 2, y + h + style.panel_border // 2),
            outline=style.fg,
            width=style.panel_border,
        )
    # Label
    lx = x + style.label_offset[0]
    ly = y + style.label_offset[1]
    if style.label_outside:
        # label sits just above the panel top-left
        ly = y - style.label_strip + 6
    draw.text((lx, ly), label, fill=style.fg, font=font)


# ---------------------------------------------------------------------------
# Figure 1: Grid Oligo composite
# ---------------------------------------------------------------------------

def build_fig1(style: Style) -> Image.Image:
    # Layout sketch:
    #   Row 1: A (half) | B (half)       — continuous surfaces, smaller rows
    #   Row 2: C (full)                  — recall discrete grid
    #   Row 3: D (full)                  — precision discrete grid
    font = load_font(style.label_size)

    target_total_w = 3200  # working canvas width in px
    pad = style.outer_pad
    gap = style.gap
    strip = style.label_strip if style.label_outside else 0

    inner_w = target_total_w - 2 * pad
    # Top row panels: each (inner_w - gap) / 2 wide
    top_panel_w = (inner_w - gap) // 2

    # Continuous surfaces are roughly square-ish; pick a modest height relative to width.
    # Use the real aspect ratio of the source to decide top row height.
    img_A = load_rgb(FIG1_SOURCES["A_recall_surface"])
    img_B = load_rgb(FIG1_SOURCES["B_precision_surface"])
    # Both should share aspect ratio; use A's.
    ar_top = img_A.width / img_A.height
    top_panel_h = max(1, int(round(top_panel_w / ar_top)))

    # Discrete grid panels — full width. Use the source aspect ratio.
    img_C = load_rgb(FIG1_SOURCES["C_recall_grid"])
    img_D = load_rgb(FIG1_SOURCES["D_precision_grid"])
    ar_grid = img_C.width / img_C.height
    grid_h = max(1, int(round(inner_w / ar_grid)))

    total_h = (
        pad
        + (strip + top_panel_h)
        + gap
        + (strip + grid_h)
        + gap
        + (strip + grid_h)
        + pad
    )
    canvas = Image.new("RGB", (target_total_w, total_h), style.bg)

    y = pad + strip
    # Row 1
    x = pad
    draw_panel(canvas, img_A, (x, y, top_panel_w, top_panel_h), "A", style, font)
    x += top_panel_w + gap
    draw_panel(canvas, img_B, (x, y, top_panel_w, top_panel_h), "B", style, font)

    # Row 2
    y += top_panel_h + gap + strip
    draw_panel(canvas, img_C, (pad, y, inner_w, grid_h), "C", style, font)

    # Row 3
    y += grid_h + gap + strip
    draw_panel(canvas, img_D, (pad, y, inner_w, grid_h), "D", style, font)

    return canvas


# ---------------------------------------------------------------------------
# Figure 1-Recall / Figure 1-Precision: split variants (Nature Comms style)
# ---------------------------------------------------------------------------

def _build_cartesian_pair(style: Style, surface_key: str, grid_key: str) -> Image.Image:
    """Two-panel composition: A = continuous surface (full width), B = discrete grid (full width).

    Used for Nature Comms-style split figures where each metric (recall/precision)
    lives in its own figure with two panels.
    """
    font = load_font(style.label_size)

    target_total_w = 3200
    pad = style.outer_pad
    gap = style.gap
    strip = style.label_strip if style.label_outside else 0
    inner_w = target_total_w - 2 * pad

    # Panel A: continuous surface — render at full inner width so the plot-area
    # width matches the discrete grid below. Upstream theme harmonisation
    # (legend key width/height and font size) in 04 and 07 keep legend strips
    # the same across the two panels.
    img_A = load_rgb(FIG1_SOURCES[surface_key])
    surface_target_w = inner_w
    ar_A = img_A.width / img_A.height
    surface_h = max(1, int(round(surface_target_w / ar_A)))

    # Panel B: discrete grid — full width so the cell values are readable.
    img_B = load_rgb(FIG1_SOURCES[grid_key])
    ar_B = img_B.width / img_B.height
    grid_h = max(1, int(round(inner_w / ar_B)))

    total_h = pad + (strip + surface_h) + gap + (strip + grid_h) + pad
    canvas = Image.new("RGB", (target_total_w, total_h), style.bg)

    # Row 1: surface, letterboxed inside the full inner width
    y = pad + strip
    draw_panel(canvas, img_A, (pad, y, inner_w, surface_h), "A", style, font)

    # Row 2: full-width grid
    y += surface_h + gap + strip
    draw_panel(canvas, img_B, (pad, y, inner_w, grid_h), "B", style, font)

    return canvas


def build_fig1_recall(style: Style) -> Image.Image:
    return _build_cartesian_pair(style, "A_recall_surface", "C_recall_grid")


def build_fig1_precision(style: Style) -> Image.Image:
    return _build_cartesian_pair(style, "B_precision_surface", "D_precision_grid")


# ---------------------------------------------------------------------------
# Figure 2: Bulk validation composite
# ---------------------------------------------------------------------------

def build_fig2(style: Style) -> Image.Image:
    """Oligo-only bulk validation: A = recovery bar chart, B = LFC correlation (NoXIST)."""
    font = load_font(style.label_size)

    target_total_w = 3200
    pad = style.outer_pad
    gap = style.gap
    strip = style.label_strip if style.label_outside else 0
    inner_w = target_total_w - 2 * pad

    # Panel A: oligo recovery stratified bar chart, full width, preserve source aspect.
    img_A = load_rgb(FIG2_SOURCE_BAR)
    ar_A = img_A.width / img_A.height
    bar_h = max(1, int(round(inner_w / ar_A)))

    # Panel B: oligo correlation plot (NoXIST), full width, preserve source aspect.
    img_B = load_rgb(FIG2_SOURCE_CORR)
    ar_B = img_B.width / img_B.height
    corr_h = max(1, int(round(inner_w / ar_B)))

    total_h = pad + (strip + bar_h) + gap + (strip + corr_h) + pad
    canvas = Image.new("RGB", (target_total_w, total_h), style.bg)

    y = pad + strip
    draw_panel(canvas, img_A, (pad, y, inner_w, bar_h), "A", style, font)

    y += bar_h + gap + strip
    draw_panel(canvas, img_B, (pad, y, inner_w, corr_h), "B", style, font)

    return canvas


# ---------------------------------------------------------------------------
# Figure 3: Cross-dataset reproducibility column (3 panels stacked)
# ---------------------------------------------------------------------------

def build_fig3(style: Style) -> Image.Image:
    """Single column: A = universe-controlled IoU, B = FDP heatmap.

    Each panel is rendered at a shared width so the matrices line up visually.
    The raw (without-universe-control) IoU panel has been moved to the
    supplementary at Chris's request — there it sits stacked with the
    universe-controlled version, both rendered with x-axis labels (see
    02b_Regenerate_Heatmaps.R) so the two matrices align column-for-column.
    """
    font = load_font(style.label_size)

    target_total_w = 1800  # compact canvas — fits within a page without filling it
    pad = style.outer_pad
    gap = style.gap
    strip = style.label_strip if style.label_outside else 0
    inner_w = target_total_w - 2 * pad

    img_A = load_rgb(FIG3_SOURCES["A_iou_controlled"])
    img_B = load_rgb(FIG3_SOURCES["B_fdp"])

    # Preserve each panel's own aspect ratio at the shared inner width.
    h_A = max(1, int(round(inner_w / (img_A.width / img_A.height))))
    h_B = max(1, int(round(inner_w / (img_B.width / img_B.height))))

    # Small header above each panel so A/B labels sit above the matrix top
    # rather than overlapping it. Only used for inside-label styles.
    header = 0 if style.label_outside else int(round(style.label_size * 0.9))

    total_h = (
        pad
        + (strip + header + h_A)
        + gap
        + (strip + header + h_B)
        + pad
    )
    canvas = Image.new("RGB", (target_total_w, total_h), style.bg)
    draw = ImageDraw.Draw(canvas)

    def place(label: str, img: Image.Image, y_top: int, h: int) -> None:
        if style.label_outside:
            draw_panel(canvas, img, (pad, y_top, inner_w, h), label, style, font)
        else:
            # Draw panel with empty label, then put label in the header above.
            draw_panel(canvas, img, (pad, y_top, inner_w, h), "", style, font)
            draw.text(
                (pad + style.label_offset[0], y_top - header + 2),
                label, fill=style.fg, font=font,
            )

    y = pad + strip + header
    place("A", img_A, y, h_A)
    y += h_A + gap + strip + header
    place("B", img_B, y, h_B)

    return canvas


# ---------------------------------------------------------------------------
# Figure 6 (FigPower): power curve overlay + practical sensitivity
# ---------------------------------------------------------------------------

def build_fig6(style: Style) -> Image.Image:
    """A/B composite for the power curve fitting figure.

    A (top): Power_Curve_Overlay (no legend) — predicted detection probability
             at median gene covariates across all five cell-type models.
    B (bottom): CellType_Sensitivity_Simple — same five models at a practical
             reference point. Carries the shared legend at the bottom, which
             also serves panel A.
    """
    font = load_font(style.label_size)

    target_total_w = 3200
    pad = style.outer_pad
    gap = style.gap
    strip = style.label_strip if style.label_outside else 0
    inner_w = target_total_w - 2 * pad

    img_A = load_rgb(FIG6_SOURCES["A_overlay"])
    img_B = load_rgb(FIG6_SOURCES["B_sensitivity"])

    h_A = max(1, int(round(inner_w / (img_A.width / img_A.height))))
    h_B = max(1, int(round(inner_w / (img_B.width / img_B.height))))

    header = 0 if style.label_outside else int(round(style.label_size * 0.9))

    total_h = (
        pad
        + (strip + header + h_A)
        + gap
        + (strip + header + h_B)
        + pad
    )
    canvas = Image.new("RGB", (target_total_w, total_h), style.bg)
    draw = ImageDraw.Draw(canvas)

    def place(label: str, img: Image.Image, y_top: int, h: int) -> None:
        if style.label_outside:
            draw_panel(canvas, img, (pad, y_top, inner_w, h), label, style, font)
        else:
            draw_panel(canvas, img, (pad, y_top, inner_w, h), "", style, font)
            draw.text(
                (pad + style.label_offset[0], y_top - header + 2),
                label, fill=style.fg, font=font,
            )

    y = pad + strip + header
    place("A", img_A, y, h_A)
    y += h_A + gap + strip + header
    place("B", img_B, y, h_B)

    return canvas


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

FIGURES: list[tuple[str, Callable[[Style], Image.Image]]] = [
    ("figure_cartesian_oligo", build_fig1),
    ("figure_cartesian_oligo_recall", build_fig1_recall),
    ("figure_cartesian_oligo_precision", build_fig1_precision),
    ("figure_bulk_validation", build_fig2),
    ("figure_reproducibility_column", build_fig3),
    ("figure_powercurve", build_fig6),
]


def verify_sources() -> None:
    missing = []
    for p in (
        list(FIG1_SOURCES.values())
        + [FIG2_SOURCE_BAR, FIG2_SOURCE_CORR]
        + list(FIG3_SOURCES.values())
        + list(FIG6_SOURCES.values())
    ):
        if not p.exists():
            missing.append(p)
    if missing:
        print("ERROR: missing source images:", file=sys.stderr)
        for p in missing:
            print(f"  {p}", file=sys.stderr)
        sys.exit(1)


def main() -> None:
    # --only fig1,fig3  builds only matching figures (substring match on name)
    only_filter: list[str] = []
    for arg in sys.argv[1:]:
        if arg.startswith("--only="):
            only_filter = [s.strip() for s in arg.split("=", 1)[1].split(",")]
        elif arg == "--only" and sys.argv.index(arg) + 1 < len(sys.argv):
            only_filter = [s.strip() for s in sys.argv[sys.argv.index(arg) + 1].split(",")]

    figures = FIGURES
    if only_filter:
        figures = [(n, b) for n, b in FIGURES if any(f in n for f in only_filter)]
        if not figures:
            print(f"ERROR: --only {only_filter} matched no figures", file=sys.stderr)
            sys.exit(1)
    else:
        verify_sources()

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    log_lines = []
    for style in STYLES:
        style_dir = OUT_ROOT / style.name
        style_dir.mkdir(parents=True, exist_ok=True)
        for fig_name, builder in figures:
            img = builder(style)
            out_path = style_dir / f"{fig_name}.png"
            img.save(out_path, optimize=True)
            size_kb = out_path.stat().st_size // 1024
            line = f"[{style.name}] {fig_name}: {img.width}x{img.height} -> {out_path} ({size_kb} KB)"
            print(line)
            log_lines.append(line)
    (OUT_ROOT / "build_log.txt").write_text("\n".join(log_lines) + "\n")

    # Copy the preferred style (v2_spaced) to the Datalink path the paper reads from.
    PUBLISH_STYLE = "v2_spaced"
    composed_dir = DATALINK / "Composed_Figures"
    composed_dir.mkdir(parents=True, exist_ok=True)
    publish_src = OUT_ROOT / PUBLISH_STYLE
    if publish_src.is_dir():
        for png in sorted(publish_src.glob("*.png")):
            dest = composed_dir / png.name
            shutil.copy2(png, dest)
            print(f"[publish] {png.name} -> {dest}")


if __name__ == "__main__":
    main()

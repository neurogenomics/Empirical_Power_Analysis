# Figure composition pipeline

Non-destructive image compositor that stitches existing paper PNGs into
Nature Comms-style faceted panels (A / B / C / D).

**Status:** trial — three style variants produced so Chris can pick one. Nothing
in the paper draft or in `Datalink/` has been modified.

## Run

```bash
cd /home/cramptc/cbert/data/workspace/packages/Power_Analysis_Paper
python3 R/figure_composition/compose_figures.py
```

Requires Python 3 + Pillow (already present on CBERT).

## Outputs

`R/figure_composition/output/<variant>/<figure>.png`

Variants:
- `v1_tight` — minimal padding, inline bold letters in each panel corner.
- `v2_spaced` — Nature-style: larger bold letters in a dedicated strip above each panel, generous gaps.
- `v3_bordered` — hairline panel borders + inline letters, for visual grouping.

Figures:
- `figure_cartesian_oligo.png` — oligo precision/recall continuous surfaces (A, B) + discrete recall grid (C) + discrete precision grid (D).
- `figure_bulk_validation.png` — stratified recovery bar charts for all three cell types (A, full width) + bulk-anchored LFC correlations for Astro/Micro/Oligo (B, C, D as thirds).

## Design notes

- No source PNGs are modified. The pipeline reads from `Datalink/` and writes only under `R/figure_composition/output/`.
- Working canvas width is 3200 px at lossless PNG. Aspect ratios of each panel are preserved via letterboxing — nothing is stretched.
- Panel A of Figure 2 (the "bulk bar chart") currently tiles the three per-cell-type stratified recovery bar charts (Astro | Micro | Oligo) horizontally into a single full-width strip. If you want a single combined bar chart instead, point me at which PNG (or generate a new one from `R/bulk_deg_processing/`) and I'll swap the source.
- To change a source image, edit `FIG1_SOURCES` / `FIG2_SOURCES_*` at the top of `compose_figures.py`. To change layout proportions, edit `build_fig1` / `build_fig2`.

## What to do once a variant is picked

See `FIGURE_TRACKING.md` for:
- The exact source PNGs feeding each panel.
- Which supplementary figure each source is currently used in.
- Suggested main-text reference edits.
- Suggested caption text for the new composites.
- A checklist for removing the absorbed panels from the supplementary file and renumbering.

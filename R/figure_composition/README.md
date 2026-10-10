# Figure composition pipeline

Composers that assemble the multi-panel main-text and supplementary figures
from outputs of the cartesian downsampling and power curve fitting pipelines.
No source PNG is modified: composers read from `Datalink/` (or rebuild panels
from checkpoints), stage under `R/figure_composition/output*/`, and publish the
paper copy to `Datalink/Composed_Figures/`.

## Use in the Reproduction Pipeline

The canonical end-to-end route is `pipeline/` (see `pipeline/README.md`).

Stage `compose` runs:

```bash
python3 R/figure_composition/compose_figures.py --only figure_cartesian_oligo
Rscript R/figure_composition/compose_recall_panels.R --cell Oligo --mode cells
# for mode in cells, reads:
Rscript R/figure_composition/compose_crossscore.R --mode <mode> --figset mega --style discrete
python3 R/figure_composition/compose_supplementary.py --mode <mode>
```

Stage `xgb_post` runs `Rscript R/figure_composition/compose_xgb_panels.R`.

Stage `ad_cartesian` runs, for each of Astro, Micro and Oligo, `compose_recall_panels.R` and
`compose_precision_panels.R` with `--mode cells --dataset Roussos_ADvsCtrl --composed-dir AD_vs_Control/Composed_Figures`.

## Scripts and the figures they produce

| Script | Output (paper copy) | Paper figure |
|--------|---------------------|--------------|
| `compose_figures.py --only figure_cartesian_oligo` | `Composed_Figures/figure_cartesian_oligo_recall.png`, `figure_cartesian_oligo_precision.png` | Fig. 2 (`FigRecall`), Fig. 4 (`FigPrecision`): Oligo continuous surface (A) + discrete grid (B) |
| `compose_recall_panels.R` | `Composed_Figures/figure_recall_panels_oligo_cells.png` | Fig. 3 (`FigRecallPanels`): recall by padj quartile (A) and by mean expression for tiny (B) and small (C) effect bins |
| `compose_xgb_panels.R` | `Power_Curve_Fitting/Models/xgboost/figure_xgb_composed.png` | Fig. 5 (`FigXGB`), panels A-E |
| `compose_crossscore.R` | `Composed_Figures/Cross_Score/crossscore_{recall,precision}_mega_{Astro,Micro,Oligo}_{cells,reads}_discrete_raw.png` | Supplementary cross-method grids |
| `compose_supplementary.py` | `Composed_Figures/Supplementary/supp_*_bymethod_{cells,reads}.png` (36 files) | Supplementary per-method surfaces and tiny/small panels |
| `interp_helpers.R` | (helper) | Local copy of the interpolation helpers from `R/cartesian_downsampling/07_Visualise_Continuous_Grid.R`, sourced by `compose_crossscore.R` |
| `compose_precision_panels.R` | `AD_vs_Control/Composed_Figures/figure_precision_panels_<cell>_cells.png` | AD vs control figures only (precision counterpart of the recall panels; placement in the manuscript not yet decided) |

Fig. 1 (`FigOverview`) uses the static `Figure1_Diagram.png` at the repository root.

Fig. 6 (`FigGeneGroups`) is produced by `R/gene_classification/`, not here.

`compose_figures.py` also defines `figure_cartesian_oligo` (4-panel combined),
`figure_bulk_validation`, `figure_reproducibility_column` and
`figure_powercurve`. These are obsolete (not used for the paper): the bulk
validation and cross-dataset reproducibility figures were removed with the
GTEx and cross-dataset analyses. The `--only figure_cartesian_oligo` filter
(substring match) builds only the oligo cartesian figures.

## Style variants (`compose_figures.py`)

Each figure is rendered in three styles under
`R/figure_composition/output/<variant>/`; `v2_spaced` is copied to
`Datalink/Composed_Figures/`.

- `v1_tight`: minimal padding, inline bold letters in each panel corner.
- `v2_spaced`: Nature-style, larger bold letters in a strip above each panel (published).
- `v3_bordered`: hairline panel borders + inline letters.

## Design notes

- Working canvas width is 3200 px at lossless PNG; panel aspect ratios are preserved by letterboxing.
- `compose_recall_panels.R` rebuilds its panels from the cartesian checkpoints rather than tiling PNGs, so facet annotations and subtitles can be dropped.
- To change a source image for `compose_figures.py`, edit the `FIG1_SOURCES` dictionary at the top of the script; layout lives in `build_fig1*`.
- Python dependencies (Pillow etc.) are pinned in `docker/python-requirements.txt`.

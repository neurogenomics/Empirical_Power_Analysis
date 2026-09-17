# Figure tracking — new composites + supplementary-figure implications

Everything needed to (a) move a chosen variant into the main draft, (b) strip
the absorbed panels out of the supplementary, and (c) re-point any main-text
references or captions. Nothing below has been edited in the draft yet — this
is a handover document, not a change log.

## Layouts currently produced (rerun with `python3 compose_figures.py`)

| Output file (per style dir) | Panels | Use case |
|---|---|---|
| `figure_cartesian_oligo.png` | A recall surface, B precision surface, C recall grid, D precision grid | Original 4-panel combined layout |
| `figure_cartesian_oligo_recall.png` | A recall surface, B recall grid | Split variant — recall-only (Nature Comms style) |
| `figure_cartesian_oligo_precision.png` | A precision surface, B precision grid | Split variant — precision-only (Nature Comms style) |
| `figure_bulk_validation.png` | A oligo recovery bar chart, B oligo LFC correlation (NoXIST) | Oligo-only bulk validation |

**Recommendation on split vs combined:** the 4-panel is dense and the discrete
grids have small cell labels. The two-figure split gives each grid more vertical
space and makes the cell values readable. Pick one path before moving to the
draft.

## Pending upstream change (staged, not yet propagated)

`R/bulk_deg_processing/05_Visualise_Bulk_vs_SC.R` has been edited in
`generate_corr_plot()` to:

- Shorten the cell-type label in the plot title: "Oligodendrocytes" → "Oligo",
  "Astrocytes" → "Astro", "Microglia" → "Micro".
- Override theme font sizes so title (12pt), legend title (9pt), legend text
  (8pt), and axis titles (11pt) / axis text (10pt) no longer crowd the plot
  area. Legend key size and spacing also reduced.

**Status:** edit is saved locally, but the ICL server's Syncthing peer
(`server-03`) was disconnected at the time of edit, so the file has not
transferred. To apply:

1. Wait for / restart Syncthing on server-03.
2. Confirm the file has arrived: `ssh icl-neuro stat
   /mnt/data/chrisc/Power_Analysis_Paper/R/bulk_deg_processing/05_Visualise_Bulk_vs_SC.R`
   — the modify time should reflect today.
3. Rerun the viz step: `ssh icl-neuro Rscript
   /mnt/data/chrisc/Power_Analysis_Paper/R/bulk_deg_processing/05_Visualise_Bulk_vs_SC.R`
   (or the appropriate pipeline entry point if the full pipeline is needed).
4. After ~60s of Syncthing backfill, rerun the composition:
   `python3 R/figure_composition/compose_figures.py`.

Until this runs, `figure_bulk_validation.png` still uses the old correlation
plot with the cramped "Oligodendrocytes" title and large legend.

---

## Figure 1 — Cartesian Oligo Composite (new main figure)

Replaces the current `FigRecall` chunk in `docs/paper_final_draft.rmd`
(currently a single-panel figure: the oligo recall continuous surface only).

### Panel sources

| Panel | Source PNG (under `Datalink/`) | Currently lives in | Proposed status |
|---|---|---|---|
| **A** (top-left, recall surface) | `Cartesian_Downsampling/Oligo/Roussos/cells/Continuous_Viz_Log/Mean_Recall_Total_Continuous_SemiLog.png` | **Main Figure `FigRecall`** (the current single-panel figure). | Keep in main, now panel A. |
| **B** (top-right, precision surface) | `Cartesian_Downsampling/Oligo/Roussos/cells/Continuous_Viz_Log/Mean_Precision_Total_Continuous_SemiLog.png` | **Not currently in paper.** (Precision continuous surface for oligo is not shown; only the three-cell discrete precision grid is in Supp 4.) | Promote to main fig. |
| **C** (mid, full width, recall grid) | `Cartesian_Downsampling/Oligo/Roussos/cells/Plots/Mean_Recall_Total.png` | **Supp `FigRecall-Supp2-DiscreteHeatmaps`** (chunk at line 545; three cell types side by side — Astro | Micro | **Oligo**). | Remove the oligo column from the supp (or remove the whole chunk if astro + micro should also move — see note below). |
| **D** (bottom, full width, precision grid) | `Cartesian_Downsampling/Oligo/Roussos/cells/Plots/Mean_Precision_Total.png` | **Supp `Figure2-Supp2-Precision-Cells`** (chunk at line 620; Astro | Micro | **Oligo**). | Remove the oligo column. |

### Main text references to update

All `\@ref(fig:FigRecall)` references continue to work — they'll now point at
the composite, which includes the recall surface they originally meant. No
rename needed. You may still want to tighten the prose since the figure now
carries more than just recall:

| Line | Current text (abridged) | Suggested edit |
|---|---|---|
| L85 | "The continuous recall surface for oligodendrocytes (Figure \\@ref(fig:FigRecall)) reveals the power landscape…" | "Figure \\@ref(fig:FigRecall) shows the oligodendrocyte power landscape: continuous recall and precision surfaces (panels A–B) and the underlying discrete grids (panels C–D)." |
| L85–86 | "Large-effect genes (|LFC| > 3) achieved 80% recall with as few as 10 donors…" | Change trailing "(Figure \\@ref(fig:FigRecall))" to "(Figure \\@ref(fig:FigRecall)A)". |
| L22 | "Precision exhibited more complex patterns … (Figure \\@ref(fig:FigRecall)—figure supplement 4)." | Point the within-oligo precision discussion at panels B and D of the main figure; keep the supp reference only for astro + micro precision grids. |
| L22 | "For large-effect oligodendrocyte genes (|LFC| > 3, n = 25), precision exceeded 80% across nearly all donor–cell configurations…" | Add "(Figure \\@ref(fig:FigRecall)B,D)". |

Astro and micro continuous surfaces in Supp 1 still exist and should stay
there. The astro + micro discrete grids (Supps 2 and 4) can remain, just with
the oligo column dropped, or you can keep the three-cell-type panels intact
in the supp for completeness and accept the small duplication — your call.

### Suggested caption

> **Figure X. Empirical power landscape for oligodendrocyte sex-biased DE.**
> (A) Continuous recall surface (semi-log scale) across donor count and
> cells per donor, stratified by |LFC| bin. (B) Matched continuous precision
> surface. Surfaces are bilinearly interpolated from the measured grid,
> averaged across six independent robustness runs. (C) Underlying discrete
> recall grid: mean recall across robustness runs at each tested
> donor × cells configuration. Red borders mark configurations achieving
> ≥ 80% recall. (D) Matching discrete precision grid. Panels A, C share the
> recall colour scale; B, D share the precision colour scale. Large-effect
> genes (|LFC| > 3) reach 80% recall at ≥ 10 donors with ≥ 1,500 cells;
> small-effect genes (|LFC| 0.2–2) remain below 80% recall across the entire
> grid. Oligodendrocytes are shown as the primary reference owing to their
> most complete ground truth (96% gene testability). Astrocyte and
> microglial surfaces are provided in Figure \\@ref(fig:FigRecall)—figure
> supplement 1; full three-cell-type discrete grids in figure
> supplements 2 and 4.

---

## Figure 2 — Bulk Validation Composite (Oligo-only, XIST excluded)

This composite does **not yet have a main-figure counterpart**. Everything
shown is presently in the supplementary. Promoting it would be a new
main-figure insertion; place it wherever the bulk validation paragraph sits
(currently in the ground-truth section around L79 of the prose and in Supps
9–11).

**Scope updated 2026-04-22:** only the oligodendrocyte panels are now in the
composite (the astrocyte and microglial bar charts + correlations stay in
the supplement). The correlation plot uses the **NoXIST** variant because
XIST is a huge outlier and visually dominates the scatter; Pearson R without
XIST is the number quoted in the caption anyway.

### Panel sources

| Panel | Source PNG | Currently lives in |
|---|---|---|
| **A** (full width, top) | `Bulk_Analysis/Oligo/Recovery_Stats_Stratified.png` | **Supp `Figure1-Supp10-Recovery-Stats`** (chunk at line 516; one of three per-cell-type bar charts). |
| **B** (full width, bottom) | `Bulk_Analysis/Oligo/Compare_Correlation_BulkAnchored_NoXIST.png` | **Supp `Figure1-Supp9-Bulk-SC-Correlation`** (line 504; the with-XIST version is what's there, the NoXIST version lives alongside it in `Datalink/` but isn't currently rendered in the supp). |

### Supplementary cleanup if this becomes a main figure

Neither source supplement is now fully redundant — the astro and micro
panels still live there. Options:

- **Minimal change:** leave Supp 9 and Supp 10 alone; they retain all three
  cell types and this new main figure simply elevates the oligo panels.
  Simplest and least disruptive. Recommended.
- **Or swap the Supp 9 oligo plot** from `Compare_Correlation_BulkAnchored.png`
  (with-XIST) to `Compare_Correlation_BulkAnchored_NoXIST.png` for consistency
  with the main figure. Your call.

No renumbering of supplements is required.

### Main text references to update

| Line | Current text | Suggested edit |
|---|---|---|
| L3 (big paragraph) | "…Figure \\@ref(fig:FigOverview)—figure supplement 10" (LFC correlations) | Add a reference to the new main figure's panel B for the oligo correlation; keep the supp pointer for astro + micro. |
| L3 | "…Figure \\@ref(fig:FigOverview)—figure supplement 11" (recovery stats) | Add a reference to panel A for oligo; keep the supp pointer for astro + micro. |

### Suggested caption

> **Figure Y. Orthogonal validation against GTEx bulk RNA-seq
> (oligodendrocytes).** (A) Recovery of the 79 tissue-consistent GTEx bulk
> sex-biased DEGs by Roussos oligodendrocyte single-cell pseudobulk,
> stratified by chromosome type (autosomal vs sex-chromosome) and bulk |LFC|
> bin. Bars show each gene's classification: Recovered (FDR < 0.05, matching
> LFC direction), Insignificant (tested but not significant or discordant
> direction), NA-Filtered (removed by DESeq2 independent filtering), or
> Absent (not in the single-cell gene set). (B) LFC concordance between
> GTEx bulk tissue and Roussos oligodendrocyte single-cell pseudobulk for
> genes significant in both modalities, **with XIST excluded** to remove the
> sex-chromosome outlier. Pearson R = 0.88. Autosomal genes in grey,
> X-chromosome genes in pink; translucent points indicate genes not
> significant in the single-cell analysis. Dashed line is the identity
> (y = x). Astrocyte and microglial bar charts and correlations are provided
> in Figure \\@ref(fig:FigOverview)—figure supplements 10 and 9.

---

## Supplementary renumbering checklist (after migration)

Current supplementary structure (abbreviated from `paper_final_draft.rmd`):

```
Supp 1   FDP (not affected)
...
Supp 8   Bulk volcano plot              <- keep
Supp 9   Bulk-SC LFC correlations        <- DELETE (moves to main Fig Y, B–D)
Supp 10  Recovery stats (stratified)     <- DELETE (moves to main Fig Y, A)
Supp 11+ ...                             <- renumber down by 2

FigRecall Supp 1  Astro+Micro surfaces   <- keep
FigRecall Supp 2  Discrete recall grids  <- drop oligo column, keep Astro+Micro
FigRecall Supp 3  Chromosome strat.       <- keep
FigRecall Supp 4  Precision grids        <- drop oligo column from cells-based only
```

You asked me not to reorder or renumber proactively in the draft itself — I
haven't touched it. The list above is the set of edits to apply once you pick
a variant and say the word.

---

## File inventory (what this pipeline wrote)

```
R/figure_composition/
├── README.md
├── FIGURE_TRACKING.md          <-- this file
├── compose_figures.py
└── output/
    ├── build_log.txt
    ├── v1_tight/
    │   ├── figure_cartesian_oligo.png
    │   └── figure_bulk_validation.png
    ├── v2_spaced/
    │   ├── figure_cartesian_oligo.png
    │   └── figure_bulk_validation.png
    └── v3_bordered/
        ├── figure_cartesian_oligo.png
        └── figure_bulk_validation.png
```

Nothing outside `R/figure_composition/` has been added, moved, or edited.

# Model Card: `predict_power.R` (XGBoost power model)

This card documents the gradient-boosted power model shipped in
`Datalink/Power_Curve_Fitting/Exports/xgboost/realistic/`. It follows the
[Model Cards for Model Reporting](https://arxiv.org/abs/1810.03993) framework.
The earlier logistic-regression model (pipeline v4.0.0) is retired and is not
the model reported in the paper.

---

## Model details

- **Model type**: XGBoost booster (gradient-boosted trees, `tree_method = "hist"`),
  objective `binary:logistic` on the per-design-point detection rate, rows
  weighted by `n_trials` (the number of iterations behind each rate).
- **Features** (12 design columns): six continuous features used on their raw
  scale (`abs_gt_lfc`, `gt_dispersion`, `gt_basemean`, `n_donors`, `n_cells`,
  `reads_per_cell`) plus one-hot `de_method` (deseq2, edger, limma) and
  `cell_type` (Astro, Micro, Oligo).
- **Monotone constraints**: predicted power is non-decreasing in donors, cells,
  reads per cell, |LFC| and baseMean, and non-increasing in dispersion.
- **Selected configuration**: config 25 of 36 (`alpha = 0`, `max_depth = 6`,
  `eta = 0.03`, `min_child_weight = 10`, `subsample = 0.8`,
  `colsample_bytree = 0.8`, 1,383 boosting rounds; grouped-CV log loss 0.1001).
- **Training code**: `R/power_curve_fitting/02_Fit_XGBoost.R` (search),
  `03_Select_XGBoost.R` (selection), `04_Export_XGBoost.R` (export).
- **Exported with**: pipeline v4.2.0, xgboost 3.2.1.1.
- **No closed-form equation**: the booster itself is the power function.

## Intended use

- **Power planning** for pseudobulk differential expression in snRNA-seq of
  human dorsolateral prefrontal cortex glia (astrocytes, microglia,
  oligodendrocytes), with a two-group donor-level contrast.
- **Triage of designs and genes** into "likely detected", "borderline" and
  "unlikely detected". The model is most reliable at the extremes of the power
  range (see Evaluation), so a three-class reading is recommended over treating
  the output as an exact probability.

### Out-of-scope uses

- Neurons, other glial subtypes, or non-brain tissues (not in the training data).
- Brain regions or cohorts other than PsychAD DLPFC (not validated).
- Single-cell-level or mixed-model DE methods; only DESeq2, edgeR and
  limma-voom on pseudobulk counts are represented.
- Designs or genes outside the training ranges (see Limitations).
- Reporting a precise power value (for example "power = 0.63") in the mid
  range as if it were exact.

---

## Training data

- **Source**: PsychAD snRNA-seq (Fullard et al. 2025), DLPFC only.
- **Donors**: 1,494 before filtering; after requiring age >= 40 years and
  >= 30 cells per donor in the target cell type, 1,293 (Astro), 1,214 (Micro)
  and 1,303 (Oligo).
- **Downsampling grid**: 12 donor counts (5 to 600) x 14 cells-per-donor
  levels (5 to 3,000, plus all available cells), 20 iterations per point,
  10 robustness runs, each with an independent sex-stratified 50/50 donor
  split into a truth half and a sampling half. A parallel mode subsamples by
  UMI budget per donor instead of cell count; both modes are pooled using the
  measured `n_cells` and `reads_per_cell` of each iteration.
- **DE methods**: DESeq2, edgeR and limma-voom, all LFCs shrunk with ashr.
- **Truth set**: genes with padj < 0.05 and |LFC| > 0.2 in the truth half
  (sex contrast). A gene counts as detected in a subsample when padj < 0.05.
- **Training rows**: detection rate per gene, cell type, DE method, design
  point and run, restricted to `n_cells >= 100` (the "realistic" sub-grid,
  `01e_Filter_Realistic.R`). Rows below 100 cells per donor are excluded
  because they do not correspond to designs anyone would plan.

## Evaluation

- **Protocol**: 5-fold cross-validation grouped by gene (all rows of a gene
  sit in one fold), early stopping on log loss. The 36 configurations
  (12 tree settings x 3 loss re-weighting exponents, `alpha` in {0, 0.5, 1})
  are ranked by out-of-fold bin-balanced MAE (`model_selection.csv`).
- **Extremes**: 97.6% (power 0 to 0.1) and 98.0% (power 0.9 to 1) of
  out-of-fold predictions fall within 0.10 of the observed detection rate.
- **Mid band**: MAE is 0.22 to 0.24 for observed power between 0.2 and 0.8,
  across all configurations tried.
- **Data balance**: only 7.8% of grid rows have observed power in (0.1, 0.9);
  detection is close to all-or-nothing for most gene and design combinations.
- **Classification view**: `09_Trinary_Classification.R` reports per-class
  precision, recall and F1 for several low / borderline / high partitions.

## Limitations and caveats

1. **Mid-band unreliability.** Between roughly 0.2 and 0.8 the prediction can
   be off by about 0.2 on average. Treat such outputs as "borderline".
2. **Extrapolation.** Trees predict a flat constant beyond the observed range
   of each feature. `predict_power()` warns when any continuous input falls
   outside the training range; designs below 100 cells per donor or above
   600 donors are not covered.
3. **Single cohort.** Trained and evaluated on PsychAD DLPFC only. There is
   no cross-cohort validation; transfer to other cohorts, protocols or
   sequencing depths is untested.
4. **Sex as the contrast.** The truth set is sex-differential expression,
   used because it gives a known, cohort-wide biological signal. Disease
   contrasts may have a different effect-size and dispersion structure, so
   predicted power for them is less certain.
5. **Cell type default.** `cell_type` defaults to Oligo, the cell type with
   the least depth-censored truth set. The ICE and swing analyses
   (`04b_CellType_ICE.R`, `04c_Swing_vs_Error.R`) show that defaulting to Oligo
   tends to under-state power for the other cell types, which is the
   conservative direction for planning.
6. **Gene inputs must be estimated.** `abs_lfc`, `basemean` and `dispersion`
   are on the DESeq2 scale of a pseudobulk fit. Values from pilot data or
   other tools need to be on the same scale.
7. **Demographic scope.** PsychAD donors are US-based; generalisation to
   other populations is untested.

## How to load and use

```r
# needs only R and the xgboost package
source("predict_power.R")   # loads power_model_bundle.rds from the same folder
predict_power(n_donors = 100, abs_lfc = 1.5, basemean = 50, dispersion = 0.5)
power_curve(donors = c(10, 25, 50, 100, 200, 500),
            abs_lfc = 1.5, basemean = 50, dispersion = 0.5)
```

- Required: `n_donors`, `abs_lfc`, `basemean`, `dispersion`.
- Optional, with defaults: `n_cells = 250`, `reads_per_cell = 5742` (PsychAD
  oligodendrocyte mean), `cell_type = "Oligo"`, `de_method = "deseq2"`.
- All arguments are vectorised and recycled to a common length.
- The bundle stores the booster as raw UBJ bytes and restores it with
  `xgb.load.raw()`. To load the booster yourself, use `booster.ubj` (or the
  bundle), not `booster.qs`, which is tied to the xgboost version that wrote it.
- `example_power_curve.csv` gives a worked power-versus-donors curve to check
  an installation against.

When citing a prediction, state the training context: PsychAD DLPFC glia,
sex contrast, the DE method used, and that values in the mid band are
approximate.

## Provenance

This card is versioned with the pipeline in git. Any refit, change of
training data or change of specification requires an update here.

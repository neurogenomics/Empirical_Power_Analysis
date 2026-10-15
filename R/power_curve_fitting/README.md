# Power Curve Fitting Pipeline (v4.2.0, XGBoost)

## 1. Purpose

Turns the cartesian downsampling grid into a predictive power model. For every
truth gene (padj < 0.05 and |LFC| > 0.2 in the truth half of a run), each
subsample records whether the gene was detected (padj < 0.05). These outcomes
are aggregated to a detection rate per gene and design point, and an XGBoost
model is trained to predict that rate from design parameters (donors, cells
per donor, reads per cell), gene characteristics (|LFC|, baseMean,
dispersion), DE method and cell type.

There is one power model, trained on the rows with `n_cells >= 100` (designs
below 100 cells per donor are not ones anyone would plan). It is exported as a
portable bundle with a `predict_power()` wrapper. See `docs/model-card-predict_power.md` for its intended use and
limitations.

Scope: PsychAD (`Roussos` shards, DLPFC) only; cell types Astro, Micro, Oligo;
DE methods DESeq2, edgeR and limma-voom.

## 2. Scripts in run order

The canonical route is the `pipeline/` stage graph (`pipeline/stages.tsv`,
`pipeline/stages/*.sh`). The stages that call this directory, in order:

| Stage | Script | What it does |
|-------|--------|--------------|
| `pcf_aggregate` | `00_Generate_GT_DDS_Checkpoints.R` | Re-fits DESeq2 on the truth half of each run to obtain per-gene dispersion and baseMean |
| | `00b_Cache_Shard_Read_Counts.R` | Caches per-cell read counts for every donor shard |
| | `01_Build_Fitting_Dataset.R` | Per-iteration detection outcomes with measured `n_cells` and `reads_per_cell`, cells and reads modes merged; run once per DE method |
| | `01b_Aggregate_Fitting_Data.R` | Collapses iterations to `detection_rate` and `n_trials` per (gene, cell type, DE method, mode, design point, run) |
| | `01e_Filter_Design.R` | Keeps rows with `n_cells >= 100`: the training table `aggregated_fitting_data.qs` |
| `xgb_search` (array of 108) | `02_Fit_XGBoost.R` | One hyperparameter configuration per task (`XGB_CONFIG_INDEX`), gene-grouped 5-fold CV, then a full fit |
| `xgb_select` | `03_Select_XGBoost.R` | Ranks the 108 configurations by out-of-fold bin-balanced MAE (per-fold values kept), promotes the winner, refits the five models shown in SHAP figures (winner per alpha, median, worst) |
| `xgb_shap` | `04_SHAP_Figures.R` | SHAP dependence and violin figures for the five selected models |
| `xgb_post` | `04_Export_XGBoost.R` | Writes the portable model bundle and `predict_power.R` |
| | `07_LowPower_Diversity.R` | Mean-power heatmaps of the low-power rows (Figure 5E) |
| | `08_Compare_Reweighted.R` | Calibration across configurations, SHAP comparison, and the hyperparameter bracketing with fold-paired standard errors (Figure 5B) |
| | `09_Trinary_Classification.R` | The selected model's out-of-fold per-class precision, recall and F1 for low / borderline / high power partitions (Figure 5C) |
| | `_regenerate_xgb_panels.R` | Figure 5 panels (then composed by `R/figure_composition/compose_xgb_panels.R`) |

The XGBoost stages run with 8 threads (`XGB_NTHREAD` and `OMP_NUM_THREADS`).
The exported bundle is then used by the `gene_class` stage
(`R/gene_classification/`, Figure 6).

### Removed

The retired logistic-regression model (v4.0.0 to v4.1.0, including the
master runner `run_pipeline.R`) and the supporting analyses that were not in
the paper (`01d`, `04b`, `04c`, `10`) are no longer in the repository; they
remain in the git history. The unit tests (`tests/test_power_curve_fitting.R`)
and toy test (`toy/run_toy_test.R`) date from the logistic model and do not
cover the XGBoost steps.

## 3. Inputs

| Input | Path |
|-------|------|
| Cartesian downsampling DE results and sampling metadata | `Datalink/Cartesian_Downsampling/{CT}/Roussos/{cells,reads}/{deseq2,edger,limma}/Run_{1-10}/` |
| Donor shards | `Datalink/Shards/{CT}/Roussos/` |
| BioMart cache | `Datalink/Full_Ground_Truth/gene_audit_biomart_cache.csv` |

## 4. Outputs

All paths are under `Datalink/Power_Curve_Fitting/`.

| Output | Path |
|--------|------|
| Truth-half DESeq2 checkpoints | `GT_DDS_Checkpoints/{CT}/GT_DDS_Run_{N}.qs` |
| Read-count cache | `Shard_Cache/{CT}/donor_{ID}_reads.qs` |
| Per-iteration detection, per DE method | `Fitting_Data/{CT}/unified_fitting_data_per_iteration__{method}.qs` |
| Aggregated training table (full grid) | `Fitting_Data/aggregated_fitting_data_all.qs` |
| Training table (`n_cells >= 100`) | `Fitting_Data/aggregated_fitting_data.qs`, plus `filter_summary.csv`, `filter_dropped_by_group.csv` |
| Power model | `Models/xgboost/` (`booster.ubj`, `feature_spec.qs`, `model_selection.csv`, `model_selection_folds.csv`, `cv_metrics.csv`, `calibration_residuals.csv`, `shap_selection.csv`, `shap_models/`, SHAP and importance tables) |
| Search | `Models/xgboost/search/config_{k}/` (one per configuration, with out-of-fold predictions) |
| Downstream analyses | `Models/xgboost/{shap_figures,comparison,lowpower_diversity,trinary_classification,supp_panels}/` |
| Figure 5 | `Models/xgboost/figure_xgb_composed.png` |
| Portable model | `Exports/xgboost/{power_model_bundle.rds,predict_power.R,example_power_curve.csv,README.md}` |

`model_selection.csv` is Supplementary Data 4. Load
boosters from `booster.ubj` (or the bundle via `xgb.load.raw()`), not from
`booster.qs`.

## 5. Model specification

- Objective `binary:logistic` on `detection_rate`, weighted by `n_trials`
  (times a normalised inverse-frequency power-bin weight when `alpha > 0`).
- Continuous features on their raw scale: `n_donors`, `n_cells`,
  `reads_per_cell`, `abs_gt_lfc`, `gt_basemean`, `gt_dispersion`. One-hot
  `cell_type` and `de_method`. 12 columns in total. `gt_basemean` is the
  gene's median raw pseudobulk count across the ground-truth-half donors
  (`median_raw` in the DE result tables).
- Monotone constraints: increasing in the first five continuous features,
  decreasing in dispersion.
- Search grid: `max_depth` {3, 4, 6, 10} x `eta` {0.01, 0.03, 0.1} x
  `min_child_weight` {1, 10, 25} x `alpha` {0, 0.5, 1} = 108 configurations,
  `subsample = 0.8`, `colsample_bytree = 0.8`, up to 20,000 rounds with early
  stopping after 50 rounds without improvement in log loss. Folds are grouped
  by gene.
- Selection: lowest out-of-fold bin-balanced MAE; the winner is the first row
  of `model_selection.csv`.
- Export defaults: `cell_type = "Oligo"`, `de_method = "deseq2"`,
  `n_cells = 250`, `reads_per_cell = 5742`.

## 6. Key configuration (`config.R`)

| Constant | Value |
|----------|-------|
| `DONOR_STEPS` | 5, 10, 15, 25, 50, 75, 100, 150, 200, 300, 500, 600 |
| `CELL_STEPS` | 5, 10, 20, 50, 75, 100, 250, 500, 750, 1000, 1500, 2000, 3000, all cells |
| `READ_STEPS` | 25k to 15M reads per donor (12 levels), plus all reads |
| `N_ROBUSTNESS_RUNS` / `N_ITERATIONS` | 10 / 20 |
| `GT_PADJ_THRESHOLD` / `GT_LFC_THRESHOLD` | 0.05 / 0.2 |
| `DETECT_PADJ_THRESHOLD` | 0.05 |
| `MIN_GENE_COUNT` | 10 |
| `MEAN_READS_PER_CELL` | Astro 9456, Micro 4426, Oligo 5742 |
| `CARTESIAN_DE_METHOD` | env var; `deseq2`, `edger` or `limma` |

`config.R` also still defines logistic-era settings (`DATASET_CONFIGS` with
depth-equalised Astro variants, `GENE_SETS`, `LRT_ALPHA`, `LOOKUP_*`), which
the XGBoost steps do not use.

XGBoost scripts are controlled by environment variables: `XGB_FITTING_DATA`
(default `aggregated_fitting_data.qs`), `XGB_MODEL_SUBDIR` (default empty, i.e.
`Models/xgboost/`; set only for side experiments), `XGB_CONFIG_INDEX` (array
task, 1 to 108), `XGB_MAX_NROUNDS` (20,000), `XGB_EARLY_STOP` (50) and
`XGB_NTHREAD` (fixed at 8 in the pipeline, since thread count can change
floating-point summation order). `XGB_MIN_N_CELLS` sets the 01e threshold
(default 100).

## 7. How to run

Only through the stage graph (see `pipeline/README.md`), on a PBS cluster or
with `--local` on one host:

```bash
pipeline/run.sh $WORK                  # submits all stages, including those above
pipeline/run.sh $WORK --status
pipeline/run.sh $WORK --only xgb_post  # rerun one stage (its dependencies must be complete)
```

`pipeline/stages/<stage>.sh` shows the exact scripts and environment variables
each stage uses.

## 8. Not tested

- No automated tests cover the XGBoost steps; the existing tests target the
  legacy logistic model.
- No cross-cohort validation: trained and evaluated on PsychAD only (the AD
  vs control test uses the same cohort).
- Extrapolation beyond the training ranges (trees predict a constant there;
  `predict_power()` warns).
- Accuracy of the truth-half DESeq2 dispersion estimates is assumed.

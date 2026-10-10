# Power Curve Fitting Pipeline (v4.2.0, XGBoost)

## 1. Purpose

Turns the cartesian downsampling grid into a predictive power model. For every
truth gene (padj < 0.05 and |LFC| > 0.2 in the truth half of a run), each
subsample records whether the gene was detected (padj < 0.05). These outcomes
are aggregated to a detection rate per gene and design point, and an XGBoost
model is trained to predict that rate from design parameters (donors, cells
per donor, reads per cell), gene characteristics (|LFC|, baseMean,
dispersion), DE method and cell type.

The production model is the "realistic" model, trained on the `n_cells >= 100`
sub-grid. It is exported as a portable bundle with a `predict_power()`
wrapper. See `docs/model-card-predict_power.md` for its intended use and
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
| | `01e_Filter_Realistic.R` | Keeps rows with `n_cells >= 100` (the realistic sub-grid) |
| `xgb_base_search` (array of 12) | `02_Fit_XGBoost.R` | Baseline model on the full grid: one hyperparameter config per task, gene-grouped 5-fold CV |
| `xgb_base_select` | `03_Select_XGBoost.R` | Picks the baseline winner by CV log loss; calibration and SHAP. Used only as a comparator in Figure 5 |
| `xgb_real_search` (array of 36) | `02_Fit_XGBoost.R` | Realistic model: 12 tree configs x 3 re-weighting exponents (`alpha` in {0, 0.5, 1}) |
| `xgb_real_select` | `03_Select_XGBoost.R` | Ranks the 36 configs by out-of-fold bin-balanced MAE, promotes the winner, refits the five models shown in SHAP figures (winner per alpha, median, worst) |
| `xgb_real_shap` | `04_SHAP_Figures.R` | SHAP dependence and violin figures for the selected realistic models |
| `xgb_post` | `08b_Baseline_OOF.R` | Out-of-fold predictions for the baseline winner |
| | `04_Export_XGBoost.R` | Writes the portable model bundle and `predict_power.R` |
| | `07_LowPower_Diversity.R` | Mean-power heatmaps of the low-power rows (Figure 5E) |
| | `08_Compare_Reweighted.R` | Calibration and SHAP comparison across the re-weighted models |
| | `09_Trinary_Classification.R` | Out-of-fold class F1 for low / borderline / high power partitions (Figure 5C) |
| | `_regenerate_xgb_panels.R` | Figure 5 panels (then composed by `R/figure_composition/compose_xgb_panels.R`) |

The exported bundle is then used by the `gene_class` stage
(`R/gene_classification/`, Figure 6).

### Supporting analyses (not called by `pipeline/`)

| Script | What it does |
|--------|--------------|
| `01d_Detection_Scatter.R` | Detection of small-effect truth genes in baseMean x dispersion space at fixed designs |
| `04b_CellType_ICE.R` | ICE / swing analysis of the `cell_type` input; justifies defaulting it to Oligo |
| `04c_Swing_vs_Error.R` | Compares the cell-type swing with the model's out-of-fold error |
| `10_PerGene_Bimodality.R` | Tests whether detection is close to all-or-nothing per gene |

### Legacy (logistic era, not used for the paper)

`02_Fit_Models.R`, `03_Visualise.R`,
`04_Cross_Validate_Modes.R`, `05_Export_Equations.R`,
`06_Compare_All_Models.R` and the master runner `run_pipeline.R` belong to the
retired logistic-regression model (v4.0.0 to v4.1.0). `run_pipeline.R` still
calls these scripts and does not run the XGBoost steps. The unit tests
(`tests/test_power_curve_fitting.R`) and toy test (`toy/run_toy_test.R`) also
cover the logistic model only.

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
| Aggregated training table (realistic) | `Fitting_Data/aggregated_fitting_data_realistic.qs`, plus `filter_summary.csv`, `filter_dropped_by_group.csv` |
| Baseline model | `Models/xgboost/` (`booster.ubj`, `feature_spec.qs`, `model_selection.csv`, `calibration*.csv`, SHAP outputs, `search/config_{k}/`) |
| Realistic (production) model | `Models/xgboost/realistic/` (`booster.ubj`, `feature_spec.qs`, `model_selection.csv`, `cv_metrics.csv`, `shap_selection.csv`, `search/config_{k}/` with per-config out-of-fold predictions) |
| Downstream analyses | `Models/xgboost/realistic/{shap_figures,comparison,lowpower_diversity,trinary_classification,supp_panels}/` |
| Figure 5 | `Models/xgboost/figure_xgb_composed.png` |
| Portable model | `Exports/xgboost/realistic/{power_model_bundle.rds,predict_power.R,example_power_curve.csv,README.md}` |

`model_selection.csv` for the realistic model is Supplementary Data 4. Load
boosters from `booster.ubj` (or the bundle via `xgb.load.raw()`), not from
`booster.qs`.

## 5. Model specification

- Objective `binary:logistic` on `detection_rate`, weighted by `n_trials`
  (times a normalised inverse-frequency power-bin weight when `alpha > 0`).
- Continuous features on their raw scale: `n_donors`, `n_cells`,
  `reads_per_cell`, `abs_gt_lfc`, `gt_basemean`, `gt_dispersion`. One-hot
  `cell_type` and `de_method`. 12 columns in total.
- Monotone constraints: increasing in the first five continuous features,
  decreasing in dispersion.
- Search grid: `max_depth` {3, 4, 6} x `eta` {0.03, 0.1} x `min_child_weight`
  {1, 10}, `subsample = 0.8`, `colsample_bytree = 0.8`, up to 5,000 rounds with
  early stopping after 50. Folds are grouped by gene.
- Selected realistic model: config 25 (`alpha = 0`, `max_depth = 6`,
  `eta = 0.03`, `min_child_weight = 10`, 1,383 rounds).
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
(`aggregated_fitting_data_realistic.qs` for the production model),
`XGB_MODEL_SUBDIR` (`realistic`), `XGB_CONFIG_INDEX` (array task), and
`XGB_NTHREAD` (fixed at 8 in the pipeline, since thread count can change
floating-point summation order). `XGB_MIN_N_CELLS` sets the 01e threshold
(default 100).

## 7. How to run

End to end, on a PBS cluster, via the stage graph (see `pipeline/README.md`):

```bash
pipeline/run.sh $WORK            # submits all stages, including those above
pipeline/run.sh $WORK --status
```

To run one script by hand from the repository root, set the same environment
variables the stage scripts set, for example the production export:

```bash
XGB_MODEL_SUBDIR=realistic \
XGB_FITTING_DATA=aggregated_fitting_data_realistic.qs \
Rscript R/power_curve_fitting/04_Export_XGBoost.R
```

## 8. Not tested

- No automated tests cover the XGBoost steps; the existing tests target the
  legacy logistic model.
- No cross-cohort validation: trained and evaluated on PsychAD only.
- Extrapolation beyond the training ranges (trees predict a constant there;
  `predict_power()` warns).
- Accuracy of the truth-half DESeq2 dispersion estimates is assumed.

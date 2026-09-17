# Portable XGBoost Power Model

Exported 2026-08-11 from the `realistic` model (pipeline v4.2.0, xgboost 3.2.1.1).

## Files
- `power_model_bundle.rds` -- self-contained: serialised booster + feature spec + defaults + training ranges + diagnostics.
- `predict_power.R` -- source() this; needs only base R + the `xgboost` package.
- `example_power_curve.csv` -- a worked power-vs-donors curve to check against.

## Quick start
```r
source("predict_power.R")
predict_power(n_donors = 100, abs_lfc = 1.5, basemean = 50, dispersion = 0.5)
power_curve(donors = c(10, 25, 50, 100, 200, 500),
            abs_lfc = 1.5, basemean = 50, dispersion = 0.5)
```

## Inputs
Required: `n_donors`, `abs_lfc`, `basemean`, `dispersion`.
Optional (defaults in brackets): `n_cells` [250], `reads_per_cell` [5742], `cell_type` [Oligo], `de_method` [deseq2].
All arguments accept vectors and are recycled to a common length.

## Notes
- Output is predicted power (detection probability) in [0, 1]; there is no closed-form equation (gradient-boosted trees).
- `cell_type`, `de_method` and `reads_per_cell` are kept in the model but default to PsychAD/Oligo settings; override to vary them.
- Predictions outside the training grid warn: tree models extrapolate as a flat constant beyond the observed range.
- Reference cell type is Oligo (least depth-censored truth set). Winner: config 25, alpha=0, CV logloss=0.1001.

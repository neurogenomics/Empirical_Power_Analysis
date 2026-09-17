# N-ary classification metrics

Generated 2026-08-06 10:10:22 UTC by 09_Trinary_Classification.R.

## Partition schemes evaluated

  binary                breaks=[0.5]  classes=unpowered / powered
  trinary_02_08         breaks=[0.2, 0.8]  classes=unpowered / contested / powered
  trinary_04_08         breaks=[0.4, 0.8]  classes=definitely_not / maybe / yes
  trinary_01_09         breaks=[0.1, 0.9]  classes=unpowered / contested / powered
  four_class            breaks=[0.2, 0.5, 0.8]  classes=unpowered / probably_not / probably_yes / powered
  quartile              breaks=[0.25, 0.5, 0.75]  classes=q1_low / q2_low_mid / q3_high_mid / q4_high

## Files

- `summary_by_scheme.csv`         one row per (model, scheme): accuracy, macro-F1, weighted-F1
- `per_class_metrics.csv`         precision / recall / F1 per (model, scheme, class)
- `confusion_winner_<scheme>.csv` confusion matrix for the realistic winner (K x K)
- `best_scheme_per_model.csv`     best-scoring scheme per model (by weighted-F1)
- `class_f1_grid.png`             faceted overview: all schemes x winner vs baseline
- `class_f1_<scheme>.png`         individual bar chart per scheme

## Mode

`XGB_TRINARY_MODE = insample`

**insample**: booster applied to its own training data. Overstates true generalisation.

## Caveats

- The 'baseline' model was trained on the unfiltered 7.2M-row grid and is
  evaluated here on the 3.4M-row realistic subset. This is not a symmetric
  comparison: the baseline saw more training data of a different distribution.
  See xgboost_section_plan.md section 7 for a note on this asymmetry.
- 'Weighted-F1' weights each class by its support (row count). For schemes
  where unpowered rows dominate (>85% of the grid), even a degenerate model
  that predicts 'unpowered' everywhere scores a high weighted-F1. Macro-F1
  is the fairer metric.
- Boundaries are set on the detection_rate (fraction of 20 iterations), not
  on the XGBoost predicted probability. Both are discretised using the same
  break-points for comparability.


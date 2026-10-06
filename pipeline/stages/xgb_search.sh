# shellcheck shell=bash
# Stage xgb_search (array:108): hyperparameter search for the power model,
# trained on the n_cells >= 100 design rows. 4 depths x 3 learning rates x
# 3 min_child_weights x 3 loss re-weighting exponents alpha in {0, 0.5, 1};
# task k fits configuration k of 02_Fit_XGBoost.R's param_grid with
# gene-grouped 5-fold CV (early stopping, cap 20,000 rounds) and saves its
# out-of-fold predictions.
source "$(dirname "${BASH_SOURCE[0]}")/_xgb_env.sh"
stage_main() {
  epa_init_ws
  epa_in pcf_aggregate fitting_data
  epa_xgb_env
  epa_env XGB_CONFIG_INDEX "$EPA_TASK"
  epa_R R/power_curve_fitting/02_Fit_XGBoost.R
  epa_out "$EPA_TASK_OUT" "Power_Curve_Fitting/Models/xgboost/search/config_$EPA_TASK"
}

# shellcheck shell=bash
# Stage xgb_real_search (array:36): hyperparameter search for the production
# ("realistic", n_cells >= 100) model: 12 tree configurations x 3 loss
# re-weighting exponents alpha in {0, 0.5, 1}; task k fits config k.
source "$(dirname "${BASH_SOURCE[0]}")/_xgb_env.sh"
stage_main() {
  epa_init_ws
  epa_in pcf_aggregate fitting_data
  epa_xgb_env realistic
  epa_env XGB_CONFIG_INDEX "$EPA_TASK"
  epa_R R/power_curve_fitting/02_Fit_XGBoost.R
  epa_out "$EPA_TASK_OUT" "Power_Curve_Fitting/Models/xgboost/realistic/search/config_$EPA_TASK"
}

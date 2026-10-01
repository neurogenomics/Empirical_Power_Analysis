# shellcheck shell=bash
# Stage xgb_base_search (array:12): hyperparameter search for the baseline
# (full-grid) XGBoost model; task k fits config k of 02_Fit_XGBoost.R's grid
# with grouped 5-fold CV. The baseline model is only used as a comparator in
# Figure 5.
source "$(dirname "${BASH_SOURCE[0]}")/_xgb_env.sh"
stage_main() {
  epa_init_ws
  epa_in pcf_aggregate fitting_data
  epa_xgb_env baseline
  epa_env XGB_CONFIG_INDEX "$EPA_TASK"
  epa_R R/power_curve_fitting/02_Fit_XGBoost.R
  epa_out "$EPA_TASK_OUT" "Power_Curve_Fitting/Models/xgboost/search/config_$EPA_TASK"
}

# shellcheck shell=bash
# Stage xgb_real_select: rank the 36 realistic configurations by out-of-fold
# bin-balanced MAE (model_selection.csv = Supplementary Data 4), promote the
# winner, and refit the five models whose SHAP values the paper shows
# (winner per alpha, median, worst).
source "$(dirname "${BASH_SOURCE[0]}")/_xgb_env.sh"
stage_main() {
  epa_init_ws
  epa_in pcf_aggregate fitting_data
  epa_in_array xgb_real_search
  epa_xgb_env realistic
  epa_R R/power_curve_fitting/03_Select_XGBoost.R
  epa_require Power_Curve_Fitting/Models/xgboost/realistic/booster.ubj \
              Power_Curve_Fitting/Models/xgboost/realistic/model_selection.csv \
              Power_Curve_Fitting/Models/xgboost/realistic/cv_metrics.csv
  epa_out model_realistic Power_Curve_Fitting/Models/xgboost/realistic
  epa_complete
}

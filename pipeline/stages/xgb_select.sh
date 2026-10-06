# shellcheck shell=bash
# Stage xgb_select: rank every search configuration by out-of-fold
# bin-balanced MAE, overall and per CV fold (model_selection.csv = Supplementary
# Data 4, model_selection_folds.csv), promote the winner to Models/xgboost/, and
# refit the five models whose SHAP values the paper shows (winner per alpha,
# median, worst).
source "$(dirname "${BASH_SOURCE[0]}")/_xgb_env.sh"
stage_main() {
  epa_init_ws
  epa_in pcf_aggregate fitting_data
  epa_in_array xgb_search
  epa_xgb_env
  epa_R R/power_curve_fitting/03_Select_XGBoost.R
  local m=Power_Curve_Fitting/Models/xgboost
  epa_require "$m/booster.ubj" "$m/model_selection.csv" "$m/model_selection_folds.csv" \
              "$m/cv_metrics.csv" "$m/shap_selection.csv"
  # The per-configuration search outputs (out-of-fold predictions for 108
  # configurations) are only needed by xgb_post, so they ship separately.
  epa_out model "$m/*.*" "$m/shap_models"
  epa_out search "$m/search"
  epa_complete
}

# shellcheck shell=bash
# Stage xgb_base_select: pick the best baseline configuration, refit it and
# compute calibration and SHAP summaries (03_Select_XGBoost.R).
source "$(dirname "${BASH_SOURCE[0]}")/_xgb_env.sh"
stage_main() {
  epa_init_ws
  epa_in pcf_aggregate fitting_data
  epa_in_array xgb_base_search
  epa_xgb_env baseline
  epa_R R/power_curve_fitting/03_Select_XGBoost.R
  epa_require Power_Curve_Fitting/Models/xgboost/booster.ubj \
              Power_Curve_Fitting/Models/xgboost/feature_spec.qs
  epa_out model_baseline \
    "Power_Curve_Fitting/Models/xgboost/*.*" \
    Power_Curve_Fitting/Models/xgboost/search
  epa_complete
}

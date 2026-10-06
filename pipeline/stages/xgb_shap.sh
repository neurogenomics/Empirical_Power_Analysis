# shellcheck shell=bash
# Stage xgb_shap: SHAP dependence / violin figures for the five selected models
# (04_SHAP_Figures.R).
source "$(dirname "${BASH_SOURCE[0]}")/_xgb_env.sh"
stage_main() {
  epa_init_ws
  epa_in pcf_aggregate fitting_data
  epa_in xgb_select model
  epa_xgb_env
  epa_env XGB_SHAP_MAX_ROWS 200000
  epa_R R/power_curve_fitting/04_SHAP_Figures.R
  local d=Power_Curve_Fitting/Models/xgboost/shap_figures
  epa_require "$d/median/shap_dependence.png" "$d/median/shap_violin.png" \
              "$d/winner_alpha_0/shap_dependence.png" "$d/winner_alpha_0.5/shap_dependence.png" \
              "$d/winner_alpha_1/shap_dependence.png"
  epa_out shap_figures "$d"
  epa_complete
}

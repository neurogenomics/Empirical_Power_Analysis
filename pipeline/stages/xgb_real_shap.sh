# shellcheck shell=bash
# Stage xgb_real_shap: SHAP dependence / violin figures for the selected
# realistic models (04_SHAP_Figures.R).
source "$(dirname "${BASH_SOURCE[0]}")/_xgb_env.sh"
stage_main() {
  epa_init_ws
  epa_in pcf_aggregate fitting_data
  epa_in xgb_real_select model_realistic
  epa_xgb_env realistic
  epa_env XGB_SHAP_MAX_ROWS 200000
  epa_env EPA_XGB_TAG reproduction
  epa_R R/power_curve_fitting/04_SHAP_Figures.R
  local d=Power_Curve_Fitting/Models/xgboost/realistic/shap_figures
  epa_require "$d/median/shap_dependence.png" "$d/median/shap_violin.png" \
              "$d/winner_alpha_0/shap_dependence.png" "$d/winner_alpha_0.5/shap_dependence.png" \
              "$d/winner_alpha_1/shap_dependence.png"
  epa_out shap_figures "$d"
  epa_complete
}

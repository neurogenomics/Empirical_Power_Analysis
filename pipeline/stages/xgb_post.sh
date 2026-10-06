# shellcheck shell=bash
# Stage xgb_post: everything downstream of the selected XGBoost model, up to
# main-text Figure 5 and the shareable model bundle.
#   04_Export_XGBoost.R         power_model_bundle.rds + predict_power.R (used by Figure 6)
#   07_LowPower_Diversity.R     mean-power heatmap (Fig. 5 panel E)
#   08_Compare_Reweighted.R     calibration / SHAP comparison tables, hyperparameter bracketing
#   09_Trinary_Classification   per-class precision / recall / F1 (Fig. 5 panel C), out-of-fold
#   _regenerate_xgb_panels.R    calibration, precision/recall/F1, SHAP and bracketing panels
#   compose_xgb_panels.R        figure_xgb_composed.png
source "$(dirname "${BASH_SOURCE[0]}")/_xgb_env.sh"
stage_main() {
  epa_init_ws
  epa_in pcf_aggregate fitting_data
  epa_in xgb_select model search
  local x=Power_Curve_Fitting/Models/xgboost
  epa_xgb_env

  epa_step "export the model bundle"
  epa_R R/power_curve_fitting/04_Export_XGBoost.R
  epa_step "low-power diversity"
  epa_R R/power_curve_fitting/07_LowPower_Diversity.R
  epa_step "compare search configurations"
  epa_R R/power_curve_fitting/08_Compare_Reweighted.R
  epa_step "per-class precision / recall / F1, out-of-fold"
  epa_env XGB_TRINARY_MODE oof
  epa_env XGB_TRINARY_CONFIGS winner
  epa_R R/power_curve_fitting/09_Trinary_Classification.R
  epa_step "Figure 5 panels"
  epa_R R/power_curve_fitting/_regenerate_xgb_panels.R

  # compose_xgb_panels.R reads its five panels from supp_panels/.
  local sp="$WS/Datalink/$x/supp_panels" pv="$WS/Datalink/$x/comparison/_panel_preview"
  mkdir -p "$sp"
  cp "$pv/calibration_curves_2x3.png"     "$sp/xgb_supp_calibration.png"
  cp "$pv/class_prf_grid.png"             "$sp/xgb_supp_prf_grid.png"
  cp "$pv/shap_importance_pretty.png"     "$sp/xgb_supp_shap.png"
  cp "$pv/hyperparameter_bracketing.png"  "$sp/xgb_supp_hp_bracketing.png"
  cp "$WS/Datalink/$x/lowpower_diversity/density2d/abs_gt_lfc_vs_n_donors_meanpower.png" \
     "$sp/xgb_supp_meanpower_heatmap.png"
  epa_step "compose Figure 5"
  epa_R R/figure_composition/compose_xgb_panels.R

  epa_require "$x/figure_xgb_composed.png" \
              Power_Curve_Fitting/Exports/xgboost/power_model_bundle.rds
  epa_out xgb_post \
    "$x/figure_xgb_composed.png" \
    "$x/comparison" "$x/lowpower_diversity" "$x/trinary_classification" "$x/supp_panels" \
    Power_Curve_Fitting/Exports/xgboost
  epa_complete
}

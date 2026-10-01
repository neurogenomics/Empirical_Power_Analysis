# shellcheck shell=bash
# Stage xgb_post: everything downstream of the two selected XGBoost models, up to
# main-text Figure 5 and the shareable model bundle.
#   08b_Baseline_OOF.R          out-of-fold predictions of the baseline winner
#   04_Export_XGBoost.R         power_model_bundle.rds + predict_power.R (used by Figure 6)
#   07_LowPower_Diversity.R     mean-power heatmap (Fig. 5 panel E)
#   08_Compare_Reweighted.R     calibration / SHAP comparison tables
#   09_Trinary_Classification   class F1 (Fig. 5 panel C), out-of-fold
#   _regenerate_xgb_panels.R    calibration, F1, SHAP and MAE panels
#   compose_xgb_panels.R        figure_xgb_composed.png
source "$(dirname "${BASH_SOURCE[0]}")/_xgb_env.sh"
stage_main() {
  epa_init_ws
  epa_in pcf_aggregate fitting_data
  epa_in xgb_base_select model_baseline
  epa_in xgb_real_select model_realistic
  local x=Power_Curve_Fitting/Models/xgboost r=Power_Curve_Fitting/Models/xgboost/realistic

  epa_step "baseline out-of-fold predictions"
  epa_xgb_env baseline
  epa_R R/power_curve_fitting/08b_Baseline_OOF.R

  epa_xgb_env realistic
  epa_step "export the realistic model bundle"
  epa_R R/power_curve_fitting/04_Export_XGBoost.R
  epa_step "low-power diversity"
  epa_R R/power_curve_fitting/07_LowPower_Diversity.R
  epa_step "compare re-weighted models"
  epa_R R/power_curve_fitting/08_Compare_Reweighted.R
  epa_step "class F1, out-of-fold"
  epa_env XGB_TRINARY_MODE oof
  epa_env XGB_TRINARY_CONFIGS winner
  epa_R R/power_curve_fitting/09_Trinary_Classification.R
  epa_step "Figure 5 panels"
  epa_R --cwd R/power_curve_fitting _regenerate_xgb_panels.R

  # compose_xgb_panels.R reads its five panels from realistic/supp_panels/.
  local sp="$WS/Datalink/$r/supp_panels" pv="$WS/Datalink/$r/comparison/_panel_preview"
  mkdir -p "$sp"
  cp "$pv/calibration_curves_2x3.png"     "$sp/xgb_supp_calibration.png"
  cp "$pv/class_f1_grid_v2.png"           "$sp/xgb_supp_f1_grid.png"
  cp "$pv/shap_importance_pretty.png"     "$sp/xgb_supp_shap.png"
  cp "$pv/bin_balanced_MAE_pretty.png"    "$sp/xgb_supp_mae_by_alpha.png"
  cp "$WS/Datalink/$r/lowpower_diversity/density2d/abs_gt_lfc_vs_n_donors_meanpower.png" \
     "$sp/xgb_supp_meanpower_heatmap.png"
  epa_step "compose Figure 5"
  epa_R R/figure_composition/compose_xgb_panels.R

  epa_require "$x/figure_xgb_composed.png" \
              Power_Curve_Fitting/Exports/xgboost/realistic/power_model_bundle.rds
  epa_out xgb_post \
    "$x/oof_predictions.qs" "$x/figure_xgb_composed.png" \
    "$r/comparison" "$r/lowpower_diversity" "$r/trinary_classification" "$r/supp_panels" \
    Power_Curve_Fitting/Exports/xgboost/realistic
  epa_complete
}

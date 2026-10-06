# shellcheck shell=bash
# Stage ad_validate: the sex-trained power model applied to the AD vs control
# grid as an external test set (R/ad_contrast/03_Validate_Power_Model.R).
source "$(dirname "${BASH_SOURCE[0]}")/_xgb_env.sh"
stage_main() {
  epa_init_ws
  epa_in pcf_aggregate fitting_data
  epa_in ad_pcf ad_fitting_data
  epa_in xgb_select model
  epa_in xgb_post
  epa_xgb_env
  epa_R R/ad_contrast/03_Validate_Power_Model.R
  local d=AD_vs_Control/Power_Model_Test
  epa_require "$d/test_metrics.csv" "$d/per_class_metrics.csv" "$d/trinary_skewed_prf.png"
  epa_out ad_validate "$d"
  epa_complete
}

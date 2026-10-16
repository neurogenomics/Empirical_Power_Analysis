# shellcheck shell=bash
# Stage supp_data: Supplementary Data 1 and 4 (R/supplementary_data/build_supplementary_data.R).
stage_main() {
  epa_init_ws
  epa_in eda
  epa_in fgt
  epa_in xgb_select model
  epa_R R/supplementary_data/build_supplementary_data.R
  epa_out supp_data Supplementary_Data
  epa_complete
}

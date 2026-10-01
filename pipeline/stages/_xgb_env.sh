# shellcheck shell=bash
# Container environment for the XGBoost stages. XGB_NTHREAD is fixed at 8 in
# every XGBoost stage (as in the paper's runs) because the thread count can
# change floating-point summation order.
epa_xgb_env() {  # baseline | realistic
  epa_env XGB_NTHREAD 8
  if [ "$1" = realistic ]; then
    epa_env XGB_MODEL_SUBDIR realistic
    epa_env XGB_FITTING_DATA aggregated_fitting_data_realistic.qs
  fi
}

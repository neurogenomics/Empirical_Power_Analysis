# shellcheck shell=bash
# Container environment for the XGBoost stages. XGB_NTHREAD is fixed at 8 in
# every XGBoost stage because the thread count can change floating-point
# summation order. XGBoost caps its threads at the OpenMP maximum, which
# epa_init_ws pins to 1 for every other stage, so OMP_NUM_THREADS must be
# raised to match or XGBoost silently runs single-threaded. The scripts'
# defaults already point at the single power model (Models/xgboost/, trained on
# Fitting_Data/aggregated_fitting_data.qs).
epa_xgb_env() {
  epa_env XGB_NTHREAD 8
  epa_env OMP_NUM_THREADS 8
}

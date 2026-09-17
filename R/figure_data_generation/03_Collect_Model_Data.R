# ==============================================================================
# STEP 03: Model Coefficients & Summary — From Power Curve Fitting Exports
# ==============================================================================
# Reads the already-exported power_equation_coefficients.csv and computes
# the model summary table with required-donors predictions.
#
# Output:
#   Paper_Source_Data/model_coefficients.csv   (copy of PCF exports, canonical)
#   Paper_Source_Data/model_summary.csv        (with Req_Donors_80pct predictions)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

if (!exists("PROJECT_ROOT")) source("config.R")

message("\n### STEP 03: MODEL COEFFICIENTS & SUMMARY ###\n")

# ==============================================================================
# 1. Copy model coefficients
# ==============================================================================

coeff_src <- file.path(PCF_BASE, "Exports", "power_equation_coefficients.csv")

if (!file.exists(coeff_src)) {
  stop("Model coefficients not found: ", coeff_src)
}

coeffs <- read.csv(coeff_src)
message(sprintf("  Read %d models from %s", nrow(coeffs), coeff_src))

# Write canonical copy
write.csv(coeffs, file.path(OUTPUT_BASE, "model_coefficients.csv"), row.names = FALSE)

# ==============================================================================
# 2. Compute model summary with required-donors predictions
# ==============================================================================

# Predict power using the logistic regression equation:
#   logit(p) = intercept + b1*log(donors) + b2*log(cells) + b3*log(rpc)
#              + b4*log(|lfc|) + b5*log(basemean) + b6*log(disp)
# Solve for donors at p = 0.80

predict_donors_for_power <- function(row, target_power = 0.80,
                                     n_cells = 250, rpc = 5000,
                                     lfc = 1.5, basemean = 50, disp = 0.5) {
  logit_target <- log(target_power / (1 - target_power))

  # logit(p) - (everything except donor term) = b1 * log(donors)
  non_donor <- row$intercept +
    row$beta_log_n_cells * log(n_cells) +
    row$beta_log_reads_per_cell * log(rpc) +
    row$beta_log_lfc * log(abs(lfc)) +
    row$beta_log_basemean * log(basemean) +
    row$beta_log_dispersion * log(disp)

  log_donors <- (logit_target - non_donor) / row$beta_log_donors
  return(round(exp(log_donors)))
}

model_summary <- coeffs %>%
  rowwise() %>%
  mutate(
    # Display name
    Model = ifelse(model_type == "combined", "Combined",
                   ifelse(grepl("Shallow_Fewer", cell_type), "Astro_Shallow_Fewer",
                          ifelse(grepl("Shallow", cell_type), "Astro_Shallow", cell_type))),
    Gene_Set = ifelse(model_type == "combined", "Combined", "Full"),
    R2 = round(pseudo_r2, 3),
    N_Obs = n_obs,
    # Moderate effect: |LFC|=1.5, baseMean=50, disp=0.5
    Req_Donors_80pct_moderate = predict_donors_for_power(
      cur_data(), lfc = 1.5, basemean = 50, disp = 0.5),
    # Large effect: |LFC|=3, baseMean=50, disp=0.5
    Req_Donors_80pct_large = predict_donors_for_power(
      cur_data(), lfc = 3, basemean = 50, disp = 0.5)
  ) %>%
  ungroup() %>%
  select(Model, Gene_Set, R2, N_Obs, Req_Donors_80pct_moderate, Req_Donors_80pct_large)

out_path <- file.path(OUTPUT_BASE, "model_summary.csv")
write.csv(model_summary, out_path, row.names = FALSE)

message(sprintf("\n[DONE] Wrote model_coefficients.csv (%d rows) and model_summary.csv (%d rows)",
                nrow(coeffs), nrow(model_summary)))

# Print summary
for (i in seq_len(nrow(model_summary))) {
  r <- model_summary[i, ]
  message(sprintf("  %s: R2=%.3f, donors@80%%: moderate=%d, large=%d",
                  r$Model, r$R2, r$Req_Donors_80pct_moderate, r$Req_Donors_80pct_large))
}

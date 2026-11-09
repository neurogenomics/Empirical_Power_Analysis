# =============================================================================
# _supp_table_index.R  --  Canonical supplementary-table ordering
# -----------------------------------------------------------------------------
# The table counterpart of _supp_fig_index.R, shared by BOTH
# `nature_comms_draft.rmd` and `nature_comms_supplementary.rmd`.
#
# TO REORDER / ADD / DROP a supplementary table:
#   1. Edit the `supp_table_keys` vector below.
#   2. Move the matching table block in the supplement to the same position
#      (the vector order MUST track the physical table order).
#   3. Re-knit. Every number in BOTH documents updates automatically.
#
# USAGE in an .Rmd:
#   ## Supplementary Table `r ST("deg_loss_equalisation")`. Title...
#   ... (Supplementary Table `r ST("deg_loss_equalisation")`) ...
#
# ST()  -> integer position of a key (errors loudly on an unknown key).
# =============================================================================

supp_table_keys <- c(
  # --- Cited from the main text, in order of first citation ------------------
  "study_survey",                     # published sc/snRNA-seq DE studies surveyed
  "deg_loss_equalisation",            # DEG loss under depth / cell-count equalisation
  "de_package_agreement",             # full-cohort agreement between DE packages
  "sex_label_validation",             # donors with discordant sex labels

  # --- Not cited from the main text ------------------------------------------
  "autosomal_vs_sexchr"               # autosomal vs sex-chromosome truth-gene properties
)

# Resolve a key to its 1-based supplementary-table number.
ST <- function(key) {
  idx <- match(key, supp_table_keys)
  if (anyNA(idx)) {
    stop(
      "ST(): unknown supplementary-table key(s): ",
      paste(key[is.na(idx)], collapse = ", "),
      "\n  Known keys are defined in _supp_table_index.R (supp_table_keys).",
      call. = FALSE
    )
  }
  idx
}

if (anyDuplicated(supp_table_keys)) {
  stop(
    "_supp_table_index.R: duplicate key(s) in supp_table_keys: ",
    paste(unique(supp_table_keys[duplicated(supp_table_keys)]), collapse = ", "),
    call. = FALSE
  )
}

# ==============================================================================
# R/shared/contrast.R -- the case-control contrast the grid is run on.
# ------------------------------------------------------------------------------
# The cartesian_downsampling and power_curve_fitting pipelines name the tested
# variable after the paper's primary contrast, sex: SEX_COL is the colData
# column, MALE_VAL the case level (log2 fold changes are MALE_VAL vs
# FEMALE_VAL) and FEMALE_VAL the reference level. The donor split and every
# downsampled cohort are stratified 50:50-proportionally on SEX_COL.
#
# EPA_CONTRAST selects which contrast those variables describe:
#   sex            (default) male vs female; covariates CATEGORICAL_VARS as set
#                  by the sourcing config (Source, ethnicity, disease)
#   ad_vs_control  AD vs control on the Roussos_ADvsCtrl shards built by
#                  R/sharding/Build_AD_Control_Shards.R (column ad_group:
#                  "AD" = AD as the only recorded diagnosis, "Control" = no
#                  recorded diagnosis). Sex becomes a covariate. Source and
#                  ethnicity are dropped: every selected donor is from one
#                  source (MSSM) and four Asian-ancestry donors are all controls.
#
# Sourced at the end of each pipeline's metadata section, so it overrides the
# sex defaults only when another contrast is requested.
# ==============================================================================

EPA_CONTRAST <- Sys.getenv("EPA_CONTRAST", "sex")

CONTRAST_SHARD_DATASETS <- c(sex = "Roussos", ad_vs_control = "Roussos_ADvsCtrl")
if (!EPA_CONTRAST %in% names(CONTRAST_SHARD_DATASETS)) {
  stop(sprintf("EPA_CONTRAST must be one of %s (got '%s')",
               paste(names(CONTRAST_SHARD_DATASETS), collapse = ", "), EPA_CONTRAST))
}

if (EPA_CONTRAST == "ad_vs_control") {
  SEX_COL          <- "ad_group"
  MALE_VAL         <- "AD"
  FEMALE_VAL       <- "Control"
  CATEGORICAL_VARS <- c("sex")
}

#' Stop if a shard dataset does not belong to EPA_CONTRAST (e.g. sex DE run on
#' the AD subset, or AD DE run on the full cohort). The sex contrast accepts any
#' dataset except the AD subset (it also runs on depth-equalised variants).
check_contrast_dataset <- function(dataset) {
  if (EPA_CONTRAST == "sex" && dataset %in% CONTRAST_SHARD_DATASETS[names(CONTRAST_SHARD_DATASETS) != "sex"]) {
    stop(sprintf("Shard dataset '%s' is for a non-sex contrast but EPA_CONTRAST='sex'", dataset))
  }
  if (EPA_CONTRAST != "sex" && !identical(dataset, CONTRAST_SHARD_DATASETS[[EPA_CONTRAST]])) {
    stop(sprintf("EPA_CONTRAST='%s' must run on shard dataset '%s' (got '%s')",
                 EPA_CONTRAST, CONTRAST_SHARD_DATASETS[[EPA_CONTRAST]], dataset))
  }
  invisible(TRUE)
}

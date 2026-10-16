#!/usr/bin/env Rscript
# ==============================================================================
# TEST: Sex balancing for the AD vs control contrast
# PURPOSE: With EPA_CONTRAST=ad_vs_control the ground-truth split is stratified
#          by AD status x sex and every subsample draws each group as close to
#          50:50 male/female as the eligible donors allow. The sex contrast must
#          sample exactly as before.
# STAGE: 1 (Unit Test, no data)
# ==============================================================================

cat("\n  UNIT TEST: AD vs control sex balancing\n\n")

Sys.setenv(EPA_CONTRAST = "ad_vs_control", SHARD_DATASET = "Roussos_ADvsCtrl")
if (!exists("PROJECT_ROOT")) {
  for (config_path in c("R/cartesian_downsampling/config.R", "../config.R", "config.R")) {
    if (file.exists(config_path)) {
      source(config_path)
      break
    }
  }
}
if (!exists("SCRIPT_DIR")) stop("Could not find config.R")
source(file.path(SCRIPT_DIR, "utils.R"))

tests_passed <- 0
tests_failed <- 0
run_test <- function(name, expr) {
  cat(sprintf("[TEST] %s... ", name))
  ok <- tryCatch(isTRUE(expr), error = function(e) {
    cat(sprintf("ERROR: %s ", e$message))
    FALSE
  })
  cat(if (ok) "PASSED\n" else "FAILED\n")
  if (ok) tests_passed <<- tests_passed + 1 else tests_failed <<- tests_failed + 1
  invisible(ok)
}

# ------------------------------------------------------------------------------
# 1. Contrast settings
# ------------------------------------------------------------------------------

run_test("AD contrast turns sex balancing on", {
  balance_sex_active() && identical(BALANCE_SEX_COL, "sex") &&
    identical(MALE_VAL, "AD") && identical(FEMALE_VAL, "Control")
})

# ------------------------------------------------------------------------------
# 2. alloc_sex_balanced()
# ------------------------------------------------------------------------------

av <- function(m, f) c(male = m, female = f)
run_test("Even n splits 50:50", identical(alloc_sex_balanced(10, av(20, 20)), av(5, 5)))
run_test("Odd n gives the extra donor to the sex with more eligible donors", {
  identical(alloc_sex_balanced(11, av(20, 30)), av(5, 6)) &&
    identical(alloc_sex_balanced(11, av(30, 20)), av(6, 5))
})
run_test("Odd n with a tie gives the extra donor to the first sex", {
  identical(alloc_sex_balanced(11, av(20, 20)), av(6, 5))
})
run_test("A short sex gives all it has and the other fills the rest", {
  identical(alloc_sex_balanced(10, av(3, 20)), av(3, 7)) &&
    identical(alloc_sex_balanced(10, av(20, 0)), av(10, 0))
})
run_test("Too few donors in the group returns NULL", is.null(alloc_sex_balanced(10, av(3, 6))))

# ------------------------------------------------------------------------------
# 3. Mock donors: AD is mostly female, controls are balanced
# ------------------------------------------------------------------------------

set.seed(42)
mk <- function(prefix, group, n_m, n_f) {
  data.frame(
    donor_id = sprintf("%s_%03d", prefix, seq_len(n_m + n_f)),
    sex = group,
    cell_count = sample(50:2000, n_m + n_f, replace = TRUE),
    total_reads = sample(1e5:1e7, n_m + n_f, replace = TRUE),
    balance_sex = rep(c("male", "female"), c(n_m, n_f)),
    stringsAsFactors = FALSE
  )
}
donor_meta <- rbind(mk("AD", "AD", 60, 120), mk("CT", "Control", 90, 90))

# ------------------------------------------------------------------------------
# 4. Ground-truth split
# ------------------------------------------------------------------------------

set.seed(1235)
gt_ids <- split_gt_donors_sex_balanced(donor_meta, 0.5)
gt_meta <- donor_meta[donor_meta$donor_id %in% gt_ids, ]
run_test("Ground truth takes half of each AD status x sex stratum", {
  identical(as.vector(table(paste(gt_meta$sex, gt_meta$balance_sex))[
    c("AD male", "AD female", "Control male", "Control female")]), c(30L, 60L, 45L, 45L))
})
run_test("Ground truth donors are unique", !anyDuplicated(gt_ids))
set.seed(1235)
run_test("Ground truth split is reproducible", identical(split_gt_donors_sex_balanced(donor_meta, 0.5), gt_ids))

val_df <- donor_meta[!donor_meta$donor_id %in% gt_ids, ]
prop_ad <- mean(donor_meta$sex == "AD")

# ------------------------------------------------------------------------------
# 5. Job plan
# ------------------------------------------------------------------------------

donor_steps <- c(5, 10, 15, 50, 100, 150)
cell_steps <- c(100, 1000, NA)
plan <- generate_job_plan(val_df, donor_steps = donor_steps, cell_steps = cell_steps,
                          mode = "cells", n_subsamples = 2, prop_male = prop_ad,
                          run_dir = tempdir(), shard_dir = "/fake")

stratum_n <- function(job, g, s) {
  sum(vapply(job$strata, function(x) if (x$group == g && x$sex == s) x$n else 0, numeric(1)))
}
eligible_n <- function(job, g, s) {
  sum(vapply(job$strata, function(x) if (x$group == g && x$sex == s) length(x$pool) else 0L, integer(1)))
}

run_test("Every job has four strata in fixed order", all(vapply(plan, function(j) {
  identical(vapply(j$strata, function(x) paste(x$group, x$sex), ""),
            c("AD male", "AD female", "Control male", "Control female"))
}, logical(1))))
run_test("Group totals match the proportional AD/control split", all(vapply(plan, function(j) {
  stratum_n(j, "AD", "male") + stratum_n(j, "AD", "female") == j$n_male &&
    stratum_n(j, "Control", "male") + stratum_n(j, "Control", "female") == j$n_female
}, logical(1))))
run_test("Each group is 50:50 (odd donor to the larger pool) unless a sex runs out", all(vapply(plan, function(j) {
  all(vapply(c("AD", "Control"), function(g) {
    n <- if (g == "AD") j$n_male else j$n_female
    avail <- c(male = eligible_n(j, g, "male"), female = eligible_n(j, g, "female"))
    got <- c(male = stratum_n(j, g, "male"), female = stratum_n(j, g, "female"))
    identical(got, as.numeric(alloc_sex_balanced(n, avail)) |> setNames(names(avail))) &&
      all(got <= avail)
  }, logical(1)))
}, logical(1))))
run_test("A short AD male pool is used in full (drift case is exercised)", any(vapply(plan, function(j) {
  stratum_n(j, "AD", "male") == eligible_n(j, "AD", "male") &&
    stratum_n(j, "AD", "male") < stratum_n(j, "AD", "female")
}, logical(1))))

key <- function(p) vapply(p, function(j) paste(j$n_donors, j$n_cells, j$iter), "")
val_nosex <- val_df[, setdiff(names(val_df), "balance_sex")]
BALANCE_SEX_COL <- NULL
plan_group_only <- generate_job_plan(val_nosex, donor_steps = donor_steps, cell_steps = cell_steps,
                                     mode = "cells", n_subsamples = 2, prop_male = prop_ad,
                                     run_dir = tempdir(), shard_dir = "/fake")
BALANCE_SEX_COL <- "sex"
run_test("Balancing keeps the same set of design points", identical(key(plan), key(plan_group_only)))

run_test("Stale validation donors without balance_sex stop the run", {
  inherits(tryCatch(generate_job_plan(val_nosex, donor_steps = 5, cell_steps = NA, mode = "cells",
                                      n_subsamples = 1, prop_male = prop_ad,
                                      run_dir = tempdir(), shard_dir = "/fake"),
                    error = function(e) e), "error")
})

# ------------------------------------------------------------------------------
# 6. Donor sampling
# ------------------------------------------------------------------------------

run_test("Sampled donors match each stratum's count and are unique", all(vapply(plan, function(j) {
  d <- sample_donors_for_job(j, seed = 99)
  m <- val_df[match(d, val_df$donor_id), ]
  length(d) == j$n_donors && !anyDuplicated(d) &&
    all(vapply(j$strata, function(x) sum(m$sex == x$group & m$balance_sex == x$sex) == x$n, logical(1)))
}, logical(1))))
j <- plan[[length(plan)]]
run_test("Same seed gives the same donors", identical(sample_donors_for_job(j, 7), sample_donors_for_job(j, 7)))
run_test("Different seeds give different donors", !identical(sample_donors_for_job(j, 7), sample_donors_for_job(j, 8)))

# ------------------------------------------------------------------------------
# 7. The sex contrast is unchanged
# ------------------------------------------------------------------------------

BALANCE_SEX_COL <- NULL
sex_val <- val_nosex
sex_val$sex <- ifelse(sex_val$sex == "AD", MALE_VAL, FEMALE_VAL)
sex_plan <- generate_job_plan(sex_val, donor_steps = donor_steps, cell_steps = cell_steps,
                              mode = "cells", n_subsamples = 2, prop_male = 0.5,
                              run_dir = tempdir(), shard_dir = "/fake")
run_test("Sex contrast jobs carry no strata", all(vapply(sex_plan, function(j) is.null(j$strata), logical(1))))
run_test("Sex contrast draws donors exactly as before", all(vapply(sex_plan, function(j) {
  set.seed(123)
  old <- c(sample(j$male_pool, j$n_male), sample(j$female_pool, j$n_female))
  identical(sample_donors_for_job(j, 123), old)
}, logical(1))))

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

cat(sprintf("\n  RESULTS: %d passed, %d failed\n\n", tests_passed, tests_failed))
quit(status = if (tests_failed > 0) 1 else 0)

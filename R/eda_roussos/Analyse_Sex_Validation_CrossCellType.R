#!/usr/bin/env Rscript
# =============================================================================
# Analyse_Sex_Validation_CrossCellType.R
# Cross-cell-type analysis of sex label discordance
#
# For each flagged donor, answers:
#   1. How severe is the discordance? (XIST/Y expression levels)
#   2. Is it seen in all three cell types, or only one?
#   3. For cell types where a donor is NOT flagged, are they present in that
#      dataset (passed thresholds) or simply absent (not enough cells)?
#
# Outputs:
#   Datalink/Sex_Validation/sex_cross_celltype_summary.csv  — one row per donor
#   Datalink/Sex_Validation/sex_cross_celltype_severity.csv — categorised flags
#   Console summary tables
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(here)
})

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

ROOT <- here::here()
DATALINK <- file.path(ROOT, "Datalink", "Sex_Validation")

CELL_TYPES <- c("Astro", "Micro", "Oligo")

# ── 1. Load all validation and exclusion tables ────────────────────────────

val <- lapply(CELL_TYPES, function(ct) {
  f <- file.path(DATALINK, paste0("sex_validation_", ct, ".csv"))
  read.csv(f, stringsAsFactors = FALSE) %>%
    mutate(cell_type = ct)
})
names(val) <- CELL_TYPES

excl <- lapply(CELL_TYPES, function(ct) {
  f <- file.path(DATALINK, paste0("sex_exclusion_", ct, ".csv"))
  read.csv(f, stringsAsFactors = FALSE) %>%
    mutate(cell_type = ct)
})
names(excl) <- CELL_TYPES

all_excl <- bind_rows(excl)
all_val  <- bind_rows(val)

flagged_donors <- unique(all_excl$donor_id)

cat("=== Flagged donors total:", length(flagged_donors), "===\n\n")

# ── 2. For each flagged donor × cell type, determine presence + flag status ──

cross <- lapply(flagged_donors, function(d) {
  lapply(CELL_TYPES, function(ct) {
    in_val  <- filter(val[[ct]], donor_id == d)
    in_excl <- filter(excl[[ct]], donor_id == d)

    present  <- nrow(in_val) > 0
    flagged  <- nrow(in_excl) > 0

    xist_pct  <- if (present) in_val$XIST_pct_cells  else NA_real_
    xist_mean <- if (present) in_val$XIST_mean_expr   else NA_real_
    n_cells   <- if (present) in_val$n_cells           else NA_integer_
    flag_text <- if (flagged) in_excl$flag             else if (present) "PASS" else "ABSENT"
    sex_label <- if (present) in_val$sex_label[1]      else NA_character_

    data.frame(
      donor_id   = d,
      cell_type  = ct,
      sex_label  = sex_label,
      present    = present,
      flag       = flag_text,
      XIST_pct   = xist_pct,
      XIST_mean  = xist_mean,
      n_cells    = n_cells,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
}) %>% bind_rows()

# ── 3. Wide summary: one row per donor ────────────────────────────────────

wide <- cross %>%
  pivot_wider(
    id_cols     = c(donor_id, sex_label),
    names_from  = cell_type,
    values_from = c(present, flag, XIST_pct, XIST_mean, n_cells)
  ) %>%
  # Count how many cell types each donor is flagged in
  mutate(
    n_flagged_ct = rowSums(
      data.frame(
        a = flag_Astro %in% c("XIST_present_in_male", "Y_present_in_female",
                               "XIST_present_in_male | Y_absent_in_male"),
        m = flag_Micro %in% c("XIST_present_in_male", "Y_present_in_female",
                               "XIST_present_in_male | Y_absent_in_male"),
        o = flag_Oligo %in% c("XIST_present_in_male", "Y_present_in_female",
                               "XIST_present_in_male | Y_absent_in_male")
      ), na.rm = TRUE
    ),
    n_present_ct = rowSums(
      data.frame(
        a = !is.na(present_Astro) & present_Astro,
        m = !is.na(present_Micro) & present_Micro,
        o = !is.na(present_Oligo) & present_Oligo
      ), na.rm = TRUE
    )
  ) %>%
  arrange(desc(n_flagged_ct), donor_id)

# ── 4. Console output ─────────────────────────────────────────────────────

cat("=== CROSS-CELL-TYPE FLAG CONCORDANCE ===\n\n")
cat(sprintf("Flagged in ALL 3 cell types:   %d donors\n",
            sum(wide$n_flagged_ct == 3)))
cat(sprintf("Flagged in 2 cell types only:  %d donors\n",
            sum(wide$n_flagged_ct == 2)))
cat(sprintf("Flagged in 1 cell type only:   %d donors\n\n",
            sum(wide$n_flagged_ct == 1)))

# Show the 1-cell-type-only donors with their presence status in other cell types
cat("=== DONORS FLAGGED IN ONLY 1 CELL TYPE — presence in others? ===\n\n")
single <- wide %>% filter(n_flagged_ct == 1)
for (i in seq_len(nrow(single))) {
  d <- single[i, ]
  astro_status <- if (!d$present_Astro) "ABSENT (no shard)" else paste0("PASS (XIST=", round(d$XIST_pct_Astro,1), "%)")
  micro_status <- if (!d$present_Micro) "ABSENT (no shard)" else paste0("PASS (XIST=", round(d$XIST_pct_Micro,1), "%)")
  oligo_status <- if (!d$present_Oligo) "ABSENT (no shard)" else paste0("PASS (XIST=", round(d$XIST_pct_Oligo,1), "%)")

  flagged_in <- c(
    if (grepl("XIST|Y_present", d$flag_Astro %||% "")) "Astro",
    if (grepl("XIST|Y_present", d$flag_Micro %||% "")) "Micro",
    if (grepl("XIST|Y_present", d$flag_Oligo %||% "")) "Oligo"
  )

  cat(sprintf("  %s (%s) — flagged in %s\n",
              d$donor_id, d$sex_label, paste(flagged_in, collapse=", ")))
  cat(sprintf("    Astro: %s | Micro: %s | Oligo: %s\n",
              if (grepl("XIST|Y_present", d$flag_Astro %||% ""))
                paste0("FLAGGED (XIST=", round(d$XIST_pct_Astro,1), "%)") else astro_status,
              if (grepl("XIST|Y_present", d$flag_Micro %||% ""))
                paste0("FLAGGED (XIST=", round(d$XIST_pct_Micro,1), "%)") else micro_status,
              if (grepl("XIST|Y_present", d$flag_Oligo %||% ""))
                paste0("FLAGGED (XIST=", round(d$XIST_pct_Oligo,1), "%)") else oligo_status
  ))
}

# ── 5. Severity classification ─────────────────────────────────────────────

cat("\n\n=== SEVERITY CLASSIFICATION (for donors in flagged cell types) ===\n\n")

severity <- all_excl %>%
  mutate(
    severity = case_when(
      # Female-labelled males / male-labelled females: strong signal
      XIST_pct_cells >= 80 ~ "EXTREME (>=80% XIST)",
      XIST_pct_cells >= 40 ~ "HIGH (40-80% XIST)",
      XIST_pct_cells >= 20 ~ "MODERATE (20-40% XIST)",
      TRUE ~ "BORDERLINE (<20% XIST)"
    )
  ) %>%
  select(donor_id, cell_type, sex_label, flag, XIST_pct_cells, XIST_mean_expr, n_cells, severity)

cat("Severity breakdown:\n")
print(table(severity$severity, severity$cell_type))

cat("\nEXTREME cases (>=80% XIST):\n")
extreme <- severity %>% filter(grepl("EXTREME", severity)) %>%
  arrange(cell_type, donor_id)
print(extreme[, c("donor_id","cell_type","sex_label","XIST_pct_cells","n_cells")])

cat("\nBORDERLINE cases (<20% XIST):\n")
borderline <- severity %>% filter(grepl("BORDERLINE", severity))
print(borderline[, c("donor_id","cell_type","sex_label","XIST_pct_cells","n_cells")])

# ── 6. Save outputs ───────────────────────────────────────────────────────

out_cross  <- file.path(DATALINK, "sex_cross_celltype_summary.csv")
out_sev    <- file.path(DATALINK, "sex_cross_celltype_severity.csv")

write.csv(wide,     out_cross, row.names = FALSE)
write.csv(severity, out_sev,   row.names = FALSE)

cat("\nOutputs written:\n")
cat("  ", out_cross, "\n")
cat("  ", out_sev,   "\n")
cat("\nDone.\n")

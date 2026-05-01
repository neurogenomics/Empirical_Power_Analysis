# In R/shared/donor_id_utils.R (new file)
roussos_pad_donor <- function(x) {
  m <- regmatches(x, regexec("^Donor_(\\d+)$", x))
  vapply(seq_along(x), function(i) {
    if (length(m[[i]]) < 2) return(NA_character_)
    sprintf("Donor_%04d", as.integer(m[[i]][2]))
  }, character(1))
}

#!/usr/bin/env Rscript
# =============================================================================
#  Hardened renv restore for Docker builds.
#
#  Public-internet R package installs are flaky. A single transient 4xx/5xx
#  from a CRAN mirror (e.g. Rcpp fails to download) cascades into every
#  package that depends on it (plyr, RcppArmadillo, RcppEigen, DESeq2, ...),
#  and renv::restore() gives up on that whole set without retrying.
#
#  This script is the defensive wrapper around renv::restore. It:
#    1. Configures multiple mirrors with binary (Posit P3M) preferred.
#    2. Sets aggressive retry/timeout options on renv's downloaders.
#    3. Pre-installs compilation-critical packages (Rcpp family) so their
#       failure can't cascade.
#    4. Retries renv::restore() up to 3 times — between attempts, any
#       package left missing gets another attempt with a different mirror
#       order.
#    5. After restore, computes the still-missing set from renv.lock and
#       brute-forces it via renv::install. Any leftover aborts the build.
# =============================================================================

# ---- retry helper -----------------------------------------------------------
with_retry <- function(expr, attempts = 3, what = "operation") {
  for (i in seq_len(attempts)) {
    res <- tryCatch(
      { force(expr); return(invisible(TRUE)) },
      error = function(e) {
        cat(sprintf("[retry] %s attempt %d/%d failed: %s\n",
                    what, i, attempts, conditionMessage(e)))
        if (i < attempts) Sys.sleep(5 * i)  # backoff: 5s, 10s
        e
      }
    )
  }
  stop(sprintf("%s failed after %d attempts", what, attempts))
}

# ---- mirror configuration ---------------------------------------------------
# Posit Public Package Manager serves precompiled Ubuntu 24.04 ("noble")
# binaries for CRAN packages. Using it as the primary repo drops build time
# from ~15 min to ~2 min AND removes most compile-error failure modes
# (headers missing, etc.).
#
# We query P3M's "latest" snapshot. For deterministic builds one could pin a
# date, but then new Bioconductor releases fall out of sync with CRAN — the
# lockfile already pins versions, so the mirror just needs to serve them.
p3m_noble <- "https://p3m.dev/cran/__linux__/noble/latest"
repos <- c(
  P3M           = p3m_noble,
  CRAN          = "https://cloud.r-project.org",
  `CRAN-UK`     = "https://www.stats.bris.ac.uk/R/",
  BioC          = "https://bioconductor.org/packages/3.21/bioc",
  `BioC-Ann`    = "https://bioconductor.org/packages/3.21/data/annotation",
  `BioC-Exp`    = "https://bioconductor.org/packages/3.21/data/experiment"
)
options(repos = repos)

# HTTPUserAgent must advertise an R version + distro for P3M to serve Linux
# binaries. Bioconductor_docker base image already sets something sensible,
# but be explicit to guarantee binary downloads work.
options(HTTPUserAgent = sprintf(
  "R/%s R (%s) on %s",
  getRversion(),
  paste(getRversion(), "x86_64-pc-linux-gnu", "x86_64", "linux-gnu"),
  "Ubuntu 24.04"
))

# ---- renv tuning ------------------------------------------------------------
options(
  renv.config.pak.enabled       = FALSE,  # pak has its own quirks; stick with renv's installer
  renv.config.connect.timeout   = 60,     # seconds
  renv.config.connect.retry     = 5,      # per-download retry
  renv.config.install.transactional = FALSE,  # don't roll back the whole batch on one failure
  renv.config.updates.check     = FALSE,
  renv.download.override        = NULL,
  Ncpus                         = max(1L, parallel::detectCores() - 1L),
  timeout                       = 600     # generous download timeout
)

# Ensure renv is loaded.
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = repos[["CRAN"]])
}

# ---- sanity check network ---------------------------------------------------
cat("=== Network reachability probe ===\n")
for (name in names(repos)) {
  url <- paste0(repos[[name]], "/src/contrib/PACKAGES")
  ok <- tryCatch({
    con <- url(url, open = "rb")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    readBin(con, "raw", n = 128L)
    TRUE
  }, error = function(e) FALSE)
  cat(sprintf("  [%s] %s %s\n", ifelse(ok, "OK ", "XX "), name, url))
}

# ---- pre-install compilation-critical packages ------------------------------
# Rcpp + friends are the most common cascade triggers. Install them FIRST
# (and retry individually) so a transient failure doesn't poison the rest
# of renv::restore. Using lockfile versions via renv::install().
critical_build_pkgs <- c("Rcpp", "RcppArmadillo", "RcppEigen", "RcppTOML")
cat("\n=== Pre-installing compilation-critical packages ===\n")
for (pkg in critical_build_pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  [skip] %s already installed (%s)\n",
                pkg, packageVersion(pkg)))
    next
  }
  with_retry(
    renv::install(pkg, prompt = FALSE),
    attempts = 4,
    what = sprintf("install %s", pkg)
  )
}

# ---- renv::restore with retries --------------------------------------------
cat("\n=== renv::restore ===\n")
restore_attempts <- 3L
for (i in seq_len(restore_attempts)) {
  cat(sprintf("\n--- restore attempt %d/%d ---\n", i, restore_attempts))
  result <- tryCatch(
    {
      renv::restore(prompt = FALSE)
      TRUE
    },
    error = function(e) {
      cat(sprintf("restore attempt %d errored: %s\n", i, conditionMessage(e)))
      FALSE
    }
  )

  # Check what's still missing against the lockfile.
  lock <- renv::lockfile_read()
  lock_pkgs <- names(lock$Packages)
  installed <- rownames(installed.packages())
  missing <- setdiff(lock_pkgs, installed)

  if (length(missing) == 0L) {
    cat("All lockfile packages installed.\n")
    break
  }

  cat(sprintf("%d package(s) still missing after attempt %d:\n  %s\n",
              length(missing), i, paste(missing, collapse = ", ")))

  if (i < restore_attempts) {
    # Rotate the mirror order so the next attempt hits a different primary.
    repos <- repos[sample(seq_along(repos))]
    options(repos = repos)
    Sys.sleep(10)
  }
}

# ---- final force-install pass ----------------------------------------------
lock <- renv::lockfile_read()
lock_pkgs <- names(lock$Packages)
installed <- rownames(installed.packages())
missing <- setdiff(lock_pkgs, installed)

if (length(missing) > 0L) {
  cat(sprintf("\n=== Force-install pass for %d remaining package(s) ===\n",
              length(missing)))
  # Install one-at-a-time with retries so one failure doesn't stop the rest.
  stubborn <- character()
  for (pkg in missing) {
    ok <- tryCatch({
      with_retry(
        renv::install(pkg, prompt = FALSE),
        attempts = 3,
        what = sprintf("force-install %s", pkg)
      )
      TRUE
    }, error = function(e) FALSE)
    if (!ok) stubborn <- c(stubborn, pkg)
  }
  if (length(stubborn) > 0L) {
    stop(sprintf(
      "Unrecoverable: %d package(s) could not be installed after all retries:\n  %s\nCheck network egress from the build host (proxy, firewall, DNS).",
      length(stubborn), paste(stubborn, collapse = ", ")
    ))
  }
}

# ---- final status ----------------------------------------------------------
cat("\n=== renv::status() ===\n")
print(renv::status())

installed <- rownames(installed.packages())
missing <- setdiff(lock_pkgs, installed)
if (length(missing) > 0L) {
  stop("Lockfile still has missing packages after all passes: ",
       paste(missing, collapse = ", "))
}
cat(sprintf("\nAll %d lockfile packages installed.\n", length(lock_pkgs)))

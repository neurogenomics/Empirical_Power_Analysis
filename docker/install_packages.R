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

# Activate the renv project so installed.packages() reflects the project
# library only, NOT the bioconductor base image's preinstalled site library.
# Without this, renv considers DESeq2 etc. already-installed and skips them,
# leaving /opt/renv-library empty.
source("renv/activate.R")

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
# binaries for CRAN packages, used here as the primary repo.
#
# NOTE on build speed: renv::restore() downloads SOURCE packages even when P3M
# has a matching binary (its retrieval path differs from install.packages(),
# which does content-negotiate for binaries). So in practice most CRAN packages
# compile from source during restore regardless of this mirror, and all Bioc
# packages compile too. Pinning a dated P3M snapshot was tested and made no
# difference to restore (verified 2026-06-05). Speeding up restore would require
# making renv itself prefer binaries — left as a future optimisation. Versions
# are unaffected either way: restore always installs the exact lockfile pins.
p3m_noble <- "https://p3m.dev/cran/__linux__/noble/latest"
repos <- c(
  P3M             = p3m_noble,
  CRAN            = "https://cloud.r-project.org",
  `CRAN-UK`       = "https://www.stats.bris.ac.uk/R/",
  BioC            = "https://bioconductor.org/packages/3.22/bioc",
  `BioC-Ann`      = "https://bioconductor.org/packages/3.22/data/annotation",
  `BioC-Exp`      = "https://bioconductor.org/packages/3.22/data/experiment"
)
options(repos = repos)

# Track packages whose installed version drifts from the lockfile pin, so the
# build can surface "we couldn't honour the strict pin for X, used current
# CRAN instead" rather than silently shipping a different version.
deviation_log <- "/tmp/lockfile-deviations.txt"
record_deviation <- function(pkg, reason) {
  line <- sprintf("%s\t%s\t%s\n",
                  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                  pkg, reason)
  cat(line, file = deviation_log, append = TRUE)
}

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
reachable <- logical(length(repos))
names(reachable) <- names(repos)
for (name in names(repos)) {
  url <- paste0(repos[[name]], "/src/contrib/PACKAGES")
  ok <- tryCatch({
    con <- url(url, open = "rb")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    readBin(con, "raw", n = 128L)
    TRUE
  }, error = function(e) FALSE)
  reachable[[name]] <- ok
  cat(sprintf("  [%s] %s %s\n", ifelse(ok, "OK ", "XX "), name, url))
}

# Drop unreachable repos from the live list so renv doesn't waste 6+ minutes
# per package timing out against a dead mirror (common: P3M outages,
# institutional proxies blocking p3m.dev). MUST keep at least CRAN reachable.
if (any(!reachable)) {
  cat("\n[prune] dropping unreachable repos from active list:\n")
  for (name in names(reachable)[!reachable]) {
    cat(sprintf("  - %s (%s)\n", name, repos[[name]]))
  }
  repos <- repos[reachable]
  if (length(repos) == 0L || !any(grepl("cran", names(repos), ignore.case = TRUE))) {
    stop("No reachable CRAN-style repo left after probe — aborting build.")
  }
  options(repos = repos)
}

# Add CRAN Archive as an additional source-only fallback so packages whose
# pinned version has been superseded on CRAN main (e.g. qs 0.27.2) can still
# be fetched. renv won't auto-look here, so this only helps the loose
# install.packages() fallback at the bottom of the script.
cran_archive <- "https://cran.r-project.org/src/contrib/Archive"

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
    locked_version <- tryCatch(lock$Packages[[pkg]]$Version, error = function(e) NA_character_)

    # Tier 1 — strict pin via renv::install (honours lockfile version).
    ok <- tryCatch({
      with_retry(
        renv::install(pkg, prompt = FALSE),
        attempts = 3,
        what = sprintf("force-install %s", pkg)
      )
      TRUE
    }, error = function(e) FALSE)

    # Tier 2 — drop the pin: take whatever version is currently on CRAN.
    # This breaks strict reproducibility for the affected package but is
    # better than failing the whole build. Logged to deviation_log.
    if (!ok) {
      cat(sprintf("\n[loose-fallback] strict install of %s failed; trying current CRAN version\n", pkg))
      ok <- tryCatch({
        install.packages(pkg, repos = "https://cloud.r-project.org",
                         dependencies = TRUE, Ncpus = max(1L, parallel::detectCores() - 1L))
        if (requireNamespace(pkg, quietly = TRUE)) {
          new_version <- as.character(packageVersion(pkg))
          if (!is.na(locked_version) && new_version != locked_version) {
            record_deviation(pkg, sprintf("strict-pin %s unavailable; installed current CRAN %s",
                                          locked_version, new_version))
            cat(sprintf("[loose-fallback] %s: lockfile pinned %s, installed %s\n",
                        pkg, locked_version, new_version))
          }
          TRUE
        } else FALSE
      }, error = function(e) {
        cat(sprintf("[loose-fallback] install.packages(%s) errored: %s\n", pkg, conditionMessage(e)))
        FALSE
      })
    }

    # Tier 3 — try CRAN Archive directly for the pinned version.
    # Last-ditch: useful when current CRAN version has a hard ABI break.
    if (!ok && !is.na(locked_version)) {
      archive_url <- sprintf("%s/%s/%s_%s.tar.gz", cran_archive, pkg, pkg, locked_version)
      cat(sprintf("[archive-fallback] trying %s\n", archive_url))
      ok <- tryCatch({
        install.packages(archive_url, repos = NULL, type = "source",
                         dependencies = TRUE, Ncpus = max(1L, parallel::detectCores() - 1L))
        requireNamespace(pkg, quietly = TRUE)
      }, error = function(e) {
        cat(sprintf("[archive-fallback] %s\n", conditionMessage(e)))
        FALSE
      })
      if (ok) record_deviation(pkg, sprintf("installed pin %s from CRAN Archive (main mirror missing it)", locked_version))
    }

    if (!ok) stubborn <- c(stubborn, pkg)
  }
  if (length(stubborn) > 0L) {
    stop(sprintf(
      "Unrecoverable: %d package(s) could not be installed after all 3 fallback tiers:\n  %s\nCheck network egress from the build host (proxy, firewall, DNS) and that the package names are correct.",
      length(stubborn), paste(stubborn, collapse = ", ")
    ))
  }
}

# ---- surface lockfile deviations ------------------------------------------
if (file.exists(deviation_log)) {
  cat("\n=== Lockfile deviations during this build ===\n")
  cat(readLines(deviation_log), sep = "\n")
  cat("\n(See /tmp/lockfile-deviations.txt inside the image.)\n")
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

# ---- bake the renv library path into R's site Renviron ---------------------
# A plain `R` / `Rscript` session does NOT add the renv project library to
# .libPaths() on its own; that normally requires renv's autoloader, which we
# deliberately disable inside the container (see docker/entrypoint.sh) for
# determinism. renv installs into a platform/version-specific leaf under
# RENV_PATHS_LIBRARY (e.g. /opt/renv-library/linux-ubuntu-noble/R-4.5/...), so
# pointing .libPaths() at /opt/renv-library itself finds nothing. Instead we
# record the *resolved* leaf in R's site Renviron, which R reads natively at
# startup, making the pinned packages visible to every R session: the
# build-time sanity checks, the smoke test, and the pipeline at runtime. We
# keep the default site-library on the path too (glmGamPoi installs there).
renv_lib <- grep("/opt/renv-library", .libPaths(), value = TRUE, fixed = TRUE)[1]
if (is.na(renv_lib)) {
  stop("Could not locate the /opt/renv-library leaf on .libPaths(); got: ",
       paste(.libPaths(), collapse = ", "))
}
site_lib <- file.path(R.home(), "site-library")
renviron_site <- file.path(R.home("etc"), "Renviron.site")
cat(sprintf("R_LIBS_SITE=%s:%s\n", renv_lib, site_lib),
    file = renviron_site, append = TRUE)
cat(sprintf("Baked R_LIBS_SITE into %s:\n  %s:%s\n",
            renviron_site, renv_lib, site_lib))

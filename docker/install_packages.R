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

# Ensure the Dockerfile-provided renv bootstrap is the exact expected version.
expected_renv <- Sys.getenv("RENV_VERSION", unset = "1.2.2")
if (!requireNamespace("renv", quietly = TRUE)) {
  stop("renv is not installed. The Dockerfile must install the pinned renv tarball before running this script.")
}
actual_renv <- as.character(packageVersion("renv"))
if (!identical(actual_renv, expected_renv)) {
  stop(sprintf("renv version mismatch: expected %s, loaded %s", expected_renv, actual_renv))
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
# be fetched without dropping the lockfile version.
cran_archive <- "https://cran.r-project.org/src/contrib/Archive"

# ---- pre-install compilation-critical packages ------------------------------
# Rcpp + friends are the most common cascade triggers. Install them FIRST
# (and retry individually) so a transient failure doesn't poison the rest
# of renv::restore. A bare renv::install(pkg) installs the LATEST version, so
# each one is requested as pkg@<lockfile version>.
lock <- renv::lockfile_read()
critical_build_pkgs <- c("Rcpp", "RcppArmadillo", "RcppEigen", "RcppTOML")
cat("\n=== Pre-installing compilation-critical packages ===\n")
for (pkg in critical_build_pkgs) {
  pin <- lock$Packages[[pkg]]$Version
  if (requireNamespace(pkg, quietly = TRUE) &&
      identical(as.character(packageVersion(pkg)), pin)) {
    cat(sprintf("  [skip] %s %s already installed\n", pkg, pin))
    next
  }
  with_retry(
    renv::install(sprintf("%s@%s", pkg, pin), prompt = FALSE),
    attempts = 4,
    what = sprintf("install %s@%s", pkg, pin)
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
        renv::install(sprintf("%s@%s", pkg, locked_version), prompt = FALSE),
        attempts = 3,
        what = sprintf("force-install %s", pkg)
      )
      TRUE
    }, error = function(e) FALSE)

    # Tier 2: try CRAN Archive directly for the pinned version.
    # Last ditch: useful when current CRAN main no longer carries the pin.
    # Never drop the lockfile version: this Docker image is the paper's
    # reproducibility boundary.
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
    }

    if (!ok) stubborn <- c(stubborn, pkg)
  }
  if (length(stubborn) > 0L) {
    stop(sprintf(
      "Unrecoverable: %d package(s) could not be installed after strict install and archive fallback:\n  %s\nCheck network egress from the build host (proxy, firewall, DNS) and that the package names are correct.",
      length(stubborn), paste(stubborn, collapse = ", ")
    ))
  }
}

# ---- surface lockfile deviations ------------------------------------------
if (file.exists(deviation_log)) {
  cat("\n=== Lockfile deviations during this build ===\n")
  cat(readLines(deviation_log), sep = "\n")
  stop("\nLockfile deviations are not permitted for the paper Docker image.")
}

# ---- final status ----------------------------------------------------------
cat("\n=== renv::status() ===\n")
print(renv::status())

# Presence is not enough: a package pulled in as a dependency, or installed
# before restore, can sit at a different version from its pin. Compare every
# lockfile package with the version in the project library, reinstall any
# mismatch at its pinned version once, then refuse to build if any remain.
#
# The repair uses base install.packages() on the exact source tarball, NOT
# renv::install(): renv::install("pkg@ver") re-resolves the package's
# dependencies and upgrades them to the newest version it can find (tested
# 2026-10-14: installing stringfish@0.16.0 replaced the pinned RcppParallel
# 5.1.11-2 with 6.2.1 from the renv cache). install.packages(repos = NULL)
# installs exactly one package and never touches its dependencies.
install_exact <- function(pkg, version, repository = "CRAN") {
  urls <- if (identical(repository, "Bioconductor")) {
    sprintf("https://bioconductor.org/packages/3.22/bioc/src/contrib/%s_%s.tar.gz", pkg, version)
  } else {
    c(sprintf("https://cloud.r-project.org/src/contrib/%s_%s.tar.gz", pkg, version),
      sprintf("%s/%s/%s_%s.tar.gz", cran_archive, pkg, pkg, version))
  }
  tarball <- file.path(tempdir(), sprintf("%s_%s.tar.gz", pkg, version))
  ok <- FALSE
  for (u in urls) {
    ok <- tryCatch(utils::download.file(u, tarball, quiet = TRUE, mode = "wb") == 0L,
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (ok) break
  }
  if (!ok) stop(sprintf("no source tarball found for %s %s", pkg, version))
  utils::install.packages(tarball, repos = NULL, type = "source",
                          lib = .libPaths()[1])
}
lock_repo <- function(pkg) {
  if (identical(lock$Packages[[pkg]]$Source, "Bioconductor")) "Bioconductor" else "CRAN"
}
version_mismatches <- function() {
  ip <- installed.packages(lib.loc = .libPaths()[1])[, "Version"]
  want <- vapply(lock$Packages, function(x) x$Version, character(1))
  have <- ip[names(want)]
  names(want)[is.na(have) | have != want]
}
bad <- version_mismatches()
for (pkg in bad) {
  pin <- lock$Packages[[pkg]]$Version
  cat(sprintf("[pin] reinstalling %s at %s\n", pkg, pin))
  try(with_retry(install_exact(pkg, pin, lock_repo(pkg)),
                 attempts = 3, what = sprintf("pin %s@%s", pkg, pin)))
}
# Packages that compile against a repaired package (LinkingTo) were built
# against the wrong headers or libraries (e.g. stringfish against TBB from
# RcppParallel 6.x fails to load once 5.1.11-2 is restored), so rebuild them.
if (length(bad) > 0L) {
  ip_all <- installed.packages(lib.loc = .libPaths()[1])
  links <- tools:::.extract_dependency_package_names
  rebuild <- rownames(ip_all)[vapply(ip_all[, "LinkingTo"], function(x)
    any(links(x) %in% bad), logical(1))]
  rebuild <- setdiff(intersect(rebuild, names(lock$Packages)), bad)
  for (pkg in rebuild) {
    pin <- lock$Packages[[pkg]]$Version
    cat(sprintf("[pin] rebuilding %s %s against repaired dependencies\n", pkg, pin))
    try(with_retry(install_exact(pkg, pin, lock_repo(pkg)),
                   attempts = 3, what = sprintf("rebuild %s@%s", pkg, pin)))
  }
}
bad <- version_mismatches()
if (length(bad) > 0L) {
  ip <- installed.packages(lib.loc = .libPaths()[1])[, "Version"]
  stop("Installed versions differ from renv.lock: ",
       paste(sprintf("%s (installed %s, locked %s)", bad, ip[bad],
                     vapply(lock$Packages[bad], function(x) x$Version, character(1))),
             collapse = "; "))
}
cat(sprintf("\nAll %d lockfile packages installed at their pinned versions.\n",
            length(lock_pkgs)))

# ---- bake the renv library path into R's site Renviron ---------------------
# A plain `R` / `Rscript` session does NOT add the renv project library to
# .libPaths() on its own; that normally requires renv's autoloader, which we
# deliberately disable inside the container (see docker/entrypoint.sh) for
# determinism. renv installs into a platform/version-specific leaf under
# RENV_PATHS_LIBRARY (e.g. /opt/renv-library/linux-ubuntu-noble/R-4.5/...), so
# pointing .libPaths() at /opt/renv-library itself finds nothing. Instead we
# record the *resolved* leaf in R's site Renviron, which R reads natively at
# startup, making the pinned packages visible to every R session: the
# build-time sanity checks, the smoke test, and the pipeline at runtime.
# Do not add the default site-library here: packages used by the paper must
# come from the lockfile-restored renv library, with base/recommended packages
# still available through R's base library.
renv_lib <- grep("/opt/renv-library", .libPaths(), value = TRUE, fixed = TRUE)[1]
if (is.na(renv_lib)) {
  stop("Could not locate the /opt/renv-library leaf on .libPaths(); got: ",
       paste(.libPaths(), collapse = ", "))
}
renviron_site <- file.path(R.home("etc"), "Renviron.site")
cat(sprintf("R_LIBS_SITE=%s\n", renv_lib), file = renviron_site, append = TRUE)
cat(sprintf("Baked R_LIBS_SITE into %s:\n  %s\n",
            renviron_site, renv_lib))

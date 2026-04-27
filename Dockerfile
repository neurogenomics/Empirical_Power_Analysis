# syntax=docker/dockerfile:1.6
# =============================================================================
#  Empirical Power Analysis — Reproducible Compute Environment
#  Base: Bioconductor 3.21 (R 4.5.1) on Ubuntu 24.04
#  Build: docker build -t ghcr.io/neurogenomics/empirical_power_analysis:1.0.0 .
#  Push:  docker push ghcr.io/neurogenomics/empirical_power_analysis:1.0.0
# =============================================================================

FROM bioconductor/bioconductor_docker:RELEASE_3_21

LABEL org.opencontainers.image.title="Empirical Power Analysis"
LABEL org.opencontainers.image.description="Reproducible environment for empirical scRNA-seq power analysis"
LABEL org.opencontainers.image.source="https://github.com/neurogenomics/Empirical_Power_Analysis"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.authors="Skene Lab, UK DRI Imperial"

# -----------------------------------------------------------------------------
#  System dependencies beyond bioconductor_docker defaults
# -----------------------------------------------------------------------------
# pandoc-citeproc was merged into pandoc itself (pandoc >= 2.11) and removed
# from Ubuntu 22.04+ repos. The Bioconductor 3.21 base is Ubuntu 24.04, so
# we use `pandoc --citeproc` instead of the deprecated separate binary.
# APT cache mounts keep downloaded .deb packages between builds, so a rebuild
# of this layer skips the re-download of ~1GB of texlive. Requires BuildKit
# (DOCKER_BUILDKIT=1 or the syntax directive at the top of this file).
# rm /etc/apt/apt.conf.d/docker-clean is needed or the bioconductor base image
# will wipe the cache after install.
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    rm -f /etc/apt/apt.conf.d/docker-clean \
 && apt-get update && apt-get install -y --no-install-recommends \
      pandoc \
      texlive-xetex \
      texlive-fonts-recommended \
      texlive-plain-generic \
      ccache \
      libglpk40 \
      jq

# -----------------------------------------------------------------------------
#  R package installation via renv — reproducible snapshot
# -----------------------------------------------------------------------------
# Copy renv artefacts first to maximise Docker layer caching
COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

# Install renv itself (rarely changes; gets its own cacheable layer).
# Cache mount on R's site-library means a rebuild of just this layer doesn't
# re-download renv from CRAN.
RUN --mount=type=cache,target=/root/.cache/R,sharing=locked \
    R -e "install.packages('renv', repos = 'https://cloud.r-project.org')"

# -----------------------------------------------------------------------------
#  Copy the hardened restore script. Keeping this in its own file means we
#  can iterate on restore logic without busting the APT/renv-install layers
#  above, and avoids shell-quoting nightmares inside a RUN.
# -----------------------------------------------------------------------------
COPY docker/install_packages.R /tmp/install_packages.R

# Restore the pinned environment with a hardened, retry-heavy policy.
#
# Why this layer is complicated: R package installation over the public
# internet is flaky — transient 4xx/5xx from CRAN mirrors, half-truncated
# downloads, DNS blips — and a single failure cascades (Rcpp fails →
# 18 dependent packages also fail). We mitigate with:
#
#   1. Posit Public Package Manager (P3M) binary repo for Ubuntu 24.04
#      ("noble"), which serves precompiled binaries. Drops total build time
#      from ~15 min of source compilation to ~2 min, AND each download is
#      one binary tarball instead of source + configure + make.
#   2. Explicit multi-mirror CRAN fallback (cloud + uk + bioc-mirror)
#      in case P3M is down or doesn't carry a package.
#   3. renv connection options: 60s timeout, 5 retries per download.
#   4. Outer retry loop around renv::restore() (3 attempts) — between
#      attempts, rebuild the "missing" list and try again.
#   5. Final pass: force-install anything still missing via renv::install
#      with the full mirror list; error hard if the set is non-empty.
#   6. Dev-tools (Rcpp, RcppArmadillo, RcppEigen, RcppTOML) installed FIRST
#      so their failure can't cascade into 18 dependents in the same attempt.
#
# Cache mounts (unchanged):
#   - /root/.cache/R/renv         : renv's own package cache (compiled binaries)
#   - /root/.cache/R              : CRAN download cache (tarballs)
#   - /tmp/ccache                 : C/C++ compilation cache
#
# RENV_CONFIG_CACHE_SYMLINKS=FALSE: the renv cache lives on a BuildKit cache
# mount which is NOT persisted into the image layer. With symlinks disabled,
# renv copies packages directly into the project library so they persist.
RUN --mount=type=cache,target=/root/.cache/R/renv,sharing=locked \
    --mount=type=cache,target=/root/.cache/R,sharing=locked \
    --mount=type=cache,target=/tmp/ccache,sharing=locked \
    mkdir -p ~/.R \
 && echo "CC=ccache gcc"   >> ~/.R/Makevars \
 && echo "CXX=ccache g++"  >> ~/.R/Makevars \
 && echo "CXX11=ccache g++" >> ~/.R/Makevars \
 && echo "CXX14=ccache g++" >> ~/.R/Makevars \
 && echo "CXX17=ccache g++" >> ~/.R/Makevars \
 && export CCACHE_DIR=/tmp/ccache \
 && export RENV_CONFIG_CACHE_SYMLINKS=FALSE \
 && Rscript /tmp/install_packages.R

# -----------------------------------------------------------------------------
#  Make renv library + rstudio home world-readable so non-root runtime users
#  (Nextflow invokes `docker run -u $(id -u):$(id -g)`) can actually load the
#  packages. Without this the image works only when executed as root.
# -----------------------------------------------------------------------------
RUN chmod -R a+rX /home/rstudio /usr/local/lib/R 2>/dev/null || true

# -----------------------------------------------------------------------------
#  Verify critical packages loaded
# -----------------------------------------------------------------------------
RUN R -e " \
    critical <- c('DESeq2', 'apeglm', 'SingleCellExperiment', 'biomaRt', \
                  'zellkonverter', 'reticulate', 'basilisk', 'qs', \
                  'Rcpp', 'RcppArmadillo', 'RcppEigen', \
                  'future.apply', 'rmarkdown', 'ggplot2', 'dplyr'); \
    missing <- character(); \
    for (p in critical) { \
        if (requireNamespace(p, quietly = TRUE)) { \
            cat(sprintf('OK %s %s\n', p, packageVersion(p))); \
        } else { \
            missing <- c(missing, p); \
            cat(sprintf('FAIL %s\n', p)); \
        } \
    }; \
    if (length(missing) > 0) stop('critical packages missing: ', paste(missing, collapse=', '))"

# -----------------------------------------------------------------------------
#  Default working directory matches POWER_ANALYSIS_ROOT convention
# -----------------------------------------------------------------------------
WORKDIR /workspace
ENV POWER_ANALYSIS_ROOT=/workspace

# -----------------------------------------------------------------------------
#  Sanity checks at container startup
# -----------------------------------------------------------------------------
CMD ["R"]

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
RUN apt-get update && apt-get install -y --no-install-recommends \
      pandoc \
      texlive-xetex \
      texlive-fonts-recommended \
      texlive-plain-generic \
      libglpk40 \
      jq \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------------
#  R package installation via renv — reproducible snapshot
# -----------------------------------------------------------------------------
# Copy renv artefacts first to maximise Docker layer caching
COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

# Install renv itself (rarely changes; gets its own cacheable layer)
RUN R -e "install.packages('renv', repos = 'https://cloud.r-project.org')"

# Restore the pinned environment. Uses a BuildKit cache mount for renv's
# package cache so failed/partial builds retain already-compiled packages —
# turns a 15-min rebuild into 1-2 min when only a few packages change.
# After restore, verify and force-install anything still missing (catches cases
# where restore silently exits early due to partial cache state).
RUN --mount=type=cache,target=/root/.cache/R/renv,sharing=locked \
    R -e "options(renv.config.pak.enabled = FALSE); \
          renv::restore(prompt = FALSE); \
          st <- renv::status(); \
          cat('--- renv::status() ---\n'); print(st); \
          missing <- c(st\$missing\$Package, st\$synchronized\$Package[!st\$synchronized\$Package %in% rownames(installed.packages())]); \
          missing <- unique(missing[!is.na(missing)]); \
          if (length(missing) > 0) { \
            cat('Force-installing missing packages:', paste(missing, collapse=', '), '\n'); \
            renv::install(missing, prompt = FALSE); \
          }"

# -----------------------------------------------------------------------------
#  Verify critical packages loaded
# -----------------------------------------------------------------------------
RUN R -e " \
    critical <- c('DESeq2', 'apeglm', 'SingleCellExperiment', 'biomaRt', \
                  'zellkonverter', 'qs', 'future.apply', 'rmarkdown'); \
    for (p in critical) { \
        stopifnot('package missing' = requireNamespace(p, quietly = TRUE)); \
        cat(sprintf('✓ %s %s\n', p, packageVersion(p))); \
    }"

# -----------------------------------------------------------------------------
#  Default working directory matches POWER_ANALYSIS_ROOT convention
# -----------------------------------------------------------------------------
WORKDIR /workspace
ENV POWER_ANALYSIS_ROOT=/workspace

# -----------------------------------------------------------------------------
#  Sanity checks at container startup
# -----------------------------------------------------------------------------
CMD ["R"]

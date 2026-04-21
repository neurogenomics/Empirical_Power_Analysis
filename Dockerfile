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
RUN apt-get update && apt-get install -y --no-install-recommends \
      pandoc \
      pandoc-citeproc \
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

# Install renv itself, then restore the pinned environment
RUN R -e "install.packages('renv', repos = 'https://cloud.r-project.org')" \
 && R -e "renv::restore(prompt = FALSE)"

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

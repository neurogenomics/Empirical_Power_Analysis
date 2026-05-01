# syntax=docker/dockerfile:1.7
# =============================================================================
#  Empirical Power Analysis — Reproducible Compute Environment
#  Base: Bioconductor 3.21 (R 4.5.1) on Ubuntu 24.04
#
#  Build:        DOCKER_BUILDKIT=1 docker build -t empirical-power-analysis:overnight .
#  Run:          ./docker-overnight.sh   (orchestrator handles everything)
#
#  Architecture:
#    * R packages installed into /opt/renv-library (NOT inside /workspace) so
#      a runtime bind-mount of the host repo on /workspace cannot shadow them.
#      Activated via RENV_PATHS_LIBRARY (which renv reads natively) and a small
#      shim in /etc/Rprofile.site.
#    * Python+anndata pre-installed via reticulate at build time so the first
#      sharding step does not have to download Miniconda over the network.
#    * Working directory is /workspace at runtime (host repo bind-mount target).
#      At build time we use /tmp/build to avoid littering /workspace with
#      installer artefacts.
#    * Runtime user is the host UID:GID via `docker run -u UID:GID`. /etc/passwd
#      is fixed up by the entrypoint so commands like whoami / Sys.info() work.
# =============================================================================

FROM bioconductor/bioconductor_docker:RELEASE_3_21

LABEL org.opencontainers.image.title="Empirical Power Analysis"
LABEL org.opencontainers.image.description="Reproducible environment for empirical scRNA-seq power analysis"
LABEL org.opencontainers.image.source="https://github.com/neurogenomics/Empirical_Power_Analysis"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.authors="Skene Lab, UK DRI Imperial"

# -----------------------------------------------------------------------------
#  System dependencies beyond bioconductor_docker defaults.
#
#  - pandoc / texlive: Rmd → PDF rendering for the paper.
#  - ccache:           cache C/C++ compiles between rebuilds.
#  - libglpk40:        igraph dependency surfaced by some Bioconductor pkgs.
#  - jq:               manifest assembly.
#  - tini:             PID 1 signal handler so SIGTERM reaches R cleanly.
#  - gosu:             drop privileges in entrypoint (used to add passwd entry).
#  - python3-venv etc: support reticulate's miniconda-free env path.
#
#  The cache mounts (BuildKit) keep apt downloads between rebuilds. Removing
#  /etc/apt/apt.conf.d/docker-clean is required because the bioconductor base
#  image installs that file, which wipes the cache after every install and
#  defeats BuildKit cache mounts.
# -----------------------------------------------------------------------------
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
      jq \
      tini \
      gosu \
      python3 \
      python3-pip \
      python3-venv \
      ca-certificates \
      curl

# -----------------------------------------------------------------------------
#  Pin the R library path OUTSIDE /workspace.
#
#  Why: at runtime we bind-mount the host repo onto /workspace. If the renv
#  library lived under /workspace/renv/library it would be shadowed by whatever
#  the host happens to have there (often nothing, or a different platform's
#  library), making the carefully-pinned packages installed at build time
#  invisible.
#
#  RENV_PATHS_LIBRARY is honoured by renv directly (no project-level config
#  needed), and we also expose it as the default site-wide library via
#  /etc/Rprofile.site so even non-renv R sessions see the packages.
# -----------------------------------------------------------------------------
ENV RENV_PATHS_LIBRARY=/opt/renv-library
ENV RENV_PATHS_CACHE=/opt/renv-cache
RUN mkdir -p /opt/renv-library /opt/renv-cache /opt/build \
 && chmod 0777 /opt/renv-library /opt/renv-cache

# Site-wide R profile that puts /opt/renv-library on .libPaths() before any
# project-level .Rprofile runs. This is the safety net when /workspace/.Rprofile
# is bind-mounted from the host: even if the host's .Rprofile is stale, the
# image's library is reachable.
RUN cat > /etc/Rprofile.site <<'EOF'
# Empirical Power Analysis — site profile
# Prepend the image's renv library so packages installed at build time are
# always reachable, even when /workspace is bind-mounted from a host that
# doesn't have its own renv library.
local({
  lib <- Sys.getenv("RENV_PATHS_LIBRARY", unset = "/opt/renv-library")
  if (dir.exists(lib)) {
    current <- .libPaths()
    if (!(lib %in% current)) {
      .libPaths(c(lib, current))
    }
  }
})
EOF

# -----------------------------------------------------------------------------
#  R package installation via renv (lockfile-pinned, retry-heavy).
# -----------------------------------------------------------------------------
WORKDIR /opt/build

COPY renv.lock           renv.lock
COPY renv/activate.R     renv/activate.R
COPY renv/settings.json  renv/settings.json
COPY docker/install_packages.R install_packages.R

# Install renv itself first (rare-changing layer; small download).
RUN --mount=type=cache,target=/root/.cache/R,sharing=locked \
    R -e "install.packages('renv', repos = 'https://cloud.r-project.org')"

# Restore the pinned environment with a hardened policy. The script lives in
# its own file so we can iterate on restore logic without busting the apt or
# renv-install layers above.
RUN --mount=type=cache,target=/root/.cache/R/renv,sharing=locked \
    --mount=type=cache,target=/root/.cache/R,sharing=locked \
    --mount=type=cache,target=/tmp/ccache,sharing=locked \
    --mount=type=cache,target=/opt/renv-cache,sharing=locked \
    mkdir -p ~/.R \
 && echo "CC=ccache gcc"     >> ~/.R/Makevars \
 && echo "CXX=ccache g++"    >> ~/.R/Makevars \
 && echo "CXX11=ccache g++"  >> ~/.R/Makevars \
 && echo "CXX14=ccache g++"  >> ~/.R/Makevars \
 && echo "CXX17=ccache g++"  >> ~/.R/Makevars \
 && export CCACHE_DIR=/tmp/ccache \
 && export RENV_CONFIG_CACHE_SYMLINKS=FALSE \
 && Rscript install_packages.R

# glmGamPoi is a Suggests of DESeq2 in renv.lock so renv::restore() skips it,
# but cartesian_downsampling/config.R uses fitType="glmGamPoi" and will hard-fail
# without it. Install via BiocManager so it lands in the same library tree.
RUN --mount=type=cache,target=/root/.cache/R/renv,sharing=locked \
    --mount=type=cache,target=/root/.cache/R,sharing=locked \
    R -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install("glmGamPoi", update=FALSE, ask=FALSE); stopifnot(requireNamespace("glmGamPoi", quietly=TRUE))'

# -----------------------------------------------------------------------------
#  Pre-bake the Python environment for reticulate (anndata).
#
#  Why: Run_Shard_Roussos.R calls reticulate::py_install("anndata") if anndata
#  isn't available. At runtime, with HOME=/tmp and an ephemeral container, that
#  install would have to download Miniconda + create the env every time the
#  container starts. We do it once at build time so the environment is in the
#  image and is found via RETICULATE_PYTHON.
#
#  We use a plain venv rather than miniconda — anndata has no native deps that
#  need conda, and pip wheels for Ubuntu 24.04 cover everything we need.
# -----------------------------------------------------------------------------
ENV RETICULATE_PYTHON=/opt/python-env/bin/python
ENV PYTHONDONTWRITEBYTECODE=1
RUN python3 -m venv /opt/python-env \
 && /opt/python-env/bin/pip install --no-cache-dir --upgrade pip wheel \
 && /opt/python-env/bin/pip install --no-cache-dir \
        anndata \
        h5py \
        numpy \
        pandas \
        scipy \
 && chmod -R a+rX /opt/python-env

# Sanity-check reticulate can drive the venv.
RUN R -e 'reticulate::use_python("/opt/python-env/bin/python", required = TRUE); reticulate::py_run_string("import anndata; print(\"anndata\", anndata.__version__)")'

# -----------------------------------------------------------------------------
#  Permissions: every path the runtime user might touch must be world-readable
#  / world-executable for directories. The runtime UID:GID is chosen by the
#  caller (docker run -u UID:GID) and won't have a /etc/passwd entry — so we
#  cannot chown to a known user. Using broad mode bits is the simplest path.
# -----------------------------------------------------------------------------
RUN chmod -R a+rX /opt/renv-library /opt/python-env /home/rstudio /usr/local/lib/R 2>/dev/null || true

# -----------------------------------------------------------------------------
#  Final critical-package check. Uses library() rather than requireNamespace()
#  so .onLoad hooks run — catches issues that requireNamespace() would miss
#  (e.g. shared library load failures, basilisk init).
# -----------------------------------------------------------------------------
RUN R -e ' \
    pkgs <- c("DESeq2","apeglm","glmGamPoi","SingleCellExperiment","biomaRt", \
              "zellkonverter","reticulate","basilisk","qs", \
              "Rcpp","RcppArmadillo","RcppEigen", \
              "future.apply","rmarkdown","ggplot2","dplyr","Matrix"); \
    bad <- character(); \
    for (p in pkgs) { \
        ok <- tryCatch({ suppressMessages(library(p, character.only = TRUE)); TRUE }, \
                       error = function(e) { cat(sprintf("FAIL %s: %s\n", p, conditionMessage(e))); FALSE }); \
        if (ok) cat(sprintf("OK   %s %s\n", p, packageVersion(p))) else bad <- c(bad, p); \
    }; \
    if (length(bad) > 0) stop(sprintf("Critical packages failed to load: %s", paste(bad, collapse=", ")))'

# -----------------------------------------------------------------------------
#  Runtime layout.
#
#  /workspace is the bind-mount target for the host repo at runtime.
#  POWER_ANALYSIS_ROOT defaults here; the wrapper scripts respect that env var.
#  Build-time scratch under /opt/build is not needed at runtime.
# -----------------------------------------------------------------------------
WORKDIR /workspace
ENV POWER_ANALYSIS_ROOT=/workspace
ENV CBERT_HEARTBEAT_FILE=/tmp/empirical-power-analysis.heartbeat

# Copy the entrypoint and healthcheck into the image (NOT bind-mounted) so
# they always reflect the image build and cannot be tampered with by host
# editors mid-run.
COPY docker/entrypoint.sh   /usr/local/bin/empirical-entrypoint
COPY docker/healthcheck.sh  /usr/local/bin/empirical-healthcheck
RUN chmod 0755 /usr/local/bin/empirical-entrypoint /usr/local/bin/empirical-healthcheck

# HEALTHCHECK reports unhealthy if the heartbeat file isn't touched within
# HEARTBEAT_STALE_MIN (default 30) minutes. The wrapper scripts inside the
# container update the heartbeat on every stage start and every per-iteration
# log line. A wedged container goes unhealthy within 30 min.
HEALTHCHECK --interval=5m --timeout=15s --start-period=10m --retries=2 \
    CMD /usr/local/bin/empirical-healthcheck

ENTRYPOINT ["/usr/bin/tini", "--", "/usr/local/bin/empirical-entrypoint"]
CMD []

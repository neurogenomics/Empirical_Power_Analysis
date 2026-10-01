# =============================================================================
#  Empirical Power Analysis — Reproducible Compute Environment
#  Base: Bioconductor 3.22 (R 4.5) on Ubuntu 24.04. The digest pins the exact
#  base image; the tag is kept only as human-readable provenance.
#
#  Build:        ./docker/docker-build.sh         (auto-detects BuildKit vs legacy)
#  Verify:       ./docker/docker-smoke.sh
#
#  This Dockerfile is intentionally BuildKit-OPTIONAL: it uses no `# syntax=`
#  directive, no `--mount=type=cache` mounts, and no RUN here-docs, so it builds
#  identically under the legacy builder and under BuildKit. That keeps it
#  "runnable anywhere", including hosts without the buildx plugin installed.
#
#  Architecture:
#    * R packages installed into /opt/renv-library (NOT inside /workspace) so
#      a runtime bind-mount of the host repo on /workspace cannot shadow them.
#      Made visible to every R session via R_LIBS_SITE, baked into
#      R_HOME/etc/Renviron.site by docker/install_packages.R (renv installs into
#      a platform-specific leaf, so we record the resolved leaf, not the root).
#    * Python+anndata pre-installed via reticulate at build time so the first
#      sharding step does not have to download Miniconda over the network.
#    * Working directory is /workspace at runtime (host repo bind-mount target).
#      At build time we use /tmp/build to avoid littering /workspace with
#      installer artefacts.
#    * Runtime user is the host UID:GID via `docker run -u UID:GID`. /etc/passwd
#      is fixed up by the entrypoint so commands like whoami / Sys.info() work.
# =============================================================================

FROM bioconductor/bioconductor_docker:RELEASE_3_22@sha256:e34177896c1579f66d48b756e4229b58237499183a089da084af72d7aba1d40a

LABEL org.opencontainers.image.title="Empirical Power Analysis"
LABEL org.opencontainers.image.description="Reproducible environment for empirical scRNA-seq power analysis"
LABEL org.opencontainers.image.source="https://github.com/neurogenomics/Empirical_Power_Analysis"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.authors="Skene Lab, UK DRI Imperial"

# -----------------------------------------------------------------------------
#  System dependencies beyond the digest-pinned bioconductor_docker base.
#
#  - pandoc:           already present in the pinned base image.
#  - texlive:          Rmd → PDF rendering for the paper.
#  - ccache:           cache C/C++ compiles between rebuilds.
#  - libglpk40:        already present in the pinned base image.
#  - jq:               manifest assembly.
#  - tini:             PID 1 signal handler so SIGTERM reaches R cleanly.
#  - gosu:             drop privileges in entrypoint (used to add passwd entry).
#  - python3/pip/venv: already present in the pinned base image.
#
# -----------------------------------------------------------------------------
COPY docker/apt-packages.txt /tmp/apt-packages.txt
RUN apt-get update \
 && grep -vE '^[[:space:]]*(#|$)' /tmp/apt-packages.txt \
      | xargs apt-get install -y --no-install-recommends \
 && rm -f /tmp/apt-packages.txt \
 && rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------------
#  Pin the R library path OUTSIDE /workspace.
#
#  Why: at runtime we bind-mount the host repo onto /workspace. If the renv
#  library lived under /workspace/renv/library it would be shadowed by whatever
#  the host happens to have there (often nothing, or a different platform's
#  library), making the carefully-pinned packages installed at build time
#  invisible.
#
#  renv installs into a platform/version-specific *leaf* under
#  RENV_PATHS_LIBRARY (e.g. /opt/renv-library/linux-ubuntu-noble/R-4.5/...).
#  To make those packages visible to every plain R/Rscript session — the
#  build-time sanity checks, the smoke test, and the pipeline at runtime —
#  install_packages.R bakes the resolved leaf into R's site Renviron
#  (R_HOME/etc/Renviron.site) as R_LIBS_SITE. That is the native, no-profile
#  mechanism R reads at startup; it does not depend on renv re-activating
#  (which we disable at runtime in docker/entrypoint.sh).
# -----------------------------------------------------------------------------
ENV RENV_PATHS_LIBRARY=/opt/renv-library
ENV RENV_PATHS_CACHE=/opt/renv-cache
RUN mkdir -p /opt/renv-library /opt/renv-cache /opt/build \
 && chmod 0777 /opt/renv-library /opt/renv-cache

# -----------------------------------------------------------------------------
#  R package installation via renv (lockfile-pinned, retry-heavy).
# -----------------------------------------------------------------------------
WORKDIR /opt/build

COPY renv.lock           renv.lock
COPY renv/activate.R     renv/activate.R
COPY renv/settings.json  renv/settings.json
COPY docker/install_packages.R install_packages.R
COPY docker/python-requirements.txt python-requirements.txt

# Install the exact renv bootstrap used by renv/activate.R and renv.lock.
ENV RENV_VERSION=1.2.2
ENV RENV_TARBALL_SHA256=fd41b9ac1ebbbaf4c2a01bc16f8676f1236db0ee94f98e97373fbc00296be350
RUN R -e ' \
    url <- sprintf("https://cran.r-project.org/src/contrib/Archive/renv/renv_%s.tar.gz", Sys.getenv("RENV_VERSION")); \
    dest <- tempfile(fileext = ".tar.gz"); \
    download.file(url, dest, mode = "wb", quiet = FALSE); \
    actual <- unname(tools::sha256sum(dest)); \
    expected <- Sys.getenv("RENV_TARBALL_SHA256"); \
    if (!identical(actual, expected)) stop(sprintf("renv tarball SHA256 mismatch: expected %s, got %s", expected, actual)); \
    install.packages(dest, repos = NULL, type = "source"); \
    stopifnot(identical(as.character(packageVersion("renv")), Sys.getenv("RENV_VERSION"))) \
'

# Restore the pinned environment with a hardened policy. The script lives in
# its own file so we can iterate on restore logic without busting the apt or
# renv-install layers above. ccache still accelerates compiles within a single
# build; cross-build caching (the former BuildKit cache mounts) is dropped to
# keep the build portable to hosts without buildx.
RUN mkdir -p ~/.R \
 && echo "CC=ccache gcc"     >> ~/.R/Makevars \
 && echo "CXX=ccache g++"    >> ~/.R/Makevars \
 && echo "CXX11=ccache g++"  >> ~/.R/Makevars \
 && echo "CXX14=ccache g++"  >> ~/.R/Makevars \
 && echo "CXX17=ccache g++"  >> ~/.R/Makevars \
 && export CCACHE_DIR=/tmp/ccache \
 && export RENV_CONFIG_CACHE_SYMLINKS=FALSE \
 && Rscript install_packages.R

# -----------------------------------------------------------------------------
#  Pre-bake the Python environment for reticulate (anndata).
#
#  Why: Run_Shard_Roussos_Full.R requires anndata and now fails if it is not
#  available. Runtime Python installs are forbidden for paper reproducibility,
#  so the environment is baked once at build time and found via RETICULATE_PYTHON.
#
#  We use a plain venv rather than miniconda — anndata has no native deps that
#  need conda, and pip wheels for Ubuntu 24.04 cover everything we need.
# -----------------------------------------------------------------------------
ENV RETICULATE_PYTHON=/opt/python-env/bin/python
ENV PYTHONDONTWRITEBYTECODE=1
RUN python3 -m venv /opt/python-env \
 && /opt/python-env/bin/pip install --no-cache-dir --only-binary=:all: --require-hashes \
        -r /opt/build/python-requirements.txt \
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
    pkgs <- c("DESeq2","apeglm","glmGamPoi","edgeR","limma","SingleCellExperiment","biomaRt", \
              "zellkonverter","reticulate","basilisk","qs", \
              "Rcpp","RcppArmadillo","RcppEigen", \
              "future.apply","rmarkdown","bookdown","ggplot2","dplyr","Matrix", \
              "xgboost","writexl","patchwork","png"); \
    bad <- character(); \
    for (p in pkgs) { \
        ok <- tryCatch({ suppressMessages(library(p, character.only = TRUE)); TRUE }, \
                       error = function(e) { cat(sprintf("FAIL %s: %s\n", p, conditionMessage(e))); FALSE }); \
        if (ok) cat(sprintf("OK   %s %s\n", p, packageVersion(p))) else bad <- c(bad, p); \
    }; \
    if (length(bad) > 0) stop(sprintf("Critical packages failed to load: %s", paste(bad, collapse=", ")))'

# The figure composers (compose_figures.py, compose_supplementary.py) run on the
# venv Python and need Pillow, matplotlib and the DejaVu TrueType fonts; fail the
# build here rather than silently producing bitmap-font panels at runtime.
RUN /opt/python-env/bin/python -c "import PIL, matplotlib, numpy, pandas; print(\"PIL\", PIL.__version__, \"matplotlib\", matplotlib.__version__)" \
 && test -f /usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf

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

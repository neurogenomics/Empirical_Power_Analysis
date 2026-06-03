# Docker Build Handover

**For:** Toby Clark
**From:** Chris Crampton
**Date:** 2026-06-03
**Repo:** `neurogenomics/Empirical_Power_Analysis` (branch: `main`)

## What this is

A Docker container that reproduces the entire Empirical Power Analysis pipeline (the scRNA-seq power analysis paper). The container packages R 4.5 + Bioconductor 3.21 + all dependencies so the full pipeline can run unattended without worrying about library versions.

## Goal

Get the Docker image to **build successfully** and pass the **smoke test**. If that works, try running the power curve fitting pipeline inside it.

## Prerequisites

- Docker with BuildKit enabled (default on Docker 20.10+)
- ~15 GB disk for the image build
- ~50 GB free on the filesystem where `Datalink/` lives (for pipeline outputs)
- The repo cloned locally with `Datalink/` populated (raw data, shards, etc.)

## Quick start

```bash
cd /path/to/Empirical_Power_Analysis

# Step 1: Build the image + run smoke test (no pipeline execution)
./docker-overnight.sh --build-only

# Step 2: If smoke passes, try the power curve fitting pipeline
./docker-overnight.sh --pipeline pcf --max-cores 24

# Step 3: If feeling brave, full DAG
nohup ./docker-overnight.sh --pipeline all --max-cores 24 \
    > nohup.docker-overnight.out 2>&1 & disown
```

## What `--build-only` does

1. **Pre-flight checks**: Docker available? BuildKit? Enough disk? Enough cores?
2. **Builds the image** (`docker build`): installs system deps, R packages from `renv.lock`, Python env (anndata), verifies critical packages load
3. **Runs a smoke test**: starts the container with `PIPELINE=smoke`, which loads every critical R package and checks reticulate/anndata work

If the smoke test passes, the environment is good.

## Key files

| File | What it does |
|------|-------------|
| `docker-overnight.sh` | Orchestrator script (build, smoke, run) |
| `Dockerfile` | Image definition (Bioconductor 3.21 base, system library install) |
| `docker/entrypoint.sh` | Container entrypoint (privilege drop, pipeline dispatch) |
| `docker/install_packages.R` | Parses `renv.lock`, installs missing packages into system library |
| `docker/healthcheck.sh` | Container health check (heartbeat file) |
| `run_all.sh` | Full pipeline driver (used inside container for `--pipeline all`) |
| `run_power_curve_fitting.sh` | PCF-only driver (used inside container for `--pipeline pcf`) |

## How the Docker image works

- **Base image**: `bioconductor/bioconductor_docker:RELEASE_3_21` (Ubuntu 24.04, R 4.5)
- **R packages**: installed into the **system library** (`/usr/local/lib/R/site-library/`), NOT via renv. `renv.lock` is used as a manifest only. This avoids path conflicts between renv's platform-suffixed paths and Docker bind mounts.
- **Python**: a plain venv at `/opt/python-env/` with anndata, h5py, numpy, pandas, scipy. Reticulate is configured to use it via `RETICULATE_PYTHON`.
- **Runtime**: the repo is bind-mounted at `/workspace`. The container runs as your host UID (via `gosu` privilege drop) so output files have correct ownership.
- **renv disabled at runtime**: `EPA_DOCKER=1` env var tells `.Rprofile` to skip renv activation.

## Useful flags

```
--build-only        Build + smoke test only, don't run a pipeline
--skip-build        Reuse existing image (skip rebuild)
--pipeline pcf      Just power curve fitting (default)
--pipeline all      Full DAG end-to-end
--max-cores N       Cores for R parallelism (default: min(24, nproc))
--no-overwrite      Skip steps whose outputs already exist (default)
--force             Regenerate everything
--dry-run           Print what would happen, run nothing
```

## Where logs go

```
logs/docker-overnight/<UTC-timestamp>/
    00_orchestrator.log    Timestamped orchestrator events
    01_build.log           Docker build output
    02_smoke.log           Smoke test output
    03_run.log             Pipeline execution log
    manifest.json          Provenance (image digest, git SHA, args, timings)
```

## Likely failure points

1. **R package compilation**: if a system library is missing, the build will fail during `install_packages.R`. Check `01_build.log` for the specific package and missing `-dev` dependency.
2. **Network issues during build**: the script downloads from CRAN, Bioconductor, and P3M mirrors. If behind a proxy, you may need to configure Docker's network.
3. **renv.lock drift**: if `renv.lock` references a package version that's been removed from CRAN/Bioc, the build will fail. The install script has retry logic and a fallback to plain `install.packages`, but truly removed packages need manual intervention.
4. **Disk space**: the Bioconductor base image is large (~4 GB), and the full build with all R packages can reach ~12-15 GB.
5. **Datalink not mounted**: if you run a pipeline without `Datalink/` present, the R scripts will fail looking for input data. `docker-overnight.sh` handles the bind mount automatically.

## Exit codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Pipeline reported partial failures |
| 2 | Pre-flight check failed |
| 3 | Build failed |
| 4 | Smoke test failed |
| 5 | Container couldn't start |
| 6 | Container crashed |

## If it breaks

- Check the log files in `logs/docker-overnight/<timestamp>/`
- The build log (`01_build.log`) is the most useful for image build failures
- The smoke log (`02_smoke.log`) shows which packages failed to load
- If a specific R package won't install, try adding the missing system dependency to the `apt-get install` block in the `Dockerfile` and rebuilding

## Questions?

Ask Chris (or the CBERT research assistant).

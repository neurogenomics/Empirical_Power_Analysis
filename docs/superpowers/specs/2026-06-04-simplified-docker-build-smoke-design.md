# Simplified Docker build + smoke — design

**Branch:** `docker-tc`
**Date:** 2026-06-04
**Author:** Toby Clark (with Claude)
**Status:** approved

## Problem

The repo owner has not been able to get Docker working. The current entry point,
`docker-overnight.sh`, is a 524-line orchestrator that bundles build, smoke test,
pipeline run, a hand-rolled `manifest.json` writer, a heartbeat/`HEALTHCHECK`
mechanism, four-layer signal forwarding, disk preflight, core clamping, and an
exit-code translation table. The essential job — *build the image and confirm the
environment is correct* — is buried under orchestration scaffolding that is hard
to debug and is the most likely place for subtle breakage to hide.

Separately, the current `DOCKER_HANDOVER.md` documents an **older, different**
image design (system-library install, an `EPA_DOCKER=1` runtime guard that exists
nowhere in the repo). The implemented design is renv.lock-pinned into
`/opt/renv-library` with an additive `.Rprofile` / `Rprofile.site` shim — which is
sound. The docs have drifted from the code.

## Goal / success bar

A host with Docker can:

1. **Build** the image with one short command.
2. **Smoke-test** it (all critical R packages + reticulate/anndata load) with one
   short command, and see `SMOKE TEST PASSED`.

This proves the environment is correct. It does **not** require raw data to be
present and does **not** run the multi-day pipeline. Verifiable in a WSL2 session
once the Docker daemon is reachable.

## Non-goals

- Running `pcf` / `all` pipelines to completion (existing `docker-overnight.sh`
  remains available for that; we do not delete or modify it).
- Changing package versions or the reproducibility model (DESeq2 version
  sensitivity is load-bearing per `CLAUDE.md`).
- Manifest provenance, heartbeat, `HEALTHCHECK`, auto-restart, signal forwarding.

## Design

### 1. Two new minimal scripts

`docker-build.sh` (~40 lines)
- Preflight: `docker` on PATH? daemon reachable (`docker info`)? `buildx`
  available? Each failure prints a specific, actionable message and exits non-zero.
- Build: `DOCKER_BUILDKIT=1 docker build -t <image> "<repo>"`. BuildKit is **forced
  on** rather than removed from the Dockerfile — BuildKit is default in Docker 23+
  and the Dockerfile's cache mounts / here-doc / `# syntax=` directive require it.
  This is the "runnable anywhere" fix: guarantee BuildKit, fail fast and clearly if
  it is genuinely unavailable.
- Flags: `--image NAME:TAG` (default `empirical-power-analysis:overnight`), `-h`.
- Foreground, plain stdout logging. No nohup/background.

`docker-smoke.sh` (~40 lines)
- Runs the image with `PIPELINE=smoke`, the **same runtime contract as a real run**:
  bind-mount repo on `/workspace`, `Datalink` on `/workspace/Datalink`,
  `TARGET_UID`/`TARGET_GID` set to the host user (so the gosu privilege drop and
  file ownership behave as in production), `HOME=/tmp`,
  `POWER_ANALYSIS_ROOT=/workspace`.
- Greps output for `SMOKE TEST PASSED`; exits non-zero (with a log tail) otherwise.
- Flags: `--image NAME:TAG`, `-h`.

### 2. Dockerfile — unchanged (for now)

Keep the renv.lock-pinned build, BuildKit cache mounts, ccache, and the pre-baked
Python venv. These provide reproducibility and network-flake hardening. The only
real risk to "build + smoke works" is whatever specifically breaks the build for
the owner — which we discover by **running it here once the daemon is up** and
fixing that concrete failure, rather than changing the image speculatively.

### 3. Entrypoint — unchanged (for now)

`docker/entrypoint.sh` already dispatches `smoke` correctly and drops privileges via
gosu. The heartbeat loop is harmless during a fast smoke run. Removing
heartbeat/healthcheck wiring would only benefit long runs we are not targeting and
risks breaking the working `pcf`/`all` paths. Leave it untouched unless the smoke
run reveals a concrete bug.

## Verification

1. Bring up a reachable Docker daemon in the test environment.
2. `./docker-build.sh` — capture and fix whatever the build actually fails on.
3. `./docker-smoke.sh` — confirm `SMOKE TEST PASSED`.
4. Iterate on the Dockerfile/entrypoint **only** in response to observed failures.

## Open question to resolve during implementation

- Where does the build currently fail for the owner? Unknown until we run it.
  Likely candidates (per `DOCKER_HANDOVER.md`): missing system `-dev` dep during R
  package compilation, network flakiness to CRAN/Bioc/P3M, or renv.lock drift. The
  hardened `install_packages.R` already mitigates the latter two.

## Follow-up (out of scope here, worth noting)

`DOCKER_HANDOVER.md` documents a stale image design. Once build+smoke works, it
should be corrected to match the implemented renv-library design.

# Changelog

All notable changes to this project will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Docker image based on `bioconductor/bioconductor_docker:RELEASE_3_22`
  (matches the Bioconductor 3.22 versions pinned in `renv.lock`)
- `docker-build.sh` / `docker-smoke.sh` / `docker-run.sh` — single-purpose
  build / verify / run scripts, plus `docker-all.sh` chaining them with
  per-phase `--skip-*` flags. See `DOCKER.md`.
- `docker/entrypoint.sh` — container entrypoint dispatching to the chosen
  wrapper (`run_power_curve_fitting.sh` by default)

### Changed
- Replaced the single 524-line `docker-overnight.sh` orchestrator with the four
  composable scripts above; made the Dockerfile BuildKit-optional; fixed the
  renv library visibility (`R_LIBS_SITE` baked into `Renviron.site`) so the
  pinned packages load in every R session; aligned the base image and Bioc repos
  to 3.22 to match `renv.lock`.
- `run_all.sh` — host-native bash driver running the strict-serial DAG with
  `--max-cores`, `--no-overwrite`/`--force`, `--mode`, `--only`, `--skip`,
  `--from`, `--sharding-datasets`, `--stop-on-fail`, `--dry-run`
- `run_power_curve_fitting.sh` — thin wrapper that runs only the DAG
  required for power-curve fitting (Roussos shards only)
- `run_all.config` — editable defaults for the bash driver
- Explicit `BiocManager::install("glmGamPoi")` in the Dockerfile (only listed
  as a `Suggests` of DESeq2 in `renv.lock`, so `renv::restore()` does not
  install it; cartesian downsampling needs it for `fitType = "glmGamPoi"`)
- GitHub Actions Docker image build + publish to GHCR on tag release
- `AGENTS.md` — canonical AI-agent entry point
- `CITATION.cff` — machine-readable citation metadata
- `docs/architecture.md` — Architecture Decision Records
- `docs/reproducing-the-paper.md` — reproduction walkthrough
- `docs/model-card-predict_power.md` — model card for `predict_power.R`
- `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, PR + issue templates

### Changed
- Paper source (`paper/paper.Rmd`) brought into the repository
- Issue / PR templates rewritten to drop Nextflow-specific fields

### Removed
- Nextflow DSL2 pipeline (`main.nf`, `nextflow.config`, `iclneuro.config`,
  `conf/`, `modules/local/`) — superseded by `run_all.sh` + Docker image
  (see ADR-001b in `docs/architecture.md`)
- `.github/workflows/ci.yml` and `.github/workflows/zenodo.yml` (Nextflow-
  specific; `.github/workflows/docker.yml` retained for image builds)

### Fixed

### Security

---

## [1.0.0] — Pending (target 2026-05-01)

Initial bioRxiv-accompanying release.

### Added
- Full empirical power analysis pipeline across PsychAD, SEA-AD, ROSMAP, Rexach,
  Johansen cohorts; GTEx bulk concordance
- `predict_power.R` standalone planning function with fitted coefficients
- `power_lookup_table.csv` pre-computed grid

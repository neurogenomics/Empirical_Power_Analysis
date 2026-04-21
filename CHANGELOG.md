# Changelog

All notable changes to this project will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Nextflow DSL2 pipeline (`main.nf`, `modules/local/**`, `conf/*`) wrapping the
  existing R-script pipelines with SLURM executor support
- Docker image based on `bioconductor/bioconductor_docker:RELEASE_3_21`
- GitHub Actions CI: R linting, Nextflow stub-run, unit tests, integration toy
  pipeline, determinism check
- GitHub Actions Docker image build + publish to GHCR on tag release
- `AGENTS.md` — canonical AI-agent entry point
- `CITATION.cff` — machine-readable citation metadata
- `docs/architecture.md` — Architecture Decision Records
- `docs/reproducing-the-paper.md` — reproduction walkthrough
- `docs/model-card-predict_power.md` — model card for `predict_power.R`
- `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, PR + issue templates
- Pre-commit hooks for R styling, YAML/Nextflow syntax, leak detection

### Changed
- Paper source (`paper/paper.Rmd`) brought into the repository

### Deprecated

### Removed

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

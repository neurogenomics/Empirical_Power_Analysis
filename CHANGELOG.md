# Changelog

All notable changes to this project will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `pipeline/`: one-command reproduction on a PBS Pro cluster (written for
  Imperial CX3). Submits every stage as a dependency-chained job inside the
  pinned Apptainer image, checks raw inputs and the offline `reference/` bundle
  against checksums, and writes per-stage result tars with `MANIFEST.sha256` and
  `PROVENANCE.txt`. `pipeline/unpack_results.sh` unpacks them into `Datalink/`
  outside the pipeline.
- Three DE methods (DESeq2, edgeR, limma-voom) run on identical samples across
  the downsampling grid, with 10 robustness runs and cross-method truth scoring.
- XGBoost power model (`R/power_curve_fitting/02_Fit_XGBoost.R`,
  `03_Select_XGBoost.R`, `04_Export_XGBoost.R`), trained on the realistic
  regime (>= 100 cells per donor), with a portable bundle and `predict_power()`
  wrapper in `Datalink/Power_Curve_Fitting/Exports/xgboost/realistic/`.
- Gene detectability classification for Figure 6 (`R/gene_classification/`).
- `R/supplementary_data/build_supplementary_data.R`: builds Supplementary
  Data 1-4 with a byte-reproducible zip.
- `_paper_values.R` (`val()`) and `_supp_fig_index.R` (`SF()`): the manuscript
  reads quoted numbers from source CSVs and supplementary figure numbers from one
  index at knit time.
- Docker image on a digest-pinned `bioconductor/bioconductor_docker:RELEASE_3_22`
  base, with apt packages pinned in `docker/apt-packages.txt`, R packages from
  `renv.lock` and hashed Python requirements; `docker/build-sif.sh` converts it
  to an Apptainer image.
- `docker/docker-build.sh` / `docker-smoke.sh` / `docker-run.sh`, chained by
  `docker-all.sh` with per-phase `--skip-*` flags. See `docs/docker.md`.
- `run_all.sh`, `run_power_curve_fitting.sh` and `run_all.config`: host-native
  drivers.
- `CITATION.cff`, `docs/architecture.md`, `docs/model-card-predict_power.md`,
  `CONTRIBUTING.md`.

### Changed
- Scope narrowed to PsychAD only. Cross-cohort validation (SEA-AD, ROSMAP,
  Johansen, Rexach) and the GTEx bulk comparison were removed from the paper.
- Log fold changes for all three DE methods are now shrunk with ashr (previously
  apeglm for DESeq2 and no shrinkage for edgeR and limma-voom).
- Ground-truth threshold set to |LFC| > 0.2 (BH-adjusted p < 0.05).
- Manuscript source is `nature_comms_draft.rmd` / `nature_comms_supplementary.rmd`.
- `glmGamPoi` and its dependencies are pinned in `renv.lock` rather than
  installed separately in the Dockerfile.

### Removed
- Logistic-regression power model and its exports (`predict_power.R` v4,
  `power_lookup_table.csv`), superseded by the XGBoost model.
- Nextflow DSL2 pipeline (`main.nf`, `nextflow.config`, `iclneuro.config`,
  `conf/`, `modules/local/`), superseded by `pipeline/` and the bash drivers
  (see ADR-001b and ADR-007 in `docs/architecture.md`).

## [1.0.0] (pending)

First release, to accompany the paper.

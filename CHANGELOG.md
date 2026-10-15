# Changelog

All notable changes to this project will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Stage `paper_values` (`R/figure_data_generation/09_Paper_Values.R`) writes
  `Paper_Source_Data/paper_values.csv`: every data-derived number quoted in
  the text (cohort sizes, grid recall/precision claims, DEG loss, power-model
  calibration, gene-group detectability), each with a written definition and
  its source file. `val()` falls back to this table for keys `_paper_values.R`
  does not register. Step 1 also writes `cohort_prefilter.csv` (donors before
  filtering), and the Figure 6 script writes `oligo_gene_group_power_summary.csv`
  (the percentages in its column headers).
- `docker/smoke_test.R`: checks every `renv.lock` and Python pin at its exact
  version; runs at the end of the image build and inside every `.sif`.

### Fixed
- The container build could install the Rcpp family at the latest CRAN version
  instead of the `renv.lock` pin (the 2026-10-01 image had Rcpp 1.1.2 and
  RcppParallel 6.2.1 against pins 1.1.1-1 and 5.1.11-2). The build now
  requests pinned versions and fails on any version mismatch.
- `renv.lock` pinned Rcpp 1.1.1-1, a GitHub development build that was never
  released to CRAN, so the lockfile could not be restored by anyone. Repinned
  to Rcpp 1.1.1-1.1, the nearest CRAN release (2026-04-24).

### Changed
- Documentation merged from the PsychAD-only refresh and aligned with the
  current tree. `pipeline/run.sh` is the only supported way to run the code,
  on a PBS cluster or with `--local` on one host. Container documentation
  moves from `docs/docker.md` to `docs/container.md`.
- One XGBoost power model, trained on the `n_cells >= 100` rows
  (`01e_Filter_Design.R`, formerly `01e_Filter_Realistic.R`; training table
  `Fitting_Data/aggregated_fitting_data.qs`). The full-grid "baseline" model
  and the "realistic" naming are gone; outputs are in `Models/xgboost/` and
  `Exports/xgboost/`.
- The hyperparameter search grows from 36 to 108 configurations
  (`max_depth` {3, 4, 6, 10} x `eta` {0.01, 0.03, 0.1} x `min_child_weight`
  {1, 10, 25} x `alpha` {0, 0.5, 1}), up to 20,000 rounds with early stopping,
  ranked by out-of-fold bin-balanced MAE with per-fold values in
  `model_selection_folds.csv`. Stages `xgb_search` (array of 108),
  `xgb_select`, `xgb_shap` and `xgb_post` replace the `xgb_base_*` and
  `xgb_real_*` stages, and `08b_Baseline_OOF.R` is removed. XGBoost stages run
  with 8 threads.
- Figure 5B shows hyperparameter bracketing with fold-paired standard errors;
  Figure 5C shows the model's out-of-fold per-class precision, recall and F1.
- Depth and cell-count equalisation (stage `deg_loss`) runs for
  oligodendrocytes as well as astrocytes; Supplementary Table 3 has a cell-type
  column.
- Container image tagged `empirical-power-analysis:v1.0.0`, built with
  `docker/docker-build.sh`, checked with `docker/docker-smoke.sh` and converted
  with `docker/build-sif.sh`; to be published as
  `ghcr.io/neurogenomics/empirical_power_analysis`.

### Added
- AD vs control contrast (`EPA_CONTRAST=ad_vs_control`, `R/shared/contrast.R`,
  `R/ad_contrast/`, `R/sharding/Build_AD_Control_Shards.R`): AD as the only
  recorded diagnosis vs no recorded diagnosis, sex as a covariate, a
  cells-mode grid with the same seeds and a full-cohort DE comparison, to show
  that recall on a disease contrast follows the same pattern as on sex. The
  power model is not tested on it. Stages `ad_shard`, `ad_gt_manifest`,
  `ad_grid`, `ad_cartesian` and `ad_de`; outputs in `Datalink/AD_vs_Control/`.

### Removed
Everything not on the reproduction path is untracked (commit `f14c924`; the
files remain in the git history):
- Superseded runners: `run_all.sh`, `run_all.config`,
  `run_power_curve_fitting.sh`, `render_pdf.R`, `docker/docker-run.sh`,
  `docker/docker-all.sh`, `hpc/`.
- Non-PsychAD code and outputs: `R/bulk_deg_processing/`, the other-cohort
  sharders, `Run_Harmonise_Gene_IDs.R`, `Audit_Roussos_PMI.R`, the
  cross-dataset diagnostics in `R/full_ground_truth/`,
  `Datalink/Bulk_Analysis/`, `Diagnostic_Results/` and `FDP_Summary.csv`.
- Logistic-model and off-paper analyses: `R/power_curve_fitting/run_pipeline.R`,
  `01d`, `02_Fit_Models.R`, `02b`, `03` to `06`, `04b`, `04c` and `10`;
  `R/figure_data_generation/` steps 2-5 and 7; `R/eda_roussos/` steps 3 and 6
  and two standalone scripts; unused `R/gene_classification/` and
  `R/figure_composition/` scripts; the matching `Paper_Source_Data/` CSVs and
  `Datalink/Dataset_Statistics_Table/`.
- One-off `_run_step*.R` helpers with hard-coded server paths.
- `Reproducibility/` and `R/shared/capture_environment.R` (a host snapshot;
  the container defines the environment).

## [Unreleased] (docs refresh, 2026-10-06)

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

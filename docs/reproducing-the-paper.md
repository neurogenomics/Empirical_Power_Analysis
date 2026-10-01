> **Superseded (2026-10-01).** This guide describes the earlier single-host
> route (DESeq2 only, logistic power model, five cohorts). It does not
> reproduce the current paper. For that, use the one-command cluster pipeline
> in [`pipeline/README.md`](../pipeline/README.md).

# Cold-Start Reproduction Guide

This document walks through reproducing every paper artefact in this
repository from scratch, as if you were a fresh researcher who has just
cloned the repo and obtained the raw data files.

It is deliberately **more operational than `README.md`** — concrete shell
commands, dependency graph, expected runtimes, logging conventions, and
how to diff the regenerated outputs against the shipped ones. For
scientific context and the figure lookup table, see `README.md`.

The guide is **modular**. Sections 0-3 are mandatory setup. Sections
A-I are the nine pipelines, each standalone with its own prereqs,
command, outputs, and verification. You can reproduce the whole paper
or cherry-pick a single pipeline.

---

## 0. Workspace

Pick a workspace directory (anywhere with a few hundred GB free). We
will refer to it as `$REPRO` throughout:

```bash
export REPRO=/path/to/your/workspace
mkdir -p "$REPRO"
cd "$REPRO"
git clone https://github.com/neurogenomics/Empirical_Power_Analysis.git
cd Empirical_Power_Analysis
export POWER_ANALYSIS_ROOT="$REPRO/Empirical_Power_Analysis"
mkdir -p logs
```

`POWER_ANALYSIS_ROOT` is the environment variable each pipeline's
`config.R` reads to locate the project root. Exporting it is optional
(the `here` package will auto-detect from `.here`) but recommended when
running inside `screen` / `tmux` or non-interactive jobs.

> The variable name is legacy (this repo was cherry-picked from a
> private `Power_Analysis_Paper` repo). The scripts still read it.

## 1. Provide the raw data

The repository ships **final paper artefacts only**. Raw inputs
(`.h5ad`, `.h5Seurat`, `.qs`, `.tar`) are not redistributed here — see
`Datalink/Dataset_Statistics_Table/Supplementary_Table_1_Dataset_Characteristics.csv`
for the original sources.

Place the raw files in `Datalink/raw_data/`. If you already have them
elsewhere on the same machine, a symlink is the cleanest option:

```bash
rm -rf "$REPRO/Empirical_Power_Analysis/Datalink/raw_data"
ln -s /path/to/your/raw_data \
      "$REPRO/Empirical_Power_Analysis/Datalink/raw_data"
```

Expected files (names the sharding scripts look for):

| Dataset | Files |
| --- | --- |
| PsychAD / Roussos | `Roussos_{Astro,Micro,Oligo}.qs` (+ `_metadata.txt`), `PsychAD-{HBCC,MSSM,RADC}.{qs,h5ad}` |
| SEA-AD / Allen    | `allen_{astrocytes,microglia,oligo}.qs` (+ `_metadata.txt`) |
| ROSMAP / DeJager  | `DeJager_{astrocytes,microglia,oligodendroglia}.h5Seurat`, `DeJagerMeta.csv` |
| Johansen          | `Johansen_{Astro,micro,oligo}.h5ad.tar` (extracted automatically) |
| Rexach            | `Rexach.h5ad` |

## 2. Install R dependencies

```bash
Rscript -e 'install.packages(c("here","qs","Matrix","dplyr","tidyr","stringr","ggplot2","patchwork","scales","RColorBrewer","viridis","UpSetR","ggnewscale","ggrepel","pheatmap","cowplot","shadowtext","akima","data.table","future","future.apply","reticulate","testthat","knitr","BiocManager"))'
Rscript -e 'BiocManager::install(c("DESeq2","apeglm","SingleCellExperiment","SummarizedExperiment","BiocParallel","biomaRt","AnnotationDbi","org.Hs.eg.db","rhdf5","zellkonverter"), update=FALSE, ask=FALSE)'
```

Check for anything missing before the long runs:

```bash
Rscript -e 'pkgs <- c("DESeq2","apeglm","SingleCellExperiment","SummarizedExperiment","BiocParallel","biomaRt","AnnotationDbi","org.Hs.eg.db","rhdf5","zellkonverter","here","qs","Matrix","dplyr","tidyr","stringr","ggplot2","patchwork","scales","RColorBrewer","viridis","UpSetR","ggnewscale","ggrepel","pheatmap","cowplot","shadowtext","akima","data.table","future","future.apply","reticulate","testthat","knitr"); miss <- pkgs[!pkgs %in% rownames(installed.packages())]; cat("Missing:", length(miss), "\n"); print(miss)'
```

## 2.5 One-shot driver: `run_all.sh`

If you just want to regenerate everything end-to-end and walk away for a few
days, use the bundled bash driver instead of running each pipeline by hand:

```bash
# From the repo root, after sections 0–2 are done:
./run_all.sh --max-cores 24 --no-overwrite

# Recommended for unattended multi-day runs:
nohup ./run_all.sh --max-cores 24 --no-overwrite > nohup.run_all.out 2>&1 &
disown
```

Defaults live in `run_all.config` (edit for permanent changes); CLI flags
override them per-invocation. Useful flags:

| Flag | Meaning |
| --- | --- |
| `--max-cores N` | cap every Rscript step at N cores (default 24; memory is the usual bottleneck) |
| `--no-overwrite` | skip work whose outputs already exist (per-pipeline `FORCE_RECALC=FALSE`) |
| `--force` | regenerate everything (per-pipeline `FORCE_RECALC=TRUE`) |
| `--mode cells\|reads\|both` | which cartesian sweep(s) to run (default `both`; required by power-curve fitting) |
| `--only A,B,C` | run only these stages |
| `--skip A,B` | skip these stages |
| `--from STAGE` | resume from this stage (skips earlier ones) |
| `--stop-on-fail` | abort on first stage failure (default: log and continue) |
| `--dry-run` | print what would run without executing |

Stages run **strictly serially** in DAG order:

1. `sharding` (5 dataset shard scripts in sequence)
2. `eda_roussos`
3. `full_ground_truth`
4. `bulk_deg_processing`
5. `cartesian_downsampling` (`cells` then `reads` if `--mode both`)
6. `depth_equalization`
7. `power_curve_fitting`
8. `de_diagnostics`
9. `figure_data_generation`

Each invocation creates a timestamped log directory under
`Datalink/Run_Output_Logs/run_all/<UTC-timestamp>/` containing:

- `00_run.log` — master log
- `<stage>.log` — per-stage stdout/stderr
- `<stage>.status` — `SUCCESS`, `FAILED`, or `SKIPPED`
- `manifest.json` — invocation context (host, git SHA, args)

A `latest` symlink points to the most recent run for convenience.
A lock file (`.run_all.lock`) prevents concurrent runs from clobbering each
other.

### How `--no-overwrite` actually works

Each pipeline's `config.R` exposes a `FORCE_RECALC` flag controlling whether
that pipeline regenerates outputs that already exist. `run_all.sh` exports
`POWER_ANALYSIS_FORCE_RECALC` as an environment variable, which every config
respects:

```r
FORCE_RECALC <- if (nzchar(Sys.getenv("POWER_ANALYSIS_FORCE_RECALC"))) {
  tolower(Sys.getenv("POWER_ANALYSIS_FORCE_RECALC")) %in% c("true", "1", "yes")
} else {
  TRUE
}
```

So `--no-overwrite` is equivalent to setting `POWER_ANALYSIS_FORCE_RECALC=false`
before launching `Rscript`. The per-pipeline skip-if-exists logic is
unchanged — `run_all.sh` is just a thin wrapper that flips the flag.

## 2.6 Fully containerised

If you don't want to install R, Bioconductor, Python or anything else on the
host, the Dockerfile in this repo bakes the entire compute environment into a
reproducible image. The workflow is four small scripts in `docker/`, run from the
repo root — `docker/docker-build.sh`, `docker/docker-smoke.sh`,
`docker/docker-run.sh`, and a `docker/docker-all.sh` wrapper that chains
build → smoke → run with per-phase skips.

```bash
# Fresh: build + smoke-test + run PCF (the default pipeline):
./docker/docker-all.sh

# Full DAG instead of just PCF:
./docker/docker-all.sh --pipeline all --max-cores 24

# Build + smoke-test only, don't run:
./docker/docker-all.sh --skip-run

# Already built? Skip ahead and just run:
./docker/docker-all.sh --skip-build --skip-smoke

# Long unattended run:
nohup ./docker/docker-all.sh --pipeline all --max-cores 24 > nohup.docker-all.out 2>&1 & disown
```

The container runs as your host UID/GID, so files written into the
bind-mounted `Datalink/` are owned by your host user, not root. `docker-run.sh`
writes a provenance `manifest.json` and a tee'd `run.log` under
`logs/docker-run/<UTC-timestamp>/` (with a `latest` symlink).

**See [`docker.md`](docker.md) for the full reference** — every script's flags,
the runtime model, and troubleshooting.

## 3. Pipeline dependency graph

```
raw_data
   │
   ▼
[A] sharding ──────────────────────────────────┐
   │                                           │
   ├──► [B] eda_roussos        (independent)   │
   ├──► [C] full_ground_truth ────────┐        │
   │        │                         │        │
   │        ├──► [G] de_diagnostics   │        │
   │        └──► [H] bulk_deg_processing       │
   │                                           │
   ├──► [D] cartesian_downsampling ◄───────────┘
   │        │     (needs Shards + GT; multi-day)
   │        └──► [F] power_curve_fitting
   │
   └──► [E] depth_equalization (independent, can run anytime after A)

All of C, D (via F), F, G, H feed into:
[I] figure_data_generation  (runs last; minutes)
```

### Suggested schedule for a full end-to-end run

| Day | Stage | Notes |
| --- | --- | --- |
| 0 | Setup (sections 0-2) | ~1 h |
| 1 (morning) | A (sharding), C (full_ground_truth) | ~3-5 h |
| 1 (afternoon) | Start D (cartesian_downsampling) in `screen`, detach | multi-day |
| 1-2 (parallel) | B, E, G, H in other windows | a few hours each |
| 3+ | F (power_curve_fitting), I (figure_data_generation) | ~1 h total |

You can run B, E, G, H in parallel with D — none of them touch the
cartesian working directory.

### `run_pipeline.R` conventions

Every pipeline under `R/*/` exposes a `run_pipeline.R` that takes two
positional arguments:

```bash
Rscript R/<pipeline>/run_pipeline.R [steps] [cores]
# Examples
Rscript R/cartesian_downsampling/run_pipeline.R all 16
Rscript R/cartesian_downsampling/run_pipeline.R 3,7 8
```

`steps` defaults to `all`; `cores` defaults to a sensible number. Reduce
`cores` if you hit OOM.

---

## A. Sharding

Convert heterogeneous raw files into standardised donor-level shards.

- **Depends on:** raw data (section 1)
- **Produces:** `Datalink/Shards/{CellType}/{Dataset}/*.qs`
- **Runtime:** 1-3 h total (I/O bound)
- **Required for:** every downstream pipeline

```bash
cd "$POWER_ANALYSIS_ROOT"
screen -S A_shard
Rscript R/sharding/Run_Shard_Roussos.R         2>&1 | tee logs/A1_roussos.log
Rscript R/sharding/Run_Shard_Allen.R           2>&1 | tee logs/A2_allen.log
Rscript R/sharding/Run_Shard_DeJager.R         2>&1 | tee logs/A3_dejager.log
Rscript R/sharding/Run_Shard_Johansen.R        2>&1 | tee logs/A4_johansen.log
Rscript R/sharding/Run_Shard_Rexach.R          2>&1 | tee logs/A5_rexach.log
Rscript R/sharding/Run_Shard_Roussos_Full.R    2>&1 | tee logs/A6_roussos_full.log
Rscript R/sharding/Run_Harmonise_Gene_IDs.R    2>&1 | tee logs/A7_harmonise.log
```

**Verify:**
```bash
ls Datalink/Shards/Astro/Roussos | head
find Datalink/Shards -name '*.qs' | wc -l   # expect 1000+ across 5 datasets
```

---

## B. EDA — Roussos

Exploratory PCA gallery and filtering summary.

- **Depends on:** A
- **Produces:** `Datalink/Exploratory_Analysis/`
- **Runtime:** minutes
- **Downstream:** none (cosmetic / sanity check)

```bash
Rscript R/eda_roussos/run_pipeline.R 2>&1 | tee logs/B_eda.log
```

**Verify:**
```bash
ls Datalink/Exploratory_Analysis/
git diff --stat Datalink/Exploratory_Analysis/
```

---

## C. Full ground truth (DE across 5 datasets × 3 cell types)

DESeq2 with dataset-specific covariates + apeglm LFC shrinkage. 15
DE analyses, plus cross-dataset diagnostics (UpSet, LFC correlation,
concordance, forest plots).

- **Depends on:** A
- **Produces:** `Datalink/Full_Ground_Truth/`
- **Runtime:** 1-2 h
- **Required for:** D, G, H, I

```bash
Rscript R/full_ground_truth/run_pipeline.R 2>&1 | tee logs/C_gt.log
```

**Verify:**
```bash
ls Datalink/Full_Ground_Truth/
head -3 Datalink/Full_Ground_Truth/FDP_Summary.csv
git diff --no-index -- \
  <(git show HEAD:Datalink/Full_Ground_Truth/FDP_Summary.csv) \
  Datalink/Full_Ground_Truth/FDP_Summary.csv | head -30
```

DEG counts should match: Astro 34, Micro 29, Oligo 44 at padj<0.05,
|LFC|>1 on Roussos.

---

## D. Cartesian downsampling (THE big one)

Empirical power grid: 12 × 14 (donors × cells) × 20 iterations ×
6 robustness runs, in cells-based and reads-based modes.

- **Depends on:** A, C
- **Produces:** `Datalink/Cartesian_Downsampling/`
- **Runtime:** **many hours to days** depending on cores. Use `screen`.
- **Required for:** F, I

```bash
screen -S D_downsample
export POWER_ANALYSIS_ROOT="$REPRO/Empirical_Power_Analysis"
cd "$POWER_ANALYSIS_ROOT"
Rscript R/cartesian_downsampling/run_pipeline.R all 16 2>&1 | tee logs/D_downsampling.log
# Ctrl-a d to detach. Reattach later: screen -r D_downsample
```

**Verify (after completion):**
```bash
find Datalink/Cartesian_Downsampling -name 'results*.csv' | wc -l
git diff --stat Datalink/Cartesian_Downsampling/
```

**While it runs,** you can check progress with:
```bash
tail -f logs/D_downsampling.log
find Datalink/Cartesian_Downsampling -newer logs/D_downsampling.log -type f | wc -l
```

---

## E. Depth equalisation (optional sensitivity)

Multinomial downsampling of Astro reads to Micro-equivalent depth,
then DE.

- **Depends on:** A
- **Produces:** `Datalink/Depth_Equalization/`
- **Runtime:** 1-2 h
- **Required for:** I (only for one supplementary figure)

```bash
Rscript R/depth_equalization/run_pipeline.R 2>&1 | tee logs/E_depth.log
```

**Verify:** Expect ShallowFewer at |LFC|>1 to show -20.6% DEG reduction
(34 → 27 for Astro).

---

## F. Power-curve fitting

Logistic regression power equations, per cell type and combined.
Produces the standalone `predict_power.R` function + coefficient /
lookup tables.

- **Depends on:** D
- **Produces:** `Datalink/Power_Curve_Fitting/`, `Datalink/Power_Curve_Fitting/Exports/predict_power.R`
- **Runtime:** ~30 min
- **Required for:** I

```bash
Rscript R/power_curve_fitting/run_pipeline.R 2>&1 | tee logs/F_pcf.log
```

**Verify:**
```bash
git diff --no-index -- \
  <(git show HEAD:Datalink/Power_Curve_Fitting/Exports/power_equation_coefficients.csv) \
  Datalink/Power_Curve_Fitting/Exports/power_equation_coefficients.csv | head -40
Rscript -e 'source("Datalink/Power_Curve_Fitting/Exports/predict_power.R"); print(predict_power(n_donors=25, n_cells=250, reads_per_cell=4000, lfc=1.5, baseMean=50, dispersion=0.5, cell_type="Microglia"))'
```

Coefficient values should agree with the shipped table to 3-4 d.p.

---

## G. DE diagnostics

Density plots, PCA validation, split recovery, chromosome grids,
shared-gene diagnostics.

- **Depends on:** C
- **Produces:** `Datalink/DE_Diagnostics/`
- **Runtime:** ~30 min
- **Required for:** nothing downstream (informational)

```bash
Rscript R/de_diagnostics/run_pipeline.R 2>&1 | tee logs/G_diag.log
```

**Verify:**
```bash
ls Datalink/DE_Diagnostics/
git diff --stat Datalink/DE_Diagnostics/
```

---

## H. Bulk DEG processing (GTEx comparison)

Processes GTEx v8 DLPFC bulk RNA-seq, computes sex-stratified DE,
concordance with single-cell.

- **Depends on:** C, plus GTEx raw files in `Datalink/raw_data/`
- **Produces:** `Datalink/Bulk_Analysis/`
- **Runtime:** ~1 h (biomaRt step can be slow if ENSEMBL is flaky)
- **Required for:** nothing downstream

```bash
Rscript R/bulk_deg_processing/run_pipeline.R 2>&1 | tee logs/H_bulk.log
```

**Verify:**
```bash
ls Datalink/Bulk_Analysis/
git diff --stat Datalink/Bulk_Analysis/
```

If `biomaRt` fails, rerun later; the pipeline caches successful
responses.

---

## I. Figure data generation

Turns the outputs of C, D, E, F into the canonical CSVs the paper
imports (`read.csv()` from `Datalink/Paper_Source_Data/`).

- **Depends on:** C, F (and optionally E for the depth supplementary)
- **Produces:** `Datalink/Paper_Source_Data/`
- **Runtime:** minutes
- **Required for:** anything that reads the paper source tables

```bash
Rscript R/figure_data_generation/run_pipeline.R 2>&1 | tee logs/I_figdata.log
```

**Verify:**
```bash
ls Datalink/Paper_Source_Data/
for f in \
  Datalink/Paper_Source_Data/dataset_statistics.csv \
  Datalink/Paper_Source_Data/fdp_summary.csv \
  Datalink/Paper_Source_Data/model_coefficients.csv \
  Datalink/Paper_Source_Data/truth_gene_characteristics.csv \
  Datalink/Paper_Source_Data/model_summary.csv ; do
  echo "=== $f ==="
  git diff --no-index -- <(git show HEAD:"$f") "$f" | head -20
done
```

---

## J. Regenerating composed figures

Builds the four composite PNGs that the main paper draft reads from
`Datalink/Composed_Figures/` (`figure_reproducibility_column.png`,
`figure_bulk_validation.png`, `figure_cartesian_oligo_recall.png`,
`figure_cartesian_oligo_precision.png`). Source PNGs are produced by
sections C, D and H.

- **Depends on:** C, D, H (the source heatmaps and bar charts must
  already exist under `Datalink/`)
- **Produces:** `R/figure_composition/output/{v1_tight,v2_spaced,v3_bordered}/`
  and the four PNGs that should then be copied into
  `Datalink/Composed_Figures/`
- **Runtime:** seconds
- **Prerequisites:** Python 3.10+, Pillow (`pip install Pillow`). Fonts
  resolve from system DejaVu / Liberation / FreeSans (bold sans-serif),
  with a fallback to PIL's default bitmap font.

```bash
cd R/figure_composition && python3 compose_figures.py
```

To promote a fresh render into the paper's data path:

```bash
cp R/figure_composition/output/v1_tight/figure_reproducibility_column.png \
   R/figure_composition/output/v1_tight/figure_bulk_validation.png \
   R/figure_composition/output/v1_tight/figure_cartesian_oligo_recall.png \
   R/figure_composition/output/v1_tight/figure_cartesian_oligo_precision.png \
   Datalink/Composed_Figures/
```

---

## K. Regenerating the three cross-dataset heatmaps

Lightweight viz-only rerun for the three pheatmaps in `FigRepro` (raw IoU,
universe-controlled IoU, FDP). Reads the intermediate CSVs that section C
already persists, redraws only the heatmaps, and skips any DE or set
operations.

- **Depends on:** C (`R/full_ground_truth/02_Cross_Dataset_Diagnostics.R`
  must have been run at least once to produce the intermediate CSVs)
- **Produces:** updated `Expanded_IoU_Matrix.png`,
  `Expanded_IoU_Matrix_Universe_Controlled.png`, `FDP_Heatmap.png` under
  `Datalink/Full_Ground_Truth/`
- **Runtime:** under 2 minutes

```bash
Rscript R/full_ground_truth/02b_Regenerate_Heatmaps.R
```

---

## 4. Whole-repo output comparison

After the full rerun:

```bash
cd "$REPRO/Empirical_Power_Analysis"
git status --short Datalink/
git diff --stat Datalink/
```

PNG figures will differ byte-wise across BLAS / font stacks but should
match in shape, axis ranges, and overlay text. CSV numeric content
should agree to 3-4 d.p.

Numeric tolerance check for the headline tables:

```bash
Rscript -e '
files <- c("Datalink/Paper_Source_Data/model_coefficients.csv",
           "Datalink/Paper_Source_Data/dataset_statistics.csv",
           "Datalink/Power_Curve_Fitting/Exports/power_equation_coefficients.csv")
for (f in files) {
  cat("---", f, "---\n")
  old <- read.csv(pipe(sprintf("git show HEAD:%s", f)), stringsAsFactors=FALSE)
  new <- read.csv(f, stringsAsFactors=FALSE)
  cat("dims old:", dim(old), "  new:", dim(new), "\n")
  if (all(dim(old) == dim(new))) {
    num <- sapply(old, is.numeric)
    diffs <- mapply(function(a,b) max(abs(a-b), na.rm=TRUE), old[num], new[num])
    print(round(diffs, 6))
  }
}'
```

## 5. Capture environment snapshot

```bash
Rscript R/shared/capture_environment.R
diff -r Reproducibility/ "$REPRO/Empirical_Power_Analysis/Reproducibility/" | head
```

---

## Partial-reproduction recipes

| I only want to... | Run sections |
| --- | --- |
| Check the DESeq2 ground-truth DEG counts | A, C |
| Rebuild the power equations and `predict_power.R` | A, C, D, F |
| Reproduce the depth-saturation sensitivity | A, E |
| Check paper source-data CSVs only (cheapest full-loop) | A, C, D, F, I |
| Reproduce the GTEx bulk concordance | A, C, H |
| Cross-dataset QC diagnostics only | A, C, G |

---

## Troubleshooting

- **`here::here()` picks the wrong root** — `export POWER_ANALYSIS_ROOT` before running.
- **`biomaRt` failing** — ENSEMBL's public server is periodically unavailable; the relevant pipeline (H) caches results once a successful run completes.
- **Out-of-memory in cartesian_downsampling** — reduce the core count argument (`Rscript ... all 4` instead of `16`).
- **Johansen `.tar` extraction errors** — ensure `tar` is on `PATH`; the sharding helper shells out.
- **Raw data missing for one cohort** — pipelines skip absent datasets with a warning; cross-dataset diagnostics will simply have fewer comparisons.
- **Screen session died** — cartesian_downsampling is idempotent per step; rerun `Rscript R/cartesian_downsampling/run_pipeline.R <remaining_steps> 16`.

---

> This file is intended to stay in sync with `README.md`. If you add a
> new pipeline or change the run order, update both.

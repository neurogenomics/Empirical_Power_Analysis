# Agent Guide

This repository is designed to be navigable by AI coding agents (Claude Code, Cursor,
OpenHands, Aider, etc.) as well as by humans. This document is the canonical entry
point for an agent arriving at this repository fresh.

This file is also the `CLAUDE.md` / `OPENAI_AGENTS.md` equivalent — rather than
fragmenting discovery across tool-specific files, all agent-facing documentation lives
here. Tool-specific symlinks (`CLAUDE.md` → `AGENTS.md`) are provided for convenience.

---

## TL;DR for agents

- **One command to reproduce the paper (containerised, unattended)**:
  `nohup ./docker-all.sh --pipeline all --max-cores 24 > nohup.docker-all.out 2>&1 & disown`
  (build → smoke → run; see [`DOCKER.md`](DOCKER.md))
- **One command for power-curve fitting only** (host-native bash):
  `./run_power_curve_fitting.sh --max-cores 24 --no-overwrite`
- **One command for the full DAG** (host-native bash):
  `./run_all.sh --max-cores 24 --no-overwrite`
- **Preferred edit surfaces**: `R/<pipeline>/` (analysis), `run_all.sh` /
  `run_power_curve_fitting.sh` (orchestration), `Dockerfile` + `docker/` +
  `docker-build.sh`/`docker-smoke.sh`/`docker-run.sh`/`docker-all.sh`
  (containerisation), `paper/paper.Rmd` (manuscript)
- **Read-only surfaces**: `Datalink/` (generated outputs — regenerate via
  pipeline), `Reproducibility/` (environment snapshots)
- **Entry points by intent**:

| If the user wants to... | Start reading here |
|--|--|
| Understand the scientific claim | `paper/paper.Rmd`, `README.md` § Key findings |
| Change an analysis step | `R/<pipeline>/run_pipeline.R` + per-pipeline `config.R` |
| Change workflow orchestration | `run_all.sh` (DAG order, stage skip/only, FORCE_RECALC plumbing) |
| Change container build | `Dockerfile`, `docker/install_packages.R`, `docker/entrypoint.sh` |
| Change container build/run flow | `docker-build.sh` / `docker-smoke.sh` / `docker-run.sh` / `docker-all.sh` (see `DOCKER.md`) |
| Add a new dataset | `R/sharding/Run_Shard_<Dataset>.R` + add its key to `SHARDING_DATASETS` in `run_all.sh` |
| Bump a dependency | `renv.lock` via `renv::update()` + rebuild Docker image |

---

## Architecture at a glance

```
raw data (h5ad, qs)              ← not in git (dataset access terms; pointers in README)
    │
    ▼
R/sharding/                      ← donor-level shards (memory-intensive, single-pass)
    │
    ├─→ R/eda_roussos/          ← PCA + filter diagnostics (parallel branch)
    │
    └─→ R/full_ground_truth/    ← DESeq2 DE per cell type × dataset
            │
            ├─→ R/cartesian_downsampling/   ← 120k DESeq2 runs (the hot loop)
            ├─→ R/depth_equalization/       ← depth-matched astro sensitivity
            ├─→ R/de_diagnostics/           ← chromosome grid + shared-gene checks
            │
            │  (parallel, independent)
            └─→ R/bulk_deg_processing/      ← GTEx bulk concordance
                    │
                    ▼
            R/power_curve_fitting/          ← logistic fit on 22M binary outcomes
                    │
                    ▼
            R/figure_data_generation/       ← canonical CSVs consumed by paper
                    │
                    ▼
            paper/paper.Rmd                 ← rendered manuscript
```

The same DAG is encoded in `run_all.sh` (strict serial, single host). Agents
editing the workflow should keep `run_all.sh`'s stage list in sync with the
R pipelines under `R/<pipeline>/`.

---

## Interface contracts

Each `R/<pipeline>/run_pipeline.R` exposes a **stable CLI**:

```bash
Rscript R/<pipeline>/run_pipeline.R [steps] [n_cores]
#   steps    = "all" or comma-separated step numbers (e.g. "1,2,5")
#   n_cores  = integer for within-pipeline future.apply parallelism
```

This contract is what `run_all.sh` depends on. If you change a
`run_pipeline.R`'s CLI, update the corresponding `stage_*` function in
`run_all.sh`. The contract is intentionally minimal to keep the coupling
shallow.

Each pipeline sources its own `config.R`, which resolves `POWER_ANALYSIS_ROOT`
via (in priority order):

1. Environment variable `POWER_ANALYSIS_ROOT`
2. `here::here()` fallback

`run_all.sh` exports `POWER_ANALYSIS_ROOT="${SCRIPT_DIR}"` and the Docker
entrypoint exports `POWER_ANALYSIS_ROOT=/workspace`, so pipelines run
identically host-native, container-native, or under any future workflow
manager.

---

## Testing conventions

Every pipeline directory has:

- `tests/` — unit tests via `testthat`, run with `Rscript -e "testthat::test_dir('tests')"`
- `toy/` — miniature synthetic fixtures exercising the full step sequence end-to-end

CI runs unit tests + a container build (`docker.yml`) on every PR. Agents
proposing changes should update tests in the same PR. **Do not merge
scientific-code changes without a corresponding test.**

---

## Scientific code safety rails

When modifying `R/<pipeline>/` code, respect these invariants:

1. **Deterministic seeding.** Every randomised step must call `set.seed()` with a
   seed derived from the robustness-run index + step index.  Breaking this silently
   makes results non-reproducible; CI has a determinism check that catches byte-level
   drift on toy data.
2. **Ground-truth / validation split integrity.** `GT_SPLIT <- 0.5` in config must
   be honoured; agents must not create analyses that leak information between splits.
3. **Filter thresholds are load-bearing.** `MIN_GENE_COUNT`, `FDR_THRESHOLD`, and
   the `LFC_INTERVALS` define the paper's claims.  Changes need a paper-text update.
4. **Cell-type order.** Astro, Micro, Oligo is the paper's reporting order.  Figures
   break if this is reshuffled.
5. **DESeq2 version.** Results are sensitive to DESeq2/apeglm updates.  The Docker
   image and `renv.lock` pin exact versions.  Don't bump without re-running the full
   Cartesian grid and comparing to the published reference.

---

## When not to use this pipeline

This is a **paper-specific analysis pipeline**, not a general-purpose scRNA-seq power
tool.  If a user wants to:

- Plan a study with their own data → point them at `predict_power.R` in
  `Datalink/Power_Curve_Fitting/Exports/` and the `power_lookup_table.csv`
- Build a new power-analysis method → fork and adapt `R/cartesian_downsampling/`
  as the empirical benchmarking harness
- Analyse their own scRNA-seq → use [nf-core/scrnaseq](https://nf-co.re/scrnaseq)
  for preprocessing, then this framework for power assessment

---

## Transparency about AI-assisted development

Portions of this repository were developed with AI coding assistance (Claude, GPT)
under human review. Specifically:

- `run_all.sh`, `run_power_curve_fitting.sh`, `docker-build.sh`,
  `docker-smoke.sh`, `docker-run.sh`, `docker-all.sh`, `docker/entrypoint.sh`
  — AI-scaffolded, human-reviewed
- `Dockerfile`, CI workflows — AI-scaffolded
- `R/**` scientific code — human-authored with AI-assisted refactoring at documented
  points (see commit history)
- Paper drafting review — see `ACKNOWLEDGEMENTS` in `paper/paper.Rmd`

The scientific claims, experimental design, and interpretation are the authors'.

We maintain this disclosure because we believe research software built in the
agentic-AI era should be transparent about tooling provenance — same way we disclose
funders and institutional affiliations.

---

## For agents making non-trivial changes

Before opening a PR, run locally:

```bash
# 1. Wrapper sanity check
./run_all.sh --dry-run
./run_power_curve_fitting.sh --dry-run

# 2. Unit tests
Rscript -e "for (d in list.dirs('R', recursive=FALSE)) if (dir.exists(file.path(d,'tests'))) testthat::test_dir(file.path(d,'tests'))"

# 3. Container build + smoke test (no full pipeline)
./docker-all.sh --skip-run

# 4. Lint
Rscript -e "lintr::lint_dir('R')"
```

The same checks run in CI — PRs will be blocked until they pass.

---

## Further reading

- `docs/reproducing-the-paper.md` — full reproduction walkthrough for humans
- `docs/architecture.md` — architecture decision records
- `docs/extending-to-new-data.md` — how to run this framework on your own cohort
- `CONTRIBUTING.md` — human-contributor guide
- `paper/paper.Rmd` — the manuscript itself

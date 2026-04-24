# Usage

## Introduction

This document describes how to run the **Empirical Power Analysis** pipeline.
For a detailed cold-start reproduction walkthrough (including dataset access,
dependency installation, and per-pipeline commands), see
[`REPRODUCE.md`](../REPRODUCE.md) in the repository root.

For the full pipeline architecture and design rationale, see
[`docs/architecture.md`](architecture.md).

---

## Quick start

### 1. Install Nextflow

```bash
# Requires Java 11+
curl -fsSL https://get.nextflow.io | bash
```

### 2. Pull the container

```bash
docker pull ghcr.io/neurogenomics/empirical_power_analysis:1.0.0
# or for Apptainer/Singularity on HPC:
apptainer pull empirical_power_analysis.sif \
    docker://ghcr.io/neurogenomics/empirical_power_analysis:1.0.0
```

### 3. Provide raw data

Place (or symlink) the raw `.h5ad`/`.qs` single-cell files in
`Datalink/raw_data/`. See [REPRODUCE.md §1](../REPRODUCE.md#1-provide-the-raw-data)
for the expected filenames per dataset.

### 4. Run

```bash
# HPC with SLURM + Apptainer (recommended)
nextflow run neurogenomics/Empirical_Power_Analysis \
    -profile slurm,apptainer \
    --outdir results

# Local with Docker (for testing)
nextflow run neurogenomics/Empirical_Power_Analysis \
    -profile test,docker

# Stub run (no compute needed — validates pipeline wiring)
nextflow run neurogenomics/Empirical_Power_Analysis \
    -profile test -stub-run
```

---

## Parameters

The pipeline parameters are documented in `nextflow_schema.json`. Key parameters:

| Parameter | Default | Description |
|---|---|---|
| `raw_data` | `Datalink/raw_data` | Path to raw single-cell files |
| `shards_dir` | `Datalink/Shards` | Where donor-level shards are written |
| `downsample_mode` | `cells` | `cells` or `reads` (Cartesian grid axis) |
| `n_robustness` | `6` | Number of independent ground-truth splits |
| `n_subsamples` | `20` | Subsamples per grid point |
| `render_paper` | `true` | Whether to render `paper/paper.Rmd` |
| `outdir` | `results` | Output directory |

Override any parameter on the CLI: `--n_robustness 3` or in a custom config file.

---

## Profiles

| Profile | Description |
|---|---|
| `standard` | Local executor, no container |
| `docker` | Local executor with Docker |
| `apptainer` | HPC executor with Apptainer/Singularity |
| `slurm` | SLURM executor (Imperial CX3 defaults) |
| `conda` | Conda environment (fallback) |
| `test` | Minimal toy fixtures, for CI and quick validation |
| `test_full` | Full data, all grid points (multi-day) |
| `debug` | Verbose logs, no cleanup |

---

## Running step subsets

Each R pipeline accepts a `steps` argument. To run only specific steps via
Nextflow, override the relevant Nextflow process input or run the R scripts
directly — see [REPRODUCE.md](../REPRODUCE.md) for per-pipeline instructions.

---

## Resource requirements

| Step | CPUs | RAM | Wall-time |
|---|---|---|---|
| Sharding | 8 | 200 GB | 3–12 h |
| Full ground truth | 16 | 64 GB | 1–2 h |
| Cartesian downsampling | 64 | 200 GB | 2–5 days |
| Power curve fitting | 16 | 64 GB | ~30 min |
| All others | 2–16 | 16–64 GB | minutes–hours |

---

## Troubleshooting

See [`REPRODUCE.md` §Troubleshooting](../REPRODUCE.md#troubleshooting) for
common issues (biomaRt failures, OOM in cartesian step, here::here() root,
etc.).

# Sharding

Converts the raw PsychAD snRNA-seq files (`Roussos` in code and paths) into a
uniform donor-level shard store. Every downstream pipeline (full ground truth,
cartesian downsampling, depth equalisation, EDA) expects shards in the format
produced here.

The directory also contains sharders for four other cohorts (SEA-AD, ROSMAP,
Johansen, Rexach) and a cross-dataset gene ID harmoniser. These were used for
the earlier cross-dataset analysis, which is no longer part of the paper; they
are obsolete (not used for the paper).

## Use in the Reproduction Pipeline

The canonical end-to-end route is `pipeline/` (see `pipeline/README.md`). Stage
`shard` first checks the raw inputs against `pipeline/inputs.sha256`, then runs
only:

```bash
Rscript R/sharding/Run_Shard_Roussos.R Astro,Micro,Oligo
```

## What each script does

| Script                         | Status   | Purpose |
|--------------------------------|----------|---------|
| `Run_Shard_Roussos.R`          | Paper    | Orchestrator: runs the three steps below for each cell type |
| `Run_Shard_Roussos_Full.R`     | Paper    | Step 1: unfiltered donor shards from `PsychAD-{HBCC,RADC,MSSM}.h5ad` (`Roussos_Full`) |
| `Annotate_Roussos_PMI.R`       | Paper    | Step 3: writes donor post-mortem interval (`pmi`, minutes) from `media-1.csv` into shard `colData` |
| `Audit_Roussos_PMI.R`          | Optional | Read-only check that every shard donor joins to a PMI value; writes `Datalink/Shards/.pmi_audit/` |
| `Run_Shard_Allen.R`            | Obsolete (not used for the paper) | SEA-AD sharder |
| `Run_Shard_DeJager.R`          | Obsolete (not used for the paper) | ROSMAP sharder |
| `Run_Shard_Johansen.R`         | Obsolete (not used for the paper) | Johansen sharder |
| `Run_Shard_Rexach.R`           | Obsolete (not used for the paper) | Rexach sharder |
| `Run_Harmonise_Gene_IDs.R`     | Obsolete (not used for the paper) | Converts gene symbols to Ensembl IDs in DeJager/Johansen shards |

`Run_Shard_Roussos.R` steps, per cell type:

1. Runs `Run_Shard_Roussos_Full.R` (skipped if the `Roussos_Full` manifest exists).
2. Filters `Roussos_Full` donors to age >= 40 years and >= 30 cells of the
   target cell type, copies passing shards to `Roussos/`, and writes the
   manifest and a `.filter_complete` sentinel.
3. Runs `Annotate_Roussos_PMI.R` (writes a `.pmi_annotated` sentinel per
   directory). PMI enters the DE design as `pmi_scaled`.

## Running

```bash
export POWER_ANALYSIS_ROOT=/path/to/this/repo     # optional; otherwise here::here()
Rscript R/sharding/Run_Shard_Roussos.R             # Astro, Micro, Oligo
Rscript R/sharding/Run_Shard_Roussos.R Astro,Oligo # subset
```

## Inputs

Under `Datalink/raw_data/` (not distributed): `PsychAD-HBCC.h5ad`,
`PsychAD-MSSM.h5ad`, `PsychAD-RADC.h5ad` and `media-1.csv`. See
`pipeline/README.md` for sources and `pipeline/inputs.sha256` for checksums.

## Outputs

Under `Datalink/Shards/{CellType}/Roussos_Full/` (all donors, used for the
unfiltered PCA) and `Datalink/Shards/{CellType}/Roussos/` (analysis cohort):

- `donor_<ID>.qs`: one `SingleCellExperiment` per donor, with `colData$donor_id`
  and `colData$pmi`.
- `00_MASTER_MANIFEST.qs`: list of
  - `manifest` (data.frame: `donor_id`, `cell_count`, `file_path`)
  - `feasibility` (donors passing each `CELL_STEPS` threshold)
  - `config` (record of filters applied)
- `.filter_complete` (filtered directory only) and `.pmi_annotated` sentinels.

## Configuration

`config.R` holds `MIN_CELLS_PER_DONOR` (30), `MIN_AGE_YEARS` (40) and
`CELL_STEPS`. `Run_Shard_Roussos.R` declares its own matching thresholds
(`MIN_CELLS_THRESHOLD`, `MIN_AGE_THRESHOLD`).

## Tests

`tests/test_sharding.R` checks existing shards: loadable SCE objects, required
metadata columns, integer counts, plausible dimensions, and manifest consistency
with shard files.

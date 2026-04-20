# Sharding

Standardises five heterogeneous snRNA-seq cohorts into a uniform donor-level
shard store. Every downstream pipeline (ground truth, cartesian downsampling,
depth equalization, glial characteristics) expects shards in the format produced
here.

## What each script does

| Script                         | Dataset            | Input (under `Datalink/raw_data/`)            |
|--------------------------------|--------------------|-----------------------------------------------|
| `Run_Shard_Allen.R`            | SEA-AD             | `allen_{astrocytes,microglia,oligo}.qs`       |
| `Run_Shard_DeJager.R`          | ROSMAP             | `dejager_*.h5Seurat`, `DeJagerMeta.csv`       |
| `Run_Shard_Johansen.R`         | Johansen           | `johansen_*.qs`                               |
| `Run_Shard_Rexach.R`           | Rexach             | `rexach_*.h5ad`                               |
| `Run_Shard_Roussos.R`          | PsychAD (subset)   | `roussos_*.h5ad`                              |
| `Run_Shard_Roussos_Full.R`     | PsychAD (full)     | `roussos_*.h5ad`                              |
| `Run_Harmonise_Gene_IDs.R`     | Cross-dataset      | Produced shards (maps Ensembl ↔ symbol)        |

Each script loads the raw object, parses and filters age, drops donors with
fewer than `MIN_CELLS_PER_DONOR` cells, writes one `.qs` per donor under
`Datalink/Shards/{CellType}/{Dataset}/`, and emits a `00_MASTER_MANIFEST.qs`
listing donors and feasibility for a range of cell-count levels.

## Running

The scripts source `config.R` automatically to pick up `PROJECT_ROOT` and the
shared filter constants.

```bash
# Activate the repo root however suits your setup:
export POWER_ANALYSIS_ROOT=/path/to/this/repo     # optional
Rscript R/sharding/Run_Shard_Allen.R
Rscript R/sharding/Run_Shard_Roussos.R
# ... repeat for the cohorts you have data for ...
Rscript R/sharding/Run_Harmonise_Gene_IDs.R        # after all shards exist
```

## Inputs

Raw single-cell files. Not distributed — see `Datalink/raw_data/README.md` for
dataset provenance and expected filenames.

## Outputs

Per dataset/cohort, under `Datalink/Shards/{CellType}/{Dataset}/`:

- `donor_<ID>.qs` — one `SingleCellExperiment` per donor, with `colData$donor_id`
  standardised across cohorts.
- `00_MASTER_MANIFEST.qs` — list of:
  - `manifest` (data.frame: `donor_id`, `cell_count`, `file_path`)
  - `feasibility` (named list of donors passing each `CELL_STEPS` threshold)
  - `config` (record of filters applied)

## Configuration

Shared filters live in `config.R` (`MIN_CELLS_PER_DONOR`, `MIN_AGE_YEARS`,
`CELL_STEPS`). Candidate donor and age column names are declared at the top of
each `Run_Shard_*.R` script, because they vary by dataset.

## Tests

`tests/test_sharding.R` exercises the shared age-parsing helper on a small
synthetic input.

# Raw data (not distributed)

This directory is intentionally empty. It is the expected location for the primary
single-nucleus RNA-seq input files consumed by `R/sharding/`.

Populate this directory from the original public sources before running the
sharding pipeline:

| Dataset       | Expected files                                            | Source |
|---------------|-----------------------------------------------------------|--------|
| PsychAD       | `roussos_<CellType>.h5ad`                                 | Fullard et al. 2025 (Synapse: syn52293417) |
| SEA-AD (Allen)| `allen_astrocytes.qs`, `allen_microglia.qs`, `allen_oligo.qs` | Allen Institute (sea-ad.org) |
| ROSMAP (DeJager)| `dejager_<CellType>.h5Seurat`, `DeJagerMeta.csv`        | Green et al. 2024 (Synapse: syn53366818) |
| Johansen      | `johansen_<CellType>.qs`                                  | Johansen et al. 2023 |
| Rexach        | `rexach_<CellType>.h5ad`                                  | Rexach et al. 2024 |

Exact filenames expected by each `Run_Shard_*.R` script are declared at the top
of that script (see the `INPUT_FILES` list).

Each dataset has its own access terms; please refer to the original study for
licensing and download instructions.

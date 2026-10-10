# Raw data (not distributed)

This directory is intentionally empty. It is the expected location for the raw
inputs when running the pipelines by hand (the `pipeline/` route stages them
itself). We do not redistribute them: they are already public, and the pipeline
regenerates every intermediate from them.

| Input | Source |
|-------|--------|
| `PsychAD-HBCC.h5ad`, `PsychAD-MSSM.h5ad`, `PsychAD-RADC.h5ad` | PsychAD snRNA-seq, [CELLxGENE collection 84ce6837](https://cellxgene.cziscience.com/collections/84ce6837-548d-4a1f-919f-0bc0d9a3952f) (Fullard et al. 2025) |
| `media-1.csv` | PsychAD donor table with post-mortem interval, a supplementary table of Fullard et al. 2025 |

Checksums for all four files are in `pipeline/inputs.sha256`. Exact filenames
expected by `R/sharding/Run_Shard_Roussos.R` are declared at the top of that
script. Please refer to the original study for dataset licensing terms.

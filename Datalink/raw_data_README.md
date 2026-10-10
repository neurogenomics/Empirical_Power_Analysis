# Raw data (not distributed)

The raw inputs are not distributed with this repository. `pipeline/run.sh`
reads them from `WORK/raw/` and mounts that directory read-only at
`Datalink/raw_data/` inside the container, which is where the sharding code
looks for them. We do not redistribute them: they are already public, and the
pipeline regenerates every intermediate from them.

| Input | Source |
|-------|--------|
| `PsychAD-HBCC.h5ad`, `PsychAD-MSSM.h5ad`, `PsychAD-RADC.h5ad` | PsychAD snRNA-seq, [CELLxGENE collection 84ce6837](https://cellxgene.cziscience.com/collections/84ce6837-548d-4a1f-919f-0bc0d9a3952f) (Fullard et al. 2025) |
| `media-1.csv` | PsychAD donor table with post-mortem interval, a supplementary table of Fullard et al. 2025 |

Checksums for all four files are in `pipeline/inputs.sha256`, and stage
`shard` checks them before it starts. The filenames are fixed in
`R/sharding/Run_Shard_Roussos_Full.R` (the `.h5ad` files) and
`R/sharding/Annotate_Roussos_PMI.R` (`media-1.csv`). Please refer to the original study for dataset licensing terms.

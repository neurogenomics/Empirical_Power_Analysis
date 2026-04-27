#!/bin/bash -ue
export POWER_ANALYSIS_ROOT=/mnt/data/chrisc/Power_Analysis_Paper/Empirical_Power_Analysis
mkdir -p Datalink/Shards

# Run each dataset's sharding script in sequence
Rscript /mnt/data/chrisc/Power_Analysis_Paper/Empirical_Power_Analysis/R/sharding/Run_Shard_Roussos.R
Rscript /mnt/data/chrisc/Power_Analysis_Paper/Empirical_Power_Analysis/R/sharding/Run_Shard_Allen.R
Rscript /mnt/data/chrisc/Power_Analysis_Paper/Empirical_Power_Analysis/R/sharding/Run_Shard_DeJager.R
Rscript /mnt/data/chrisc/Power_Analysis_Paper/Empirical_Power_Analysis/R/sharding/Run_Shard_Johansen.R
Rscript /mnt/data/chrisc/Power_Analysis_Paper/Empirical_Power_Analysis/R/sharding/Run_Shard_Rexach.R

cat <<-END_VERSIONS > versions.yml
"SHARDING":
    R: $(R --version | head -1 | sed 's/R version //')
    DESeq2: $(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
    zellkonverter: $(Rscript -e 'cat(as.character(packageVersion("zellkonverter")))')
END_VERSIONS

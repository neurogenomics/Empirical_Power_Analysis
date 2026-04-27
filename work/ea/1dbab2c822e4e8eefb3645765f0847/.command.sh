#!/bin/bash -ue
export POWER_ANALYSIS_ROOT=/mnt/data/chrisc/Power_Analysis_Paper/Empirical_Power_Analysis
Rscript /mnt/data/chrisc/Power_Analysis_Paper/Empirical_Power_Analysis/R/bulk_deg_processing/run_pipeline.R all 16

cat <<-END_VERSIONS > versions.yml
"BULK_DEG_PROCESSING":
    R: $(R --version | head -1 | sed 's/R version //')
    biomaRt: $(Rscript -e 'cat(as.character(packageVersion("biomaRt")))')
END_VERSIONS

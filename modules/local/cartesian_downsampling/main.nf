process CARTESIAN_DOWNSAMPLING {
    tag "cartesian:${params.downsample_mode}"
    label 'process_heavy'

    input:
    path shards
    path ground_truth

    output:
    path "Datalink/Cartesian_Downsampling/**", emit: results
    path "versions.yml"                      , emit: versions

    script:
    """
    export POWER_ANALYSIS_ROOT=${projectDir}
    export DOWNSAMPLE_MODE=${params.downsample_mode}

    # The master runner internally parallelises 120k DESeq2 runs via future.apply.
    # On SLURM, request a full node (see process label process_heavy) and let the
    # R-level future scheduler manage within-node parallelism.
    Rscript ${projectDir}/R/cartesian_downsampling/run_pipeline.R all ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //')
        DESeq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
        future: \$(Rscript -e 'cat(as.character(packageVersion("future")))')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p Datalink/Cartesian_Downsampling/Astro Datalink/Cartesian_Downsampling/Micro Datalink/Cartesian_Downsampling/Oligo
    touch Datalink/Cartesian_Downsampling/Astro/astro_cartesian_stub.csv
    touch Datalink/Cartesian_Downsampling/Micro/micro_cartesian_stub.csv
    touch Datalink/Cartesian_Downsampling/Oligo/oligo_cartesian_stub.csv
    echo "stub" > versions.yml
    """
}

process SHARDING {
    tag "sharding"
    label 'process_memory_intensive'

    input:
    path raw_data

    output:
    path "Datalink/Shards/**"           , emit: shards
    path "versions.yml"                 , emit: versions

    script:
    """
    export POWER_ANALYSIS_ROOT=\$(pwd)
    mkdir -p Datalink/Shards

    # Run each dataset's sharding script in sequence
    Rscript ${projectDir}/R/sharding/Run_Shard_Roussos.R
    Rscript ${projectDir}/R/sharding/Run_Shard_Allen.R
    Rscript ${projectDir}/R/sharding/Run_Shard_DeJager.R
    Rscript ${projectDir}/R/sharding/Run_Shard_Johansen.R
    Rscript ${projectDir}/R/sharding/Run_Shard_Rexach.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //')
        DESeq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
        zellkonverter: \$(Rscript -e 'cat(as.character(packageVersion("zellkonverter")))')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p Datalink/Shards/Astro/Roussos Datalink/Shards/Micro/Roussos Datalink/Shards/Oligo/Roussos
    touch Datalink/Shards/Astro/Roussos/sample_shard.qs
    touch Datalink/Shards/Micro/Roussos/sample_shard.qs
    touch Datalink/Shards/Oligo/Roussos/sample_shard.qs
    echo "stub" > versions.yml
    """
}

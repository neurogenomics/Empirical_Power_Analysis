process EDA_ROUSSOS {
    tag "eda_roussos"
    label 'process_medium'

    input:
    path shards

    output:
    path "Datalink/Exploratory_Analysis/**", emit: eda
    path "versions.yml"                    , emit: versions

    script:
    """
    export POWER_ANALYSIS_ROOT=\$(pwd)
    Rscript ${projectDir}/R/eda_roussos/run_pipeline.R all ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p Datalink/Exploratory_Analysis
    touch Datalink/Exploratory_Analysis/eda_stub.csv
    echo "stub" > versions.yml
    """
}

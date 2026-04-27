process DE_DIAGNOSTICS {
    tag "de_diagnostics"
    label 'process_medium'

    input:
    path ground_truth

    output:
    path "Datalink/DE_Diagnostics/**", emit: diagnostics
    path "versions.yml"              , emit: versions

    script:
    """
    export POWER_ANALYSIS_ROOT=${projectDir}
    Rscript ${projectDir}/R/de_diagnostics/run_pipeline.R all ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p Datalink/DE_Diagnostics
    touch Datalink/DE_Diagnostics/de_diagnostics_stub.csv
    echo "stub" > versions.yml
    """
}

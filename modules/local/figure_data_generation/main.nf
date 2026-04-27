process FIGURE_DATA_GENERATION {
    tag "figure_data_generation"
    label 'process_low'

    input:
    path ground_truth
    path cartesian_results
    path coefficients
    path bulk_results

    output:
    path "Datalink/Paper_Source_Data/**", emit: source_data
    path "versions.yml"                 , emit: versions

    script:
    """
    export POWER_ANALYSIS_ROOT=${projectDir}
    Rscript ${projectDir}/R/figure_data_generation/run_pipeline.R all ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p Datalink/Paper_Source_Data
    touch Datalink/Paper_Source_Data/dataset_statistics.csv
    echo "stub" > versions.yml
    """
}

process DEPTH_EQUALIZATION {
    tag "depth_equalization"
    label 'process_high'

    input:
    path shards
    path ground_truth

    output:
    path "Datalink/Depth_Equalization/**", emit: results
    path "versions.yml"                  , emit: versions

    script:
    """
    export POWER_ANALYSIS_ROOT=${projectDir}
    Rscript ${projectDir}/R/depth_equalization/run_pipeline.R all ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p Datalink/Depth_Equalization
    touch Datalink/Depth_Equalization/depth_eq_stub.csv
    echo "stub" > versions.yml
    """
}

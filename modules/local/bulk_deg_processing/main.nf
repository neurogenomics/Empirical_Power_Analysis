process BULK_DEG_PROCESSING {
    tag "bulk_deg_processing"
    label 'process_medium'

    output:
    path "Datalink/Bulk_Analysis/**", emit: bulk_results
    path "versions.yml"             , emit: versions

    script:
    """
    export POWER_ANALYSIS_ROOT=\$(pwd)
    Rscript ${projectDir}/R/bulk_deg_processing/run_pipeline.R all ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //')
        biomaRt: \$(Rscript -e 'cat(as.character(packageVersion("biomaRt")))')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p Datalink/Bulk_Analysis
    touch Datalink/Bulk_Analysis/bulk_stub.csv
    echo "stub" > versions.yml
    """
}

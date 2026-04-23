process FULL_GROUND_TRUTH {
    tag "full_ground_truth"
    label 'process_high'

    input:
    path shards

    output:
    path "Datalink/Full_Ground_Truth/**", emit: ground_truth
    path "versions.yml"                 , emit: versions

    script:
    """
    export POWER_ANALYSIS_ROOT=\$(pwd)
    Rscript ${projectDir}/R/full_ground_truth/run_pipeline.R all ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //')
        DESeq2: \$(Rscript -e 'cat(as.character(packageVersion("DESeq2")))')
        apeglm: \$(Rscript -e 'cat(as.character(packageVersion("apeglm")))')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p Datalink/Full_Ground_Truth/Astro Datalink/Full_Ground_Truth/Micro Datalink/Full_Ground_Truth/Oligo
    touch Datalink/Full_Ground_Truth/Astro/astro_gt_stub.csv
    touch Datalink/Full_Ground_Truth/Micro/micro_gt_stub.csv
    touch Datalink/Full_Ground_Truth/Oligo/oligo_gt_stub.csv
    echo "stub" > versions.yml
    """
}

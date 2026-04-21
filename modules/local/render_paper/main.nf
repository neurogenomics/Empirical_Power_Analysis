process RENDER_PAPER {
    tag "render_paper"
    label 'process_low'
    publishDir "${params.outdir}/paper", mode: params.publish_mode

    input:
    path source_data

    output:
    path "paper.html", optional: true, emit: html
    path "paper.pdf" , optional: true, emit: pdf
    path "versions.yml"               , emit: versions

    script:
    """
    export POWER_ANALYSIS_ROOT=\$(pwd)
    # Render the manuscript from Datalink/ — figures and tables are read by relative path
    Rscript -e "rmarkdown::render('${projectDir}/paper_final_draft.rmd', output_dir=getwd())"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //')
        rmarkdown: \$(Rscript -e 'cat(as.character(packageVersion("rmarkdown")))')
    END_VERSIONS
    """

    stub:
    """
    echo "<html><body>stub paper</body></html>" > paper.html
    echo "stub" > versions.yml
    """
}

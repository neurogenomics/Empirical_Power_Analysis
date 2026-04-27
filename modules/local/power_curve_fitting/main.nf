process POWER_CURVE_FITTING {
    tag "power_curve_fitting"
    label 'process_high'

    input:
    path cartesian_results
    path depth_eq_results

    output:
    path "Datalink/Power_Curve_Fitting/**"                              , emit: fits
    path "Datalink/Power_Curve_Fitting/Exports/predict_power.R"         , emit: predict_fn
    path "Datalink/Power_Curve_Fitting/Exports/power_equation_coefficients.csv", emit: coefficients
    path "Datalink/Power_Curve_Fitting/Exports/power_lookup_table.csv"  , emit: lookup_table
    path "versions.yml"                                                 , emit: versions

    script:
    """
    export POWER_ANALYSIS_ROOT=${projectDir}
    Rscript ${projectDir}/R/power_curve_fitting/run_pipeline.R all ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //')
        lmtest: \$(Rscript -e 'cat(as.character(packageVersion("lmtest")))')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p Datalink/Power_Curve_Fitting/Exports
    touch Datalink/Power_Curve_Fitting/Exports/predict_power.R
    touch Datalink/Power_Curve_Fitting/Exports/power_equation_coefficients.csv
    touch Datalink/Power_Curve_Fitting/Exports/power_lookup_table.csv
    echo "stub" > versions.yml
    """
}

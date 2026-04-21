#!/usr/bin/env nextflow
/*
========================================================================================
    Empirical Power Analysis — Main Workflow
    https://github.com/neurogenomics/Empirical_Power_Analysis
----------------------------------------------------------------------------------------
    Reproduces every figure, table, and source-data file in Crampton et al.
    "Current Single-Cell Study Designs Are Systematically Underpowered for
    Small-Effect Genes: Empirical Power Analysis Across 22 Million Detection Outcomes"

    Stages form a DAG rooted at sharding. Run with:
        nextflow run main.nf -profile slurm,apptainer
        nextflow run main.nf -profile test,docker          # toy pipeline
        nextflow run main.nf -profile standard -resume     # single-node + caching

    See docs/reproducing-the-paper.md for full instructions.
========================================================================================
*/

nextflow.enable.dsl = 2

// ======================================================================================
//  Parameters
// ======================================================================================
params {
    // Paths
    raw_data       = "${projectDir}/Datalink/raw_data"
    outdir         = "${projectDir}/Datalink"
    shards_dir     = "${params.outdir}/Shards"

    // Compute
    cores_default  = 16
    cores_heavy    = 64

    // Analysis parameters (override via config or CLI)
    downsample_mode = "cells"   // "cells" or "reads"
    n_robustness    = 6
    n_subsamples    = 20

    // Rendering
    render_paper   = true

    // Reporting
    publish_mode   = "copy"      // "copy" for final runs, "symlink" during development
}

// ======================================================================================
//  Module imports
// ======================================================================================
include { SHARDING }                from './modules/local/sharding'
include { EDA_ROUSSOS }             from './modules/local/eda_roussos'
include { FULL_GROUND_TRUTH }       from './modules/local/full_ground_truth'
include { DE_DIAGNOSTICS }          from './modules/local/de_diagnostics'
include { CARTESIAN_DOWNSAMPLING }  from './modules/local/cartesian_downsampling'
include { DEPTH_EQUALIZATION }      from './modules/local/depth_equalization'
include { POWER_CURVE_FITTING }     from './modules/local/power_curve_fitting'
include { BULK_DEG_PROCESSING }     from './modules/local/bulk_deg_processing'
include { FIGURE_DATA_GENERATION }  from './modules/local/figure_data_generation'
include { RENDER_PAPER }            from './modules/local/render_paper'

// ======================================================================================
//  Main workflow — encodes the README DAG
// ======================================================================================
workflow {

    // Header banner
    log.info """
    ================================================================================
      Empirical Power Analysis v${workflow.manifest.version ?: 'dev'}
      Profile:      ${workflow.profile}
      Project dir:  ${projectDir}
      Output dir:   ${params.outdir}
      Workflow run: ${workflow.runName}
    ================================================================================
    """.stripIndent()

    // Stage 1 — donor-level shards from raw h5ad/qs inputs
    raw_data_ch = Channel.fromPath(params.raw_data, checkIfExists: true)
    SHARDING(raw_data_ch)

    // Stage 2 — exploratory analyses (PCA, read-depth, filter rationale)
    //           runs in parallel with Stage 3 once sharding completes
    EDA_ROUSSOS(SHARDING.out.shards)

    // Stage 3 — ground truth DE across all datasets + cross-dataset diagnostics
    FULL_GROUND_TRUTH(SHARDING.out.shards)

    // Stage 4 — bulk GTEx DEG processing (independent of sharding)
    BULK_DEG_PROCESSING()

    // Stage 5 — Cartesian downsampling (the hot loop: 120k DESeq2 runs)
    CARTESIAN_DOWNSAMPLING(
        SHARDING.out.shards,
        FULL_GROUND_TRUTH.out.ground_truth
    )

    // Stage 6 — depth-equalised sensitivity analysis (astrocyte only)
    DEPTH_EQUALIZATION(
        SHARDING.out.shards,
        FULL_GROUND_TRUTH.out.ground_truth
    )

    // Stage 7 — power curve fitting (22M obs logistic regression)
    POWER_CURVE_FITTING(
        CARTESIAN_DOWNSAMPLING.out.results,
        DEPTH_EQUALIZATION.out.results
    )

    // Stage 8 — DE diagnostics (chromosome grid, shared-gene)
    DE_DIAGNOSTICS(FULL_GROUND_TRUTH.out.ground_truth)

    // Stage 9 — canonical CSVs consumed by the paper
    FIGURE_DATA_GENERATION(
        FULL_GROUND_TRUTH.out.ground_truth,
        CARTESIAN_DOWNSAMPLING.out.results,
        POWER_CURVE_FITTING.out.coefficients,
        BULK_DEG_PROCESSING.out.bulk_results
    )

    // Stage 10 — render the paper (optional)
    if (params.render_paper) {
        RENDER_PAPER(FIGURE_DATA_GENERATION.out.source_data)
    }
}

// ======================================================================================
//  Completion handler
// ======================================================================================
workflow.onComplete {
    log.info """
    ================================================================================
      Pipeline complete: ${workflow.success ? '✓ success' : '✗ failed'}
      Duration:  ${workflow.duration}
      Work dir:  ${workflow.workDir}
      Results:   ${params.outdir}
      Report:    ${params.outdir}/pipeline_info/execution_report.html
      DAG:       ${params.outdir}/pipeline_info/dag.svg
    ================================================================================
    """.stripIndent()

    if (!workflow.success) {
        log.info "  To resume from last checkpoint: nextflow run main.nf -profile ${workflow.profile} -resume"
    }
}

workflow.onError {
    log.error "Pipeline execution failed: ${workflow.errorMessage}"
}

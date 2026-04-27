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
//  Parameters — declared in nextflow.config (Nextflow 25.10+ requires it there;
//  script-level params{} blocks are no longer recognised by the strict parser)
// ======================================================================================

// ======================================================================================
//  Pre-flight validation — fail fast with a useful error before any work starts
// ======================================================================================
def preflight() {
    def raw = file(params.raw_data)
    if (!raw.exists()) {
        exit 1, """
        ERROR: raw_data path does not exist:
          ${params.raw_data}

        Expected a directory containing the raw h5ad / qs / tar inputs.
        See REPRODUCE.md for download instructions, or run with
          -profile test   (skips real inputs via stub blocks)
        """.stripIndent()
    }
    // For a real run, the directory must be non-empty (not just .gitkeep)
    if (raw.isDirectory() && workflow.profile?.tokenize(',')?.every { it != 'test' }) {
        def files = raw.list().findAll { !it.startsWith('.') }
        if (files.isEmpty()) {
            exit 1, """
            ERROR: raw_data directory is empty:
              ${params.raw_data}

            Populate it with the dataset tarballs listed in REPRODUCE.md,
            or use -profile test for a stub-only dry run.
            """.stripIndent()
        }
    }
}

// True only if the user asked to render AND the source file actually exists.
def shouldRender() {
    if (!params.render_paper) return false
    def rmd = file("${projectDir}/paper_final_draft.rmd")
    if (!rmd.exists()) {
        log.warn "render_paper=true but ${rmd} is missing — skipping rendering"
        return false
    }
    return true
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

    // Validate inputs before any process runs
    preflight()

    // Stage 1 — donor-level shards from raw h5ad/qs inputs.
    // Channel emits the raw_data directory itself; the sharding R scripts glob
    // inside it. type:'dir' + checkIfExists gives a clean error instead of a
    // cryptic "no files match" if the directory is missing.
    raw_data_ch = Channel.fromPath(params.raw_data, type: 'dir', checkIfExists: true)
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

    // Stage 10 — render the paper (optional; gated on file existence)
    if (shouldRender()) {
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

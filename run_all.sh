#!/usr/bin/env bash
# =============================================================================
# run_all.sh
#
# One-shot driver for the Empirical Power Analysis pipeline. Replaces Nextflow
# for the case where you just want to leave a single shell script churning for
# a few days and come back to a complete set of power curves.
#
# Stages run STRICTLY SERIALLY in DAG order:
#
#     1. sharding              (5 dataset shard scripts in sequence)
#     2. eda_roussos
#     3. full_ground_truth
#     4. bulk_deg_processing
#     5. cartesian_downsampling   ("cells" then "reads" if mode=both)
#     6. depth_equalization
#     7. power_curve_fitting
#     8. de_diagnostics
#     9. figure_data_generation
#
# Per-stage logs and a SUCCESS/FAILED status file are written under
# Datalink/Run_Output_Logs/run_all/<timestamp>/. With --no-overwrite, a stage
# is skipped if its previous run completed SUCCESS *and* the per-pipeline
# FORCE_RECALC check inside R also short-circuits work.
#
# Usage:
#   ./run_all.sh                       # use defaults from run_all.config
#   ./run_all.sh --max-cores 16
#   ./run_all.sh --no-overwrite        # skip work where outputs already exist
#   ./run_all.sh --force               # regenerate everything
#   ./run_all.sh --mode cells          # only the cells cartesian sweep
#   ./run_all.sh --mode both           # cells THEN reads (default)
#   ./run_all.sh --only full_ground_truth,cartesian_downsampling
#   ./run_all.sh --skip eda_roussos,de_diagnostics
#   ./run_all.sh --from cartesian_downsampling
#   ./run_all.sh --stop-on-fail
#   ./run_all.sh --dry-run
#
# Recommended invocation for an unattended multi-day run:
#   nohup ./run_all.sh --max-cores 24 --no-overwrite \
#       > nohup.run_all.out 2>&1 &
#   disown
# =============================================================================

set -u
set -o pipefail
# NB: NOT using `set -e` — we handle stage failures explicitly so the run
# continues across independent downstream stages unless --stop-on-fail is set.

# -----------------------------------------------------------------------------
# Locate ourselves
# -----------------------------------------------------------------------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="${SCRIPT_DIR}"
CONFIG_FILE="${PROJECT_ROOT}/run_all.config"

# -----------------------------------------------------------------------------
# Defaults (overridable by run_all.config and CLI flags)
# -----------------------------------------------------------------------------
MAX_CORES=24
FORCE_RECALC=false
CARTESIAN_MODE=both
SKIP_STAGES=""
ONLY_STAGES=""
FROM_STAGE=""
STOP_ON_FAIL=false
LOG_BASE="Datalink/Run_Output_Logs/run_all"
DRY_RUN=false
# Comma-separated list of dataset keys to shard (subset of: roussos,allen,dejager,johansen,rexach).
# The PCF wrapper sets this to "roussos" so unrelated datasets aren't re-sharded.
SHARDING_DATASETS="${SHARDING_DATASETS:-roussos,allen,dejager,johansen,rexach}"

# Load config file if present (override defaults)
if [[ -f "${CONFIG_FILE}" ]]; then
    # shellcheck source=/dev/null
    source "${CONFIG_FILE}"
fi

# -----------------------------------------------------------------------------
# CLI parsing (override config file)
# -----------------------------------------------------------------------------
print_usage() {
    sed -n '2,40p' "$0"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --max-cores)        MAX_CORES="$2"; shift 2 ;;
        --no-overwrite)     FORCE_RECALC=false; shift ;;
        --force)            FORCE_RECALC=true; shift ;;
        --mode)             CARTESIAN_MODE="$2"; shift 2 ;;
        --only)             ONLY_STAGES="$2"; shift 2 ;;
        --skip)             SKIP_STAGES="$2"; shift 2 ;;
        --from)             FROM_STAGE="$2"; shift 2 ;;
        --sharding-datasets) SHARDING_DATASETS="$2"; shift 2 ;;
        --stop-on-fail)     STOP_ON_FAIL=true; shift ;;
        --dry-run)          DRY_RUN=true; shift ;;
        -h|--help)          print_usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; print_usage; exit 2 ;;
    esac
done

# Validate
case "${CARTESIAN_MODE}" in
    cells|reads|both) ;;
    *) echo "ERROR: --mode must be cells|reads|both (got: ${CARTESIAN_MODE})" >&2; exit 2 ;;
esac

if ! [[ "${MAX_CORES}" =~ ^[0-9]+$ ]] || (( MAX_CORES < 1 )); then
    echo "ERROR: --max-cores must be a positive integer (got: ${MAX_CORES})" >&2
    exit 2
fi

# -----------------------------------------------------------------------------
# Lock — prevent two concurrent invocations from clobbering each other
# -----------------------------------------------------------------------------
LOCK_FILE="${PROJECT_ROOT}/.run_all.lock"
if [[ -e "${LOCK_FILE}" ]]; then
    EXISTING_PID=$(cat "${LOCK_FILE}" 2>/dev/null || echo "")
    if [[ -n "${EXISTING_PID}" ]] && kill -0 "${EXISTING_PID}" 2>/dev/null; then
        echo "ERROR: run_all.sh already running as PID ${EXISTING_PID}" >&2
        echo "       (lock file: ${LOCK_FILE})" >&2
        echo "       If you are sure no run is in progress, remove the lock file." >&2
        exit 1
    else
        echo "WARN: stale lock file from PID ${EXISTING_PID:-unknown} — removing"
        rm -f "${LOCK_FILE}"
    fi
fi
echo "$$" > "${LOCK_FILE}"
trap 'rm -f "${LOCK_FILE}"' EXIT

# -----------------------------------------------------------------------------
# Set up log directory + manifest
# -----------------------------------------------------------------------------
TIMESTAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_DIR_REL="${LOG_BASE}/${TIMESTAMP}"
LOG_DIR="${PROJECT_ROOT}/${LOG_DIR_REL}"
mkdir -p "${LOG_DIR}"
MASTER_LOG="${LOG_DIR}/00_run.log"

# Latest-pointer for easy access ("logs/run_all/latest -> <timestamp>")
LATEST_LINK="${PROJECT_ROOT}/${LOG_BASE}/latest"
mkdir -p "$(dirname "${LATEST_LINK}")"
ln -sfn "${TIMESTAMP}" "${LATEST_LINK}" || true

# Manifest captures the full invocation context so future-you can reproduce
GIT_SHA="$(git -C "${PROJECT_ROOT}" rev-parse --short HEAD 2>/dev/null || echo "unknown")"
GIT_DIRTY="$(git -C "${PROJECT_ROOT}" status --porcelain 2>/dev/null | wc -l || echo "0")"
cat > "${LOG_DIR}/manifest.json" <<EOF
{
  "timestamp": "${TIMESTAMP}",
  "host": "$(hostname)",
  "user": "$(whoami)",
  "project_root": "${PROJECT_ROOT}",
  "git_sha": "${GIT_SHA}",
  "git_dirty_files": ${GIT_DIRTY},
  "args": {
    "max_cores": ${MAX_CORES},
    "force_recalc": ${FORCE_RECALC},
    "cartesian_mode": "${CARTESIAN_MODE}",
    "only_stages": "${ONLY_STAGES}",
    "skip_stages": "${SKIP_STAGES}",
    "from_stage": "${FROM_STAGE}",
    "stop_on_fail": ${STOP_ON_FAIL},
    "dry_run": ${DRY_RUN},
    "sharding_datasets": "${SHARDING_DATASETS}"
  }
}
EOF

# -----------------------------------------------------------------------------
# Logging helpers
# -----------------------------------------------------------------------------
log() {
    local msg
    msg="[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"
    echo "${msg}" | tee -a "${MASTER_LOG}"
}

banner() {
    local title="$1"
    {
        echo ""
        echo "================================================================================"
        echo "  ${title}"
        echo "================================================================================"
    } | tee -a "${MASTER_LOG}"
}

# -----------------------------------------------------------------------------
# Environment for child Rscript processes
# -----------------------------------------------------------------------------
export POWER_ANALYSIS_ROOT="${PROJECT_ROOT}"
export POWER_ANALYSIS_FORCE_RECALC="${FORCE_RECALC}"
# Make the R process well-behaved for long unattended runs
export R_MAX_NUM_DLLS=200
# Avoid auto-loading user .Rprofile/.Renviron that might change behaviour
# (renv is honoured because it lives in the project root .Rprofile).

# -----------------------------------------------------------------------------
# Stage definitions
# -----------------------------------------------------------------------------
# Every stage is a function. They are invoked from the STAGE_ORDER array.
# A stage must:
#   - return 0 on success, non-zero on failure
#   - write its own per-stage log under ${LOG_DIR}/<stage>.log
#   - write ${LOG_DIR}/<stage>.status containing SUCCESS|FAILED|SKIPPED
# -----------------------------------------------------------------------------

STAGE_ORDER=(
    sharding
    eda_roussos
    full_ground_truth
    bulk_deg_processing
    cartesian_downsampling
    depth_equalization
    power_curve_fitting
    de_diagnostics
    figure_data_generation
)

# Run an Rscript with logging and a summary line.
# Args:  <stage_label>  <log_file>  <Rscript args...>
run_rscript() {
    local label="$1"; shift
    local logfile="$1"; shift
    local rc

    log "    -> Rscript $* (cores=${MAX_CORES}, FORCE_RECALC=${FORCE_RECALC})"
    if [[ "${DRY_RUN}" == "true" ]]; then
        echo "[DRY-RUN] would run: Rscript $*" >> "${logfile}"
        return 0
    fi

    # Run with a clean cwd at PROJECT_ROOT so here::here() resolves correctly
    # for the public-release configs.
    (
        cd "${PROJECT_ROOT}"
        Rscript "$@" >> "${logfile}" 2>&1
    )
    rc=$?
    if (( rc == 0 )); then
        log "    [${label}] OK"
    else
        log "    [${label}] FAILED (exit=${rc})"
    fi
    return ${rc}
}

# --- Stage: sharding (no run_pipeline.R; runs N shard scripts in sequence) ---
# Honours SHARDING_DATASETS (comma-separated keys: roussos,allen,dejager,johansen,rexach).
# Anything not in that list is skipped silently. Invalid keys are an error.
stage_sharding() {
    local logfile="${LOG_DIR}/sharding.log"
    : > "${logfile}"
    declare -A script_for=(
        [roussos]=R/sharding/Run_Shard_Roussos.R
        [allen]=R/sharding/Run_Shard_Allen.R
        [dejager]=R/sharding/Run_Shard_DeJager.R
        [johansen]=R/sharding/Run_Shard_Johansen.R
        [rexach]=R/sharding/Run_Shard_Rexach.R
    )
    local order=(roussos allen dejager johansen rexach)

    # Validate user's selection against the known keys
    IFS=',' read -ra requested <<< "${SHARDING_DATASETS}"
    for k in "${requested[@]}"; do
        if [[ -z "${script_for[$k]:-}" ]]; then
            log "    [sharding] ERROR: unknown dataset key '${k}' (allowed: ${!script_for[*]})"
            return 2
        fi
    done

    log "    [sharding] datasets: ${SHARDING_DATASETS}"
    local failed=0
    for k in "${order[@]}"; do
        if [[ ",${SHARDING_DATASETS}," != *",${k},"* ]]; then
            log "    [sharding] skip ${k} (not in SHARDING_DATASETS)"
            continue
        fi
        local s="${script_for[$k]}"
        run_rscript "sharding:${k}" "${logfile}" "${s}" || failed=1
        if (( failed == 1 )) && [[ "${STOP_ON_FAIL}" == "true" ]]; then
            return 1
        fi
    done
    return ${failed}
}

stage_eda_roussos()           { run_rscript "eda_roussos"           "${LOG_DIR}/eda_roussos.log"           R/eda_roussos/run_pipeline.R           all "${MAX_CORES}"; }
stage_full_ground_truth()     { run_rscript "full_ground_truth"     "${LOG_DIR}/full_ground_truth.log"     R/full_ground_truth/run_pipeline.R     all "${MAX_CORES}"; }
stage_bulk_deg_processing()   { run_rscript "bulk_deg_processing"   "${LOG_DIR}/bulk_deg_processing.log"   R/bulk_deg_processing/run_pipeline.R   all "${MAX_CORES}"; }
stage_depth_equalization()    { run_rscript "depth_equalization"    "${LOG_DIR}/depth_equalization.log"    R/depth_equalization/run_pipeline.R    all "${MAX_CORES}"; }
stage_power_curve_fitting()   { run_rscript "power_curve_fitting"   "${LOG_DIR}/power_curve_fitting.log"   R/power_curve_fitting/run_pipeline.R   all "${MAX_CORES}"; }
stage_de_diagnostics()        { run_rscript "de_diagnostics"        "${LOG_DIR}/de_diagnostics.log"        R/de_diagnostics/run_pipeline.R        all "${MAX_CORES}"; }
stage_figure_data_generation(){ run_rscript "figure_data_generation" "${LOG_DIR}/figure_data_generation.log" R/figure_data_generation/run_pipeline.R all "${MAX_CORES}"; }

# Cartesian respects --mode (cells|reads|both). For "both" we run twice in the
# same stage so the SUCCESS/FAILED bookkeeping covers both modes.
stage_cartesian_downsampling() {
    local logfile="${LOG_DIR}/cartesian_downsampling.log"
    : > "${logfile}"
    local modes
    case "${CARTESIAN_MODE}" in
        cells) modes=(cells) ;;
        reads) modes=(reads) ;;
        both)  modes=(cells reads) ;;
    esac
    local m rc=0
    for m in "${modes[@]}"; do
        log "  cartesian_downsampling mode=${m}"
        DOWNSAMPLE_MODE="${m}" run_rscript "cartesian:${m}" "${logfile}" \
            R/cartesian_downsampling/run_pipeline.R all "${MAX_CORES}" || rc=$?
        if (( rc != 0 )) && [[ "${STOP_ON_FAIL}" == "true" ]]; then
            return ${rc}
        fi
    done
    return ${rc}
}

# -----------------------------------------------------------------------------
# Stage filter logic — decide which stages to actually run
# -----------------------------------------------------------------------------
should_run_stage() {
    local stage="$1"

    # ONLY_STAGES wins over everything else if set
    if [[ -n "${ONLY_STAGES}" ]]; then
        [[ ",${ONLY_STAGES}," == *",${stage},"* ]] && return 0 || return 1
    fi

    # Check explicit skip list
    if [[ -n "${SKIP_STAGES}" ]] && [[ ",${SKIP_STAGES}," == *",${stage},"* ]]; then
        return 1
    fi

    # FROM_STAGE: skip everything before it
    if [[ -n "${FROM_STAGE}" ]]; then
        if [[ -z "${REACHED_FROM:-}" ]]; then
            if [[ "${stage}" == "${FROM_STAGE}" ]]; then
                REACHED_FROM=1
            else
                return 1
            fi
        fi
    fi

    return 0
}

# -----------------------------------------------------------------------------
# Banner
# -----------------------------------------------------------------------------
banner "EMPIRICAL POWER ANALYSIS — run_all.sh"
log "Project root:    ${PROJECT_ROOT}"
log "Log dir:         ${LOG_DIR_REL}"
log "Git SHA:         ${GIT_SHA} ($(if (( GIT_DIRTY > 0 )); then echo "${GIT_DIRTY} dirty file(s)"; else echo clean; fi))"
log "Max cores:       ${MAX_CORES}"
log "Force recalc:    ${FORCE_RECALC}"
log "Cartesian mode:  ${CARTESIAN_MODE}"
log "Stop on fail:    ${STOP_ON_FAIL}"
log "Only stages:     ${ONLY_STAGES:-<all>}"
log "Skip stages:     ${SKIP_STAGES:-<none>}"
log "From stage:      ${FROM_STAGE:-<start>}"
log "Dry run:         ${DRY_RUN}"

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------
overall_start=$SECONDS
declare -A stage_status
declare -A stage_seconds
n_run=0
n_ok=0
n_fail=0
n_skipped=0
exit_code=0

for stage in "${STAGE_ORDER[@]}"; do
    if ! should_run_stage "${stage}"; then
        banner "STAGE: ${stage}  (SKIPPED)"
        echo "SKIPPED" > "${LOG_DIR}/${stage}.status"
        stage_status[${stage}]="SKIPPED"
        ((n_skipped++))
        continue
    fi

    banner "STAGE: ${stage}"
    stage_start=$SECONDS

    # Dispatch to the matching stage_<name> function
    fn="stage_${stage}"
    if ! declare -F "${fn}" >/dev/null; then
        log "ERROR: no implementation for stage '${stage}' (function ${fn} not defined)"
        echo "FAILED" > "${LOG_DIR}/${stage}.status"
        stage_status[${stage}]="FAILED"
        ((n_fail++))
        ((n_run++))
        if [[ "${STOP_ON_FAIL}" == "true" ]]; then exit_code=1; break; fi
        continue
    fi

    if "${fn}"; then
        echo "SUCCESS" > "${LOG_DIR}/${stage}.status"
        stage_status[${stage}]="SUCCESS"
        ((n_ok++))
    else
        echo "FAILED" > "${LOG_DIR}/${stage}.status"
        stage_status[${stage}]="FAILED"
        ((n_fail++))
        if [[ "${STOP_ON_FAIL}" == "true" ]]; then
            log "STOP_ON_FAIL=true — aborting after ${stage} failure"
            exit_code=1
            ((n_run++))
            break
        fi
    fi
    stage_seconds[${stage}]=$((SECONDS - stage_start))
    ((n_run++))
done

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
overall_seconds=$((SECONDS - overall_start))
banner "RUN COMPLETE"
log "Total wall time: $((overall_seconds / 3600))h $(((overall_seconds % 3600) / 60))m $((overall_seconds % 60))s"
log "Stages run:      ${n_run}   ok=${n_ok}  failed=${n_fail}  skipped=${n_skipped}"
log ""
log "Per-stage status:"
for stage in "${STAGE_ORDER[@]}"; do
    s="${stage_status[${stage}]:-NOT-REACHED}"
    secs="${stage_seconds[${stage}]:-0}"
    printf "  %-26s %-12s  %dm %ds\n" "${stage}" "${s}" $((secs / 60)) $((secs % 60)) | tee -a "${MASTER_LOG}"
done
log ""
log "Logs:    ${LOG_DIR_REL}"
log "Latest:  ${LOG_BASE}/latest -> ${TIMESTAMP}"

if (( n_fail > 0 )) && [[ "${exit_code}" == "0" ]]; then
    exit_code=2  # non-fatal failures, but distinguish from clean success
fi

exit ${exit_code}

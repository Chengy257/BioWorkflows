#!/usr/bin/env bash
# seclip-seq workflow launcher.
# Provides backward-compatible positional arguments plus a safer option-based CLI.

set -Eeuo pipefail

SCRIPT_VERSION="0.1.0"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$REPO_DIR/workflow"
SNAKEFILE="$WORKFLOW_DIR/Snakefile"
DEFAULT_CONFIG="$REPO_DIR/config/config.yaml"
DEFAULT_SOFTWARE="$REPO_DIR/config/software.yaml"
RUNTIME_HELPER="$WORKFLOW_DIR/scripts/runtime_config.py"
# shellcheck disable=SC2034  # used by the sourced shared launcher (resolve_path)
CALL_DIR="$PWD"

# Shared cross-project launcher helpers (logging + path resolution).
BIO_WORKFLOWS_SHARED="${BIO_WORKFLOWS_SHARED:-$(cd "$REPO_DIR/.." && pwd)/shared}"
export BIO_WORKFLOWS_SHARED
source "$BIO_WORKFLOWS_SHARED/lib/launcher.sh"

PROJECT_DIR=""
CONFIG_PATH=""
SOFTWARE_PATH="${SECLIP_SOFTWARE_CONFIG:-}"
EXTRA_CONFIG=""
JOBS="${SECLIP_JOBS:-10}"
PROFILE_REQUEST="${RUN_PROFILE:-auto}"
LOG_PATH="${SECLIP_LOG:-snakemake.logs.txt}"
DRY_RUN=false
VALIDATE_ONLY=false
CHECK_SOFTWARE_ONLY=false
SKIP_SOFTWARE_CHECK=false
UNLOCK=false
QUIET=false
QUEUE="${SECLIP_QUEUE:-}"
PARTITION="${SECLIP_PARTITION:-}"
MEMORY="${SECLIP_MEMORY:-}"
RUNTIME_MIN="${SECLIP_RUNTIME_MIN:-}"
SGE_MEMORY_RESOURCE="${SECLIP_SGE_MEMORY_RESOURCE:-h_vmem}"
SCHEDULER_EXTRA="${SECLIP_SCHEDULER_EXTRA:-}"
RETRIES="${SECLIP_RETRIES:-0}"
LATENCY_WAIT="${SECLIP_LATENCY_WAIT:-}"
MAX_JOBS_PER_SECOND="${SECLIP_MAX_JOBS_PER_SECOND:-}"
MAX_STATUS_CHECKS_PER_SECOND="${SECLIP_MAX_STATUS_CHECKS_PER_SECOND:-}"

POSITIONAL=()
SNAKEMAKE_EXTRA_ARGS=()




usage() {
    cat <<'EOF'
seclip-seq workflow launcher

Usage:
  bash run.sh <project_dir> [config.yaml] [jobs]
  bash run.sh [options]

Single pipeline:
  QC, UMI extraction, trimming, repeats filter, unique alignment,
  UMI deduplication, peak calling (PureCLIP / CLIPper)

Core options:
  -P, --project DIR         Project directory; created if it does not exist
  -c, --config FILE         Project analysis configuration file
      --software FILE       Software/runtime configuration (software.yaml)
  -l, --extra-config FILE   Extra config layered last; config.local.yaml in the
                            project directory is picked up automatically
  -j, --jobs N              Maximum parallel jobs or local cores (default: 10)
      --profile NAME        auto, default, local, pbs, sge, or slurm
      --queue NAME          SGE queue name (SGE only)
      --partition NAME      SLURM partition name (SLURM only)
      --memory VALUE        Override per-rule memory for all cluster jobs, e.g. 16G
      --runtime MIN         Override per-rule walltime for all cluster jobs (minutes)
      --sge-mem-resource N  SGE memory resource name (default: h_vmem)
      --scheduler-extra S   Extra text appended to qsub/sbatch submit command
      --retries N           Retry failed jobs up to N times (default: 0)
      --latency-wait SEC    Override profile filesystem latency wait
      --max-jobs-per-sec N  Limit job submission rate
      --max-status-per-sec N Limit scheduler status checks
  -n, --dry-run             Build the DAG and show planned jobs without executing
      --validate-only       Run sample and software/runtime validation only
      --check-software      Check required executables and databases; then exit
      --skip-software-check Skip runtime preflight before a real workflow run
      --unlock              Remove a stale Snakemake working-directory lock
      --log FILE            Launcher/Snakemake log path (default: snakemake.logs.txt)
  -q, --quiet               Reduce launcher output; Snakemake output is unchanged
  -h, --help                Show this help message
      --version             Show launcher version

Scheduler selection:
  auto selects SGE when qsub is available, then SLURM when sbatch is available,
  otherwise the local default profile. RUN_PROFILE can set the same preference.

Any arguments after "--" are passed directly to Snakemake.
EOF
    cat <<'EOF'

Examples:
  bash run.sh /data/project /data/project/config.yaml 20
  bash run.sh --project /data/project --config config.yaml --jobs 20
  bash run.sh -P . -n
  bash run.sh -P . --software software.yaml --check-software
  RUN_PROFILE=slurm bash run.sh -P . -c config.yaml -j 40
  bash run.sh -P . --profile slurm --partition compute --memory 32G
  bash run.sh -P . --profile sge --queue all.q --runtime 720 --retries 2
  bash run.sh -P . -- --rerun-triggers mtime
EOF
}

print_version() {
    printf 'run.sh %s\n' "$SCRIPT_VERSION"
}


while (($#)); do
    case "$1" in
        --project=*|--project-dir=*) PROJECT_DIR="${1#*=}"; shift ;;
        --config=*) CONFIG_PATH="${1#*=}"; shift ;;
        --software=*) SOFTWARE_PATH="${1#*=}"; shift ;;
        --extra-config=*) EXTRA_CONFIG="${1#*=}"; shift ;;
        --jobs=*) JOBS="${1#*=}"; shift ;;
        --profile=*) PROFILE_REQUEST="${1#*=}"; shift ;;
        --queue=*) QUEUE="${1#*=}"; shift ;;
        --partition=*) PARTITION="${1#*=}"; shift ;;
        --memory=*|--mem=*) MEMORY="${1#*=}"; shift ;;
        --runtime=*) RUNTIME_MIN="${1#*=}"; shift ;;
        --sge-mem-resource=*) SGE_MEMORY_RESOURCE="${1#*=}"; shift ;;
        --scheduler-extra=*) SCHEDULER_EXTRA="${1#*=}"; shift ;;
        --retries=*) RETRIES="${1#*=}"; shift ;;
        --latency-wait=*) LATENCY_WAIT="${1#*=}"; shift ;;
        --max-jobs-per-sec=*) MAX_JOBS_PER_SECOND="${1#*=}"; shift ;;
        --max-status-per-sec=*) MAX_STATUS_CHECKS_PER_SECOND="${1#*=}"; shift ;;
        --log=*) LOG_PATH="${1#*=}"; shift ;;
        -P|--project|--project-dir) [[ $# -ge 2 ]] || die "$1 requires a value"; PROJECT_DIR="$2"; shift 2 ;;
        -c|--config) [[ $# -ge 2 ]] || die "$1 requires a value"; CONFIG_PATH="$2"; shift 2 ;;
        --software) [[ $# -ge 2 ]] || die "$1 requires a value"; SOFTWARE_PATH="$2"; shift 2 ;;
        -l|--extra-config) [[ $# -ge 2 ]] || die "$1 requires a value"; EXTRA_CONFIG="$2"; shift 2 ;;
        -j|--jobs) [[ $# -ge 2 ]] || die "$1 requires a value"; JOBS="$2"; shift 2 ;;
        --profile) [[ $# -ge 2 ]] || die "$1 requires a value"; PROFILE_REQUEST="$2"; shift 2 ;;
        --queue) [[ $# -ge 2 ]] || die "$1 requires a value"; QUEUE="$2"; shift 2 ;;
        --partition) [[ $# -ge 2 ]] || die "$1 requires a value"; PARTITION="$2"; shift 2 ;;
        --memory|--mem) [[ $# -ge 2 ]] || die "$1 requires a value"; MEMORY="$2"; shift 2 ;;
        --runtime) [[ $# -ge 2 ]] || die "$1 requires a value"; RUNTIME_MIN="$2"; shift 2 ;;
        --sge-mem-resource) [[ $# -ge 2 ]] || die "$1 requires a value"; SGE_MEMORY_RESOURCE="$2"; shift 2 ;;
        --scheduler-extra) [[ $# -ge 2 ]] || die "$1 requires a value"; SCHEDULER_EXTRA="$2"; shift 2 ;;
        --retries) [[ $# -ge 2 ]] || die "$1 requires a value"; RETRIES="$2"; shift 2 ;;
        --latency-wait) [[ $# -ge 2 ]] || die "$1 requires a value"; LATENCY_WAIT="$2"; shift 2 ;;
        --max-jobs-per-sec) [[ $# -ge 2 ]] || die "$1 requires a value"; MAX_JOBS_PER_SECOND="$2"; shift 2 ;;
        --max-status-per-sec) [[ $# -ge 2 ]] || die "$1 requires a value"; MAX_STATUS_CHECKS_PER_SECOND="$2"; shift 2 ;;
        -n|--dry-run) DRY_RUN=true; shift ;;
        --validate-only) VALIDATE_ONLY=true; shift ;;
        --check-software) CHECK_SOFTWARE_ONLY=true; shift ;;
        --skip-software-check) SKIP_SOFTWARE_CHECK=true; shift ;;
        --unlock) UNLOCK=true; shift ;;
        --log) [[ $# -ge 2 ]] || die "$1 requires a value"; LOG_PATH="$2"; shift 2 ;;
        -q|--quiet) QUIET=true; shift ;;
        -h|--help) usage; exit 0 ;;
        --version) print_version; exit 0 ;;
        --) shift; SNAKEMAKE_EXTRA_ARGS+=("$@"); break ;;
        -*) die "Unknown option: $1. Use --help for usage." ;;
        *) POSITIONAL+=("$1"); shift ;;
    esac
done

# Backward-compatible positional form:
#   run.sh <project_dir> [config.yaml] [jobs]
[[ -n "$PROJECT_DIR" ]] || PROJECT_DIR="${POSITIONAL[0]:-}"
[[ -n "$CONFIG_PATH" ]] || CONFIG_PATH="${POSITIONAL[1]:-}"
if [[ ${#POSITIONAL[@]} -ge 3 && "$JOBS" == "${SECLIP_JOBS:-10}" ]]; then
    JOBS="${POSITIONAL[2]}"
fi
[[ ${#POSITIONAL[@]} -le 3 ]] || die "Too many positional arguments. Use --help for usage."
[[ -n "$PROJECT_DIR" ]] || { usage >&2; die "Project directory is required."; }
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || die "--jobs must be a positive integer, got '$JOBS'."
[[ "$RETRIES" =~ ^[0-9]+$ ]] || die "--retries must be a non-negative integer, got '$RETRIES'."
[[ -z "$LATENCY_WAIT" || "$LATENCY_WAIT" =~ ^[0-9]+$ ]] || die "--latency-wait must be a non-negative integer."
[[ -z "$MAX_JOBS_PER_SECOND" || "$MAX_JOBS_PER_SECOND" =~ ^[0-9]+([.][0-9]+)?$ ]] || die "--max-jobs-per-sec must be numeric."
[[ -z "$MAX_STATUS_CHECKS_PER_SECOND" || "$MAX_STATUS_CHECKS_PER_SECOND" =~ ^[0-9]+([.][0-9]+)?$ ]] || die "--max-status-per-sec must be numeric."
[[ -z "$MEMORY" || "$MEMORY" =~ ^[1-9][0-9]*([KMGTP]i?[Bb]?|[KMGTP])?$ ]] || die "--memory must look like 8000, 8G, 16GB, or 32GiB."
[[ -z "$RUNTIME_MIN" || "$RUNTIME_MIN" =~ ^[1-9][0-9]*$ ]] || die "--runtime must be a positive integer number of minutes."
[[ "$SGE_MEMORY_RESOURCE" =~ ^[A-Za-z_][A-Za-z0-9_.-]*$ ]] || die "Invalid SGE memory resource name: $SGE_MEMORY_RESOURCE"

PROJECT_DIR="$(resolve_path "$PROJECT_DIR")"
mkdir -p "$PROJECT_DIR"
PROJECT_DIR="$(cd "$PROJECT_DIR" && pwd)"

if [[ -n "$CONFIG_PATH" ]]; then
    CONFIG_PATH="$(resolve_path "$CONFIG_PATH")"
elif [[ -f "$PROJECT_DIR/config.yaml" ]]; then
    CONFIG_PATH="$PROJECT_DIR/config.yaml"
fi

if [[ -n "$CONFIG_PATH" ]]; then
    [[ -f "$CONFIG_PATH" ]] || die "Configuration file not found: $CONFIG_PATH"
    CONFIG_PATH="$(cd "$(dirname "$CONFIG_PATH")" && pwd)/$(basename "$CONFIG_PATH")"
    CONFIG_FOR_PARSE="$CONFIG_PATH"
else
    [[ -f "$DEFAULT_CONFIG" ]] || die "Default configuration file not found: $DEFAULT_CONFIG"
    CONFIG_FOR_PARSE="$DEFAULT_CONFIG"
fi

if [[ -n "$SOFTWARE_PATH" ]]; then
    SOFTWARE_PATH="$(resolve_path "$SOFTWARE_PATH")"
elif [[ -f "$PROJECT_DIR/software.yaml" ]]; then
    SOFTWARE_PATH="$PROJECT_DIR/software.yaml"
else
    SOFTWARE_PATH="$DEFAULT_SOFTWARE"
fi
[[ -f "$SOFTWARE_PATH" ]] || die "Software configuration file not found: $SOFTWARE_PATH"
SOFTWARE_PATH="$(cd "$(dirname "$SOFTWARE_PATH")" && pwd)/$(basename "$SOFTWARE_PATH")"

# Per-rule scheduler resources: project-local resources.yaml wins, otherwise
# the repository default (config/resources.yaml). Injected as SECLIP_RESOURCES_CONFIG.
if [[ -n "${SECLIP_RESOURCES_CONFIG:-}" ]]; then
    RESOURCES_PATH="$(resolve_path "$SECLIP_RESOURCES_CONFIG")"
elif [[ -f "$PROJECT_DIR/resources.yaml" ]]; then
    RESOURCES_PATH="$PROJECT_DIR/resources.yaml"
else
    RESOURCES_PATH="$REPO_DIR/config/resources.yaml"
fi
[[ -f "$RESOURCES_PATH" ]] || die "Resources configuration file not found: $RESOURCES_PATH"
RESOURCES_PATH="$(cd "$(dirname "$RESOURCES_PATH")" && pwd)/$(basename "$RESOURCES_PATH")"
export SECLIP_RESOURCES_CONFIG="$RESOURCES_PATH"

# Extra config layer: explicit -l/--extra-config wins; otherwise pick up
# config.local.yaml from the project directory automatically.
if [[ -n "$EXTRA_CONFIG" ]]; then
    EXTRA_CONFIG="$(resolve_path "$EXTRA_CONFIG")"
    [[ -f "$EXTRA_CONFIG" ]] || die "Extra configuration file not found: $EXTRA_CONFIG"
    EXTRA_CONFIG="$(cd "$(dirname "$EXTRA_CONFIG")" && pwd)/$(basename "$EXTRA_CONFIG")"
elif [[ -f "$PROJECT_DIR/config.local.yaml" ]]; then
    EXTRA_CONFIG="$PROJECT_DIR/config.local.yaml"
    info "Found config.local.yaml in the project directory; layering it last"
fi
if [[ -n "$EXTRA_CONFIG" ]]; then
    export SECLIP_EXTRA_CONFIG="$EXTRA_CONFIG"
fi

[[ -f "$SNAKEFILE" ]] || die "Snakefile not found: $SNAKEFILE"
[[ -f "$RUNTIME_HELPER" ]] || die "Runtime helper not found: $RUNTIME_HELPER"
command -v python3 >/dev/null 2>&1 || die "python3 is required to resolve software.yaml."
RUNTIME_PYTHON="$(command -v python3)"
"$RUNTIME_PYTHON" -c 'import yaml' >/dev/null 2>&1 || die "PyYAML is required by the launcher runtime resolver."
RUNTIME_EXPORTS="$("$RUNTIME_PYTHON" "$RUNTIME_HELPER" export --config "$SOFTWARE_PATH")" || die "Unable to resolve software runtime: $SOFTWARE_PATH"
eval "$RUNTIME_EXPORTS"
export SECLIP_SOFTWARE_CONFIG="$SOFTWARE_PATH"

command -v snakemake >/dev/null 2>&1 || die "snakemake was not found after applying software.yaml. Add it to the main environment or PATH."

if command -v timeout >/dev/null 2>&1; then
    SNAKEMAKE_VERSION="$(timeout 8s snakemake --version 2>/dev/null | head -1 || true)"
else
    SNAKEMAKE_VERSION="$(snakemake --version 2>/dev/null | head -1 || true)"
fi
if [[ -z "$SNAKEMAKE_VERSION" ]]; then
    warn "Unable to determine the Snakemake version; continuing anyway."
fi

select_profile() {
    local requested
    requested="$(printf '%s' "$PROFILE_REQUEST" | tr '[:upper:]' '[:lower:]')"
    [[ "$requested" == "local" ]] && requested="default"
    if [[ "$requested" == "auto" ]]; then
        if command -v sbatch >/dev/null 2>&1; then
            requested="slurm"
        elif command -v qsub >/dev/null 2>&1; then
            # qsub is ambiguous: an SGE installation always sets SGE_ROOT; otherwise treat as PBS
            if [[ -n "${SGE_ROOT:-}" ]]; then requested="sge"; else requested="pbs"; fi
        else
            requested="default"
        fi
    fi

    case "$requested" in
        default|pbs|sge|slurm) ;;
        *) die "Unknown profile '$PROFILE_REQUEST'. Choose auto, default, pbs, sge, or slurm." ;;
    esac

    PROFILE="$requested"
    PROFILE_DIR="$WORKFLOW_DIR/profile/$PROFILE"
    [[ -d "$PROFILE_DIR" ]] || die "Profile directory not found: $PROFILE_DIR"

    case "$PROFILE" in
        pbs|sge) command -v qsub >/dev/null 2>&1 || die "$PROFILE profile selected but qsub was not found in PATH." ;;
        slurm) command -v sbatch >/dev/null 2>&1 || die "SLURM profile selected but sbatch was not found in PATH." ;;
    esac
}

runtime_check() {
    local scope="$1"
    "$RUNTIME_PYTHON" "$RUNTIME_HELPER" check \
        --config "$SOFTWARE_PATH" --scope "$scope" \
        --analysis-config "$CONFIG_FOR_PARSE"
}

if [[ "$CHECK_SOFTWARE_ONLY" == true ]]; then
    runtime_check all
    info "Software/runtime check completed successfully."
    exit 0
fi

select_profile

if [[ -n "$QUEUE" && "$PROFILE" != "pbs" && "$PROFILE" != "sge" ]]; then
    die "--queue is only valid with the PBS or SGE profile."
fi
if [[ -n "$PARTITION" && "$PROFILE" != "slurm" ]]; then
    die "--partition is only valid with the SLURM profile."
fi
if [[ -n "$MEMORY" && "$PROFILE" == "default" ]]; then
    die "--memory is a cluster scheduler request and is not valid with the local profile."
fi
if [[ -n "$RUNTIME_MIN" && "$PROFILE" == "default" ]]; then
    die "--runtime is a cluster scheduler request and is not valid with the local profile."
fi
if [[ -n "$SCHEDULER_EXTRA" && "$PROFILE" == "default" ]]; then
    die "--scheduler-extra is only valid with SGE or SLURM."
fi

CLUSTER_CMD=""
case "$PROFILE" in
    sge)
        if [[ -n "$MEMORY" ]]; then
            SGE_MEM_REQUEST="$MEMORY"
        else
            SGE_MEM_REQUEST='{resources.mem_mb}M'
        fi
        if [[ -n "$RUNTIME_MIN" ]]; then
            SGE_RUNTIME_REQUEST="$((RUNTIME_MIN * 60))"
        else
            SGE_RUNTIME_REQUEST='{resources.runtime_sec}'
        fi
        CLUSTER_CMD="qsub -V -N {rule} -l ncpus={threads} -l ${SGE_MEMORY_RESOURCE}=${SGE_MEM_REQUEST} -l h_rt=${SGE_RUNTIME_REQUEST} -j oe"
        [[ -n "$QUEUE" ]] && CLUSTER_CMD+=" -q $QUEUE"
        ;;
    slurm)
        if [[ -n "$MEMORY" ]]; then
            SLURM_MEM_REQUEST="$MEMORY"
        else
            SLURM_MEM_REQUEST='{resources.mem_mb}M'
        fi
        if [[ -n "$RUNTIME_MIN" ]]; then
            SLURM_RUNTIME_REQUEST="$RUNTIME_MIN"
        else
            SLURM_RUNTIME_REQUEST='{resources.runtime_min}'
        fi
        CLUSTER_CMD="sbatch --parsable -J {rule} -c {threads} --mem=${SLURM_MEM_REQUEST} --time=${SLURM_RUNTIME_REQUEST} -o slurm-{rule}-%j.out"
        [[ -n "$PARTITION" ]] && CLUSTER_CMD+=" -p $PARTITION"
        ;;
esac
[[ -n "$CLUSTER_CMD" && -n "$SCHEDULER_EXTRA" ]] && CLUSTER_CMD+=" $SCHEDULER_EXTRA"

if [[ -n "$CONFIG_PATH" ]]; then
    export SECLIP_CONFIG="$CONFIG_PATH"
else
    unset SECLIP_CONFIG || true
fi

cd "$PROJECT_DIR"
if [[ "$SKIP_SOFTWARE_CHECK" != true && "$DRY_RUN" != true ]]; then
    runtime_check all
elif [[ "$DRY_RUN" == true && "$SKIP_SOFTWARE_CHECK" != true ]]; then
    info "Software/runtime preflight is skipped for dry-run; use --check-software for an explicit check."
fi

if [[ ! -d "$PROJECT_DIR/1.rawdata" ]]; then
    warn "Raw-data directory is missing: $PROJECT_DIR/1.rawdata. This is acceptable only when all required upstream outputs already exist."
fi

if [[ "$QUIET" != true ]]; then
    info "seclip-seq launcher version: $SCRIPT_VERSION"
    info "Project: $PROJECT_DIR"
    info "Config: ${CONFIG_PATH:-$DEFAULT_CONFIG}"
    info "Software config: $SOFTWARE_PATH"
    info "Resources config: $RESOURCES_PATH"
    info "Software environment: ${SECLIP_SOFTWARE_TYPE:-system}${SECLIP_ENV_PREFIX:+ ($SECLIP_ENV_PREFIX)}"
    [[ -n "${R_LIBS_USER:-}" ]] && info "R libraries: $R_LIBS_USER"
    info "Profile: $PROFILE ($PROFILE_DIR)"
    info "Jobs: $JOBS"
    [[ "$RETRIES" != "0" ]] && info "Retries: $RETRIES"
    [[ -n "$LATENCY_WAIT" ]] && info "Latency wait: ${LATENCY_WAIT}s"
    [[ -n "$QUEUE" ]] && info "SGE queue: $QUEUE"
    [[ -n "$PARTITION" ]] && info "SLURM partition: $PARTITION"
    if [[ -n "$MEMORY" ]]; then
        info "Scheduler memory override: $MEMORY"
    elif [[ "$PROFILE" != "default" ]]; then
        info "Scheduler memory: per-rule resources.mem_mb"
    fi
    if [[ -n "$RUNTIME_MIN" ]]; then
        info "Scheduler runtime override: ${RUNTIME_MIN} min"
    elif [[ "$PROFILE" != "default" ]]; then
        info "Scheduler runtime: per-rule resources.runtime_min"
    fi
    [[ -n "$CLUSTER_CMD" ]] && info "Submit command: $CLUSTER_CMD"
    [[ -n "$SNAKEMAKE_VERSION" ]] && info "Snakemake: $SNAKEMAKE_VERSION"
fi
if [[ "$VALIDATE_ONLY" == true ]]; then
    info "Validation completed successfully."
    exit 0
fi

SNAKEMAKE_ARGS=(
    -s "$SNAKEFILE"
    --profile "$PROFILE_DIR"
    -j "$JOBS"
)

[[ -n "$CLUSTER_CMD" ]] && SNAKEMAKE_ARGS+=(--cluster "$CLUSTER_CMD")
[[ "$RETRIES" != "0" ]] && SNAKEMAKE_ARGS+=(--retries "$RETRIES")
[[ -n "$LATENCY_WAIT" ]] && SNAKEMAKE_ARGS+=(--latency-wait "$LATENCY_WAIT")
[[ -n "$MAX_JOBS_PER_SECOND" ]] && SNAKEMAKE_ARGS+=(--max-jobs-per-second "$MAX_JOBS_PER_SECOND")
[[ -n "$MAX_STATUS_CHECKS_PER_SECOND" ]] && SNAKEMAKE_ARGS+=(--max-status-checks-per-second "$MAX_STATUS_CHECKS_PER_SECOND")

if [[ "$DRY_RUN" == true ]]; then
    SNAKEMAKE_ARGS+=(--dry-run)
fi

if [[ "$UNLOCK" == true ]]; then
    info "Unlocking Snakemake working directory."
    snakemake "${SNAKEMAKE_ARGS[@]}" --unlock
    exit 0
fi

if ((${#SNAKEMAKE_EXTRA_ARGS[@]})); then
    SNAKEMAKE_ARGS+=("${SNAKEMAKE_EXTRA_ARGS[@]}")
fi

if [[ "$LOG_PATH" != /* ]]; then
    LOG_PATH="$PROJECT_DIR/$LOG_PATH"
fi
mkdir -p "$(dirname "$LOG_PATH")"

START_EPOCH="$(date +%s)"
on_exit() {
    local status=$?
    local end elapsed
    end="$(date +%s)"
    elapsed=$((end - START_EPOCH))
    if (( status == 0 )); then
        info "Workflow finished successfully in ${elapsed}s."
    else
        warn "Workflow exited with status $status after ${elapsed}s. Check: $LOG_PATH"
    fi
}
trap on_exit EXIT
trap 'warn "Interrupted by user."; exit 130' INT TERM

info "Log: $LOG_PATH"
if [[ "$DRY_RUN" == true ]]; then
    info "Dry-run mode enabled; no workflow jobs will be executed."
fi

set +e
snakemake "${SNAKEMAKE_ARGS[@]}" 2>&1 | tee "$LOG_PATH"
status=${PIPESTATUS[0]}
set -e
exit "$status"

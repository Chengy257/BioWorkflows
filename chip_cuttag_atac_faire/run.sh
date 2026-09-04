#!/usr/bin/env bash
# chip_cuttag_atac_faire workflow launcher (ChIP-seq / CUT&Tag / ATAC-seq / FAIRE-seq).
# Provides backward-compatible positional arguments plus a safer option-based CLI.

set -Eeuo pipefail

SCRIPT_VERSION="0.4.0"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$REPO_DIR/workflow"
SNAKEFILE="$WORKFLOW_DIR/Snakefile"
DEFAULT_CONFIG="$REPO_DIR/config/config.yaml"
DEFAULT_SOFTWARE="$REPO_DIR/config/software.yaml"
RUNTIME_HELPER="$WORKFLOW_DIR/scripts/runtime_config.py"
CALL_DIR="$PWD"

PROJECT_DIR=""
CONFIG_PATH=""
EXTRA_CONFIG=""
SOFTWARE_PATH="${CHIP_SOFTWARE_CONFIG:-}"
JOBS="${CHIP_JOBS:-3}"
PROFILE_REQUEST="${RUN_PROFILE:-auto}"
LOG_PATH="${CHIP_LOG:-snakemake.logs.txt}"
DRY_RUN=false
VALIDATE_ONLY=false
CHECK_SOFTWARE_ONLY=false
CHECK_R_ONLY=false
SKIP_SOFTWARE_CHECK=false
UNLOCK=false
QUIET=false
RENAME=false
QUEUE="${CHIP_QUEUE:-}"
PARTITION="${CHIP_PARTITION:-}"
MEMORY="${CHIP_MEMORY:-}"
RUNTIME_MIN="${CHIP_RUNTIME_MIN:-}"
SGE_MEMORY_RESOURCE="${CHIP_SGE_MEMORY_RESOURCE:-h_vmem}"
SGE_MEM_RESOURCE_EXPLICIT=false
SCHEDULER_EXTRA="${CHIP_SCHEDULER_EXTRA:-}"
RETRIES="${CHIP_RETRIES:-0}"
LATENCY_WAIT="${CHIP_LATENCY_WAIT:-}"
MAX_JOBS_PER_SECOND="${CHIP_MAX_JOBS_PER_SECOND:-}"
MAX_STATUS_CHECKS_PER_SECOND="${CHIP_MAX_STATUS_CHECKS_PER_SECOND:-}"

POSITIONAL=()
SNAKEMAKE_EXTRA_ARGS=()

timestamp() {
    date '+%Y-%m-%d %H:%M:%S'
}

info() {
    printf '[%s] [INFO] %s\n' "$(timestamp)" "$*"
}

warn() {
    printf '[%s] [WARN] %s\n' "$(timestamp)" "$*" >&2
}

die() {
    printf '[%s] [ERROR] %s\n' "$(timestamp)" "$*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
chip_cuttag_atac_faire workflow launcher (ChIP-seq / CUT&Tag / ATAC-seq / FAIRE-seq)

Usage:
  bash run.sh <project_dir>
  bash run.sh [options]

Assay routing:
  chip / cuttag / atac / faire are routed per sample by the seqtype column of
  samples.csv; mixed-assay projects need no launcher option.

Core options:
  -P, --project DIR         Project working directory (required; created if missing;
                            holds 1.rawdata/ and samples.csv)
  -w, --workdir DIR         Alias of --project (kept for compatibility)
  -c, --config FILE         Project analysis configuration file
      --software FILE       Software/runtime configuration (software.yaml)
  -j, --jobs N              Maximum parallel jobs (default: 3)
      --profile NAME        auto, default, pbs, sge, or slurm
      --queue NAME          Queue name (pbs or sge profile)
      --partition NAME      SLURM partition name (slurm only)
      --memory VALUE        Override per-rule memory for all cluster jobs, e.g. 16G
      --runtime MIN         Override per-rule walltime for all cluster jobs (minutes)
      --sge-mem-resource N  SGE memory resource name (default: h_vmem; sge only)
      --scheduler-extra S   Extra text appended to qsub/sbatch submit command
      --retries N           Retry failed jobs up to N times (default: 0)
      --latency-wait SEC    Override profile filesystem latency wait
      --max-jobs-per-sec N  Limit job submission rate
      --max-status-per-sec N Limit scheduler status checks
  -r, --rename              Normalize common R1/R2 fastq suffixes in 1.rawdata/ to
                            {sample}_1.fq.gz / {sample}_2.fq.gz before running
  -l, --extra-config FILE   Extra config layered last; config.local.yaml in the
                            project directory is picked up automatically
  -n, --dry-run             Build the DAG and show planned jobs without executing
      --validate-only       Parse workflow and configs (snakemake --list-rules) only
      --check-software      Check required executables, R, and R packages; then exit
      --check-r             Check configured Rscript, R version, libraries, and packages; then exit
      --skip-software-check Skip runtime preflight before a real workflow run
      --unlock              Remove a stale Snakemake working-directory lock
      --log FILE            Launcher/Snakemake log path (default: snakemake.logs.txt)
  -q, --quiet               Reduce launcher output; Snakemake output is unchanged
  -h, --help                Show this help message
      --version             Show launcher version

Config layering:
  Repo config/config.yaml, then the project config (-c or <project>/config.yaml),
  then config.local.yaml; later files override earlier keys.

Scheduler selection:
  auto selects SLURM when sbatch is available; a found qsub resolves to SGE when
  SGE_ROOT is set and to PBS otherwise; none means the local default profile.
  RUN_PROFILE can set the same preference.

Any arguments after "--" are passed directly to Snakemake.
EOF
    cat <<'EOF'

Examples:
  bash run.sh /data/project -j 3 --profile pbs
  bash run.sh -P /data/project -c config.yaml --check-software
  bash run.sh -P . -n
  bash run.sh -P . --profile pbs --queue workq --memory 16G --runtime 600
  bash run.sh -P . --unlock
  bash run.sh -P . -- --rerun-triggers mtime
EOF
}

print_version() {
    printf 'run.sh %s\n' "$SCRIPT_VERSION"
}

resolve_path() {
    local value="$1"
    if [[ "$value" = /* ]]; then
        printf '%s\n' "$value"
    else
        printf '%s\n' "$CALL_DIR/$value"
    fi
}

while (($#)); do
    case "$1" in
        --project=*|--project-dir=*|--workdir=*) PROJECT_DIR="${1#*=}"; shift ;;
        --config=*) CONFIG_PATH="${1#*=}"; shift ;;
        --extra-config=*) EXTRA_CONFIG="${1#*=}"; shift ;;
        --software=*) SOFTWARE_PATH="${1#*=}"; shift ;;
        --jobs=*) JOBS="${1#*=}"; shift ;;
        --profile=*) PROFILE_REQUEST="${1#*=}"; shift ;;
        --queue=*) QUEUE="${1#*=}"; shift ;;
        --partition=*) PARTITION="${1#*=}"; shift ;;
        --memory=*|--mem=*) MEMORY="${1#*=}"; shift ;;
        --runtime=*) RUNTIME_MIN="${1#*=}"; shift ;;
        --sge-mem-resource=*) SGE_MEMORY_RESOURCE="${1#*=}"; SGE_MEM_RESOURCE_EXPLICIT=true; shift ;;
        --scheduler-extra=*) SCHEDULER_EXTRA="${1#*=}"; shift ;;
        --retries=*) RETRIES="${1#*=}"; shift ;;
        --latency-wait=*) LATENCY_WAIT="${1#*=}"; shift ;;
        --max-jobs-per-sec=*) MAX_JOBS_PER_SECOND="${1#*=}"; shift ;;
        --max-status-per-sec=*) MAX_STATUS_CHECKS_PER_SECOND="${1#*=}"; shift ;;
        --log=*) LOG_PATH="${1#*=}"; shift ;;
        -P|--project|--project-dir|-w|--workdir) [[ $# -ge 2 ]] || die "$1 requires a value"; PROJECT_DIR="$2"; shift 2 ;;
        -c|--config) [[ $# -ge 2 ]] || die "$1 requires a value"; CONFIG_PATH="$2"; shift 2 ;;
        -l|--extra-config) [[ $# -ge 2 ]] || die "$1 requires a value"; EXTRA_CONFIG="$2"; shift 2 ;;
        --software) [[ $# -ge 2 ]] || die "$1 requires a value"; SOFTWARE_PATH="$2"; shift 2 ;;
        -j|--jobs) [[ $# -ge 2 ]] || die "$1 requires a value"; JOBS="$2"; shift 2 ;;
        --profile) [[ $# -ge 2 ]] || die "$1 requires a value"; PROFILE_REQUEST="$2"; shift 2 ;;
        --queue) [[ $# -ge 2 ]] || die "$1 requires a value"; QUEUE="$2"; shift 2 ;;
        --partition) [[ $# -ge 2 ]] || die "$1 requires a value"; PARTITION="$2"; shift 2 ;;
        --memory|--mem) [[ $# -ge 2 ]] || die "$1 requires a value"; MEMORY="$2"; shift 2 ;;
        --runtime) [[ $# -ge 2 ]] || die "$1 requires a value"; RUNTIME_MIN="$2"; shift 2 ;;
        --sge-mem-resource) [[ $# -ge 2 ]] || die "$1 requires a value"; SGE_MEMORY_RESOURCE="$2"; SGE_MEM_RESOURCE_EXPLICIT=true; shift 2 ;;
        --scheduler-extra) [[ $# -ge 2 ]] || die "$1 requires a value"; SCHEDULER_EXTRA="$2"; shift 2 ;;
        --retries) [[ $# -ge 2 ]] || die "$1 requires a value"; RETRIES="$2"; shift 2 ;;
        --latency-wait) [[ $# -ge 2 ]] || die "$1 requires a value"; LATENCY_WAIT="$2"; shift 2 ;;
        --max-jobs-per-sec) [[ $# -ge 2 ]] || die "$1 requires a value"; MAX_JOBS_PER_SECOND="$2"; shift 2 ;;
        --max-status-per-sec) [[ $# -ge 2 ]] || die "$1 requires a value"; MAX_STATUS_CHECKS_PER_SECOND="$2"; shift 2 ;;
        --log) [[ $# -ge 2 ]] || die "$1 requires a value"; LOG_PATH="$2"; shift 2 ;;
        -r|--rename) RENAME=true; shift ;;
        -n|--dry-run) DRY_RUN=true; shift ;;
        --validate-only) VALIDATE_ONLY=true; shift ;;
        --check-software) CHECK_SOFTWARE_ONLY=true; shift ;;
        --check-r) CHECK_R_ONLY=true; shift ;;
        --skip-software-check) SKIP_SOFTWARE_CHECK=true; shift ;;
        --unlock) UNLOCK=true; shift ;;
        -q|--quiet) QUIET=true; shift ;;
        -h|--help) usage; exit 0 ;;
        --version) print_version; exit 0 ;;
        --) shift; SNAKEMAKE_EXTRA_ARGS+=("$@"); break ;;
        -*) die "Unknown option: $1. Use --help for usage." ;;
        *) POSITIONAL+=("$1"); shift ;;
    esac
done

# Environment-provided SGE memory resource counts as an explicit request.
if [[ -n "${CHIP_SGE_MEMORY_RESOURCE:-}" ]]; then
    SGE_MEM_RESOURCE_EXPLICIT=true
fi

# Backward-compatible positional form (old main_run.sh habit):
#   run.sh <project_dir>
# All other settings must be passed as options.
[[ -n "$PROJECT_DIR" ]] || PROJECT_DIR="${POSITIONAL[0]:-}"
[[ ${#POSITIONAL[@]} -le 1 ]] || die "Too many positional arguments. Only the project directory is accepted; pass the rest as options (see --help)."
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

if [[ -n "$EXTRA_CONFIG" ]]; then
    EXTRA_CONFIG="$(resolve_path "$EXTRA_CONFIG")"
    [[ -f "$EXTRA_CONFIG" ]] || die "Extra configuration file not found: $EXTRA_CONFIG"
    EXTRA_CONFIG="$(cd "$(dirname "$EXTRA_CONFIG")" && pwd)/$(basename "$EXTRA_CONFIG")"
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

[[ -f "$SNAKEFILE" ]] || die "Snakefile not found: $SNAKEFILE"
[[ -f "$RUNTIME_HELPER" ]] || die "Runtime helper not found: $RUNTIME_HELPER"
command -v python3 >/dev/null 2>&1 || die "python3 is required to resolve software.yaml."
RUNTIME_PYTHON="$(command -v python3)"
"$RUNTIME_PYTHON" -c 'import yaml' >/dev/null 2>&1 || die "PyYAML is required by the launcher runtime resolver."
RUNTIME_EXPORTS="$("$RUNTIME_PYTHON" "$RUNTIME_HELPER" export --config "$SOFTWARE_PATH")" || die "Unable to resolve software runtime: $SOFTWARE_PATH"
eval "$RUNTIME_EXPORTS"
export CHIP_SOFTWARE_CONFIG="$SOFTWARE_PATH"

command -v snakemake >/dev/null 2>&1 || die "snakemake was not found after applying software.yaml. Add it to the main environment or PATH."

if command -v timeout >/dev/null 2>&1; then
    SNAKEMAKE_VERSION="$(timeout 8s snakemake --version 2>/dev/null | head -1 || true)"
else
    SNAKEMAKE_VERSION="$(snakemake --version 2>/dev/null | head -1 || true)"
fi
if [[ -z "$SNAKEMAKE_VERSION" ]]; then
    warn "Unable to determine the Snakemake version; continuing anyway."
elif [[ "$SNAKEMAKE_VERSION" =~ ^[0-9]+ ]]; then
    SNAKEMAKE_MAJOR="${SNAKEMAKE_VERSION%%.*}"
    if (( SNAKEMAKE_MAJOR >= 8 )); then
        warn "snakemake 8+ 改用 executor 插件体系（--cluster 语义有变化），集群 profile 实跑前请先实测"
    fi
fi

# CHIP_CONFIG 供 Snakefile 的 configfile: 指令在解析期读取（缺省回落仓库默认配置）。
# 最终生效配置以 --configfile 链为准：snakemake 会把命令行 --configfile 的合并结果
# 再次覆盖到 configfile: 指令加载的值之上，因此链上靠后的文件（config.local.yaml）
# 始终优先；CHIP_CONFIG 只承担解析期种子与缺失回落两个角色，故只指向链上第一个
# 项目级 config。
export CHIP_CONFIG="$CONFIG_FOR_PARSE"

select_profile() {
    local requested
    requested="$(printf '%s' "$PROFILE_REQUEST" | tr '[:upper:]' '[:lower:]')"
    [[ "$requested" == "local" ]] && requested="default"
    if [[ "$requested" == "auto" ]]; then
        if command -v sbatch >/dev/null 2>&1; then
            requested="slurm"
        elif command -v qsub >/dev/null 2>&1; then
            # qsub 二义性：SGE 环境必有 SGE_ROOT；否则按 PBS 处理
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
        pbs) command -v qsub >/dev/null 2>&1 || die "PBS profile selected but qsub was not found in PATH." ;;
        sge) command -v qsub >/dev/null 2>&1 || die "SGE profile selected but qsub was not found in PATH." ;;
        slurm) command -v sbatch >/dev/null 2>&1 || die "SLURM profile selected but sbatch was not found in PATH." ;;
    esac
}

runtime_check() {
    local scope="$1"
    "$RUNTIME_PYTHON" "$RUNTIME_HELPER" check \
        --config "$SOFTWARE_PATH" --scope "$scope"
}

if [[ "$CHECK_SOFTWARE_ONLY" == true ]]; then
    runtime_check all
    info "Software/runtime check completed successfully."
    exit 0
fi
if [[ "$CHECK_R_ONLY" == true ]]; then
    runtime_check r
    info "R runtime check completed successfully."
    exit 0
fi

select_profile

if [[ -n "$QUEUE" && "$PROFILE" != "pbs" && "$PROFILE" != "sge" ]]; then
    die "--queue is only valid with the PBS or SGE profile."
fi
if [[ -n "$PARTITION" && "$PROFILE" != "slurm" ]]; then
    die "--partition is only valid with the SLURM profile."
fi
if [[ "$SGE_MEM_RESOURCE_EXPLICIT" == true && "$PROFILE" != "sge" ]]; then
    die "--sge-mem-resource is only valid with the SGE profile."
fi
if [[ -n "$MEMORY" && "$PROFILE" == "default" ]]; then
    die "--memory is a cluster scheduler request and is not valid with the local profile."
fi
if [[ -n "$RUNTIME_MIN" && "$PROFILE" == "default" ]]; then
    die "--runtime is a cluster scheduler request and is not valid with the local profile."
fi
if [[ -n "$SCHEDULER_EXTRA" && "$PROFILE" == "default" ]]; then
    die "--scheduler-extra is only valid with cluster profiles (pbs, sge, or slurm)."
fi

CLUSTER_CMD=""
case "$PROFILE" in
    pbs)
        if [[ -n "$MEMORY" ]]; then
            PBS_MEM_REQUEST="$MEMORY"
        else
            PBS_MEM_REQUEST='{resources.mem_mb}mb'
        fi
        if [[ -n "$RUNTIME_MIN" ]]; then
            PBS_RUNTIME_REQUEST="$((RUNTIME_MIN * 60))"
        else
            PBS_RUNTIME_REQUEST='{resources.runtime_sec}'
        fi
        CLUSTER_CMD="qsub -V -N {rule} -l select=1:ncpus={threads}:mem=${PBS_MEM_REQUEST} -l walltime=${PBS_RUNTIME_REQUEST} -j oe"
        [[ -n "$QUEUE" ]] && CLUSTER_CMD+=" -q $QUEUE"
        ;;
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

cd "$PROJECT_DIR"
mkdir -p logs

# 可选：统一原始数据命名（perl rename 语法；系统无 rename 则跳过并提示）
if [[ "$RENAME" == true && -d "1.rawdata" ]]; then
    if command -v rename >/dev/null 2>&1; then
        (
            cd 1.rawdata
            rename _R1.fastq.gz _1.fq.gz ./*gz 2>/dev/null || true
            rename _R2.fastq.gz _2.fq.gz ./*gz 2>/dev/null || true
            rename _1.fastq.gz _1.fq.gz ./*gz 2>/dev/null || true
            rename _2.fastq.gz _2.fq.gz ./*gz 2>/dev/null || true
        )
    else
        warn "未找到 rename 命令，跳过重命名；请确保 fastq 命名为 {sample}_1.fq.gz / {sample}_2.fq.gz"
    fi
fi

# 额外配置：显式 -l/--extra-config 优先；否则自动检测项目目录下的 config.local.yaml
if [[ -z "$EXTRA_CONFIG" && -f "config.local.yaml" ]]; then
    EXTRA_CONFIG="$PROJECT_DIR/config.local.yaml"
    info "检测到 config.local.yaml，将叠加覆盖默认配置"
fi

# 配置链：仓库默认 + 项目 config + config.local.yaml，依序传给 snakemake（后者覆盖前者）。
CONFIGFILE_ARGS=(--configfile "$DEFAULT_CONFIG")
[[ -n "$CONFIG_PATH" ]] && CONFIGFILE_ARGS+=(--configfile "$CONFIG_PATH")
[[ -n "$EXTRA_CONFIG" ]] && CONFIGFILE_ARGS+=(--configfile "$EXTRA_CONFIG")

if [[ ! -d "$PROJECT_DIR/1.rawdata" ]]; then
    warn "Raw-data directory is missing: $PROJECT_DIR/1.rawdata. This is acceptable only when all required outputs already exist."
fi

if [[ "$SKIP_SOFTWARE_CHECK" != true && "$DRY_RUN" != true ]]; then
    runtime_check all
elif [[ "$DRY_RUN" == true && "$SKIP_SOFTWARE_CHECK" != true ]]; then
    info "Software/runtime preflight is skipped for dry-run; use --check-software for an explicit check."
fi

if [[ "$QUIET" != true ]]; then
    info "chip_cuttag_atac_faire launcher version: $SCRIPT_VERSION"
    info "Project: $PROJECT_DIR"
    info "Config: ${CONFIG_PATH:-$DEFAULT_CONFIG}"
    [[ -n "$EXTRA_CONFIG" ]] && info "Extra config: $EXTRA_CONFIG (layered last)"
    info "Software config: $SOFTWARE_PATH"
    info "Software environment: ${CHIP_SOFTWARE_TYPE:-system}${CHIP_ENV_PREFIX:+ ($CHIP_ENV_PREFIX)}"
    info "Rscript: ${CHIP_RSCRIPT:-Rscript}"
    [[ -n "${R_LIBS_USER:-}" ]] && info "R libraries: $R_LIBS_USER"
    info "Profile: $PROFILE ($PROFILE_DIR)"
    info "Jobs: $JOBS"
    [[ "$RETRIES" != "0" ]] && info "Retries: $RETRIES"
    [[ -n "$LATENCY_WAIT" ]] && info "Latency wait: ${LATENCY_WAIT}s"
    [[ -n "$QUEUE" ]] && info "Queue: $QUEUE"
    [[ -n "$PARTITION" ]] && info "SLURM partition: $PARTITION"
    [[ "$RENAME" == true ]] && info "Raw-data rename: enabled"
    if [[ -n "$MEMORY" ]]; then
        info "Scheduler memory override: $MEMORY"
    elif [[ "$PROFILE" != "default" ]]; then
        info "Scheduler memory: per-rule resources.mem_mb"
    fi
    if [[ -n "$RUNTIME_MIN" ]]; then
        info "Scheduler runtime override: ${RUNTIME_MIN} min"
    elif [[ "$PROFILE" != "default" ]]; then
        info "Scheduler runtime: per-rule resources.runtime_*"
    fi
    [[ -n "$CLUSTER_CMD" ]] && info "Submit command: $CLUSTER_CMD"
    [[ -n "$SNAKEMAKE_VERSION" ]] && info "Snakemake: $SNAKEMAKE_VERSION"
fi
if [[ "$VALIDATE_ONLY" == true ]]; then
    # 样本表与 config 校验集中在 Snakefile 解析期执行：--list-rules 触发完整解析后即退出。
    info "Validate-only: parsing workflow and configs (snakemake --list-rules)."
    snakemake -s "$SNAKEFILE" "${CONFIGFILE_ARGS[@]}" --profile "$PROFILE_DIR" --list-rules
    info "Validation completed successfully."
    exit 0
fi

SNAKEMAKE_ARGS=(
    -s "$SNAKEFILE"
    --profile "$PROFILE_DIR"
    -j "$JOBS"
    "${CONFIGFILE_ARGS[@]}"
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
        if [[ "$PROFILE" == "pbs" ]]; then
            # PBS 集群模式：回收落在项目目录根的 qsub 输出日志（存在才移动）
            mv ./[a-zA-Z]*.o* ./logs/ 2>/dev/null || true
        fi
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

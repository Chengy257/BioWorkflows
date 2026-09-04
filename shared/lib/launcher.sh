# Shared launcher helpers for the BioWorkflows run.sh scripts.
# Source from a project's run.sh after REPO_DIR and CALL_DIR are defined:
#   source "$BIO_WORKFLOWS_SHARED/lib/launcher.sh"
# Provides: timestamp, info, warn, die, resolve_path.
# Scheduler auto-detection and cluster submit strings stay per project
# (they encode project-specific profile sets and detection policies).

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

resolve_path() {
    local value="$1"
    if [[ "$value" = /* ]]; then
        printf '%s\n' "$value"
    else
        printf '%s\n' "$CALL_DIR/$value"
    fi
}

#!/bin/bash
# ============================================================
# Shared logging and utility functions for CLI pipeline
#
# Source in every run script:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#   source "$SCRIPT_DIR/../lib.sh"
#
# Requires REPO_DIR, CONTAINER, RUNTIME to be set by the caller.
# Recommended: set -euo pipefail in each run script.
# Lint with: shellcheck -x scripts/cli/*/run_*.sh
# ============================================================

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

_PIPELINE_START=$(date +%s)

# ---- Banner and section headers --------------------------------

log_header() {
    local title="$1"
    echo -e "\n${BOLD}══════════════════════════════════════════════════════════${NC}" >&2
    printf "${BOLD}  %-56s${NC}\n" "$title" >&2
    echo -e "${BOLD}══════════════════════════════════════════════════════════${NC}" >&2
}

log_step() {
    local n="$1" msg="$2"
    echo -e "\n${BLUE}${BOLD}▶  Step ${n}: ${msg}${NC}" >&2
    echo -e "${CYAN}   $(date '+%Y-%m-%d %H:%M:%S')${NC}" >&2
    echo -e "${CYAN}   ────────────────────────────────────────${NC}" >&2
}

# ---- Status messages -------------------------------------------

log_ok()   { echo -e "${GREEN}   ✓  $1${NC}" >&2; }
log_warn() { echo -e "${YELLOW}   ⚠  $1${NC}" >&2; }
log_err()  { echo -e "${RED}   ✗  $1${NC}" >&2; }
log_info() { echo -e "   $1" >&2; }

die() {
    log_err "$*"
    log_fail "Pipeline aborted"
    exit 1
}

# ---- Elapsed time ----------------------------------------------

_elapsed() {
    local end
    end=$(date +%s)
    local s=$(( end - _PIPELINE_START ))
    printf "%dm %02ds" $(( s / 60 )) $(( s % 60 ))
}

log_done() {
    local title="$1"
    echo -e "\n${CYAN}   Elapsed: $(_elapsed)${NC}" >&2
    echo -e "\n${GREEN}${BOLD}══════════════════════════════════════════════════════════${NC}" >&2
    printf "${GREEN}${BOLD}  ✓  %-54s${NC}\n" "$title" >&2
    echo -e "${GREEN}${BOLD}══════════════════════════════════════════════════════════${NC}\n" >&2
}

log_fail() {
    local title="$1"
    echo -e "\n${CYAN}   Elapsed: $(_elapsed)${NC}" >&2
    echo -e "\n${RED}${BOLD}══════════════════════════════════════════════════════════${NC}" >&2
    printf "${RED}${BOLD}  ✗  %-54s${NC}\n" "$title" >&2
    echo -e "${RED}${BOLD}══════════════════════════════════════════════════════════${NC}\n" >&2
}

# ---- Prerequisite checks ---------------------------------------

require_file() {
    local f="$1" hint="${2:-}"
    if [ ! -f "$f" ]; then
        log_err "Required file not found: $f"
        [ -n "$hint" ] && log_info "  → $hint"
        exit 1
    fi
}

require_dir() {
    local d="$1" hint="${2:-}"
    if [ ! -d "$d" ]; then
        log_err "Required directory not found: $d"
        [ -n "$hint" ] && log_info "  → $hint"
        exit 1
    fi
}

# ---- Docker/Podman runners -------------------------------------
# container_run "command args..."         — workdir /project
# container_run_in "/project/subdir" "command args..."

container_run() {
    "$RUNTIME" run --rm \
        -v "$REPO_DIR":/project \
        --workdir /project \
        --entrypoint bash \
        "$CONTAINER" \
        --norc -c "eval \"\$(pixi shell-hook)\" && $*"
}

container_run_in() {
    local workdir="$1"; shift
    "$RUNTIME" run --rm \
        -v "$REPO_DIR":/project \
        --workdir "$workdir" \
        --entrypoint bash \
        "$CONTAINER" \
        --norc -c "eval \"\$(pixi shell-hook)\" && $*"
}

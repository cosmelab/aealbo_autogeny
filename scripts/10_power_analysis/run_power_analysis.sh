#!/bin/bash
# ============================================================
# Power Analysis for Selection Scan Methods
# ============================================================
# Usage: bash scripts/cli/10_power_analysis/run_power_analysis.sh [repo_dir]
# CONTAINER_RUNTIME=podman bash scripts/cli/10_power_analysis/run_power_analysis.sh
#
# Note: Power analysis is provided as supplementary context and
#       was not included in the final manuscript.
# Standalone — no data dependencies.
# ============================================================
set -euo pipefail

REPO_DIR="${1:-$HOME/projects/albopictus-autogeny}"
CONTAINER="ghcr.io/cosmelab/albopictus-autogeny:latest"
RUNTIME="${CONTAINER_RUNTIME:-docker}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/cli/lib.sh
source "$SCRIPT_DIR/../lib.sh"

cd "$REPO_DIR"
mkdir -p logs output/power_analysis

log_header "Power Analysis  |  Runtime: $RUNTIME"
log_warn "Note: This analysis was not included in the final manuscript."

log_step 1 "Simulate and compute detection power"
container_run "Rscript scripts/cli/10_power_analysis/power_analysis.R /project"
log_ok "Power analysis → output/power_analysis/"

log_done "Power Analysis Complete"

#!/bin/bash
# ============================================================
# Pairwise pcadapt Analysis
# ============================================================
# Usage: bash scripts/cli/02_selection_scans/run_pairwise_pcadapt.sh [repo_dir]
# CONTAINER_RUNTIME=podman bash scripts/cli/02_selection_scans/run_pairwise_pcadapt.sh
#
# Note: Pairwise pcadapt between NON-AUTO and NON-AUTO-FIELD
#       was not used in the final manuscript.
# Prerequisites: output/quality_control/file7  (from 01_qc)
# ============================================================
set -euo pipefail

REPO_DIR="${1:-$HOME/projects/albopictus-autogeny}"
CONTAINER="ghcr.io/cosmelab/albopictus-autogeny:latest"
RUNTIME="${CONTAINER_RUNTIME:-docker}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/cli/lib.sh
source "$SCRIPT_DIR/../lib.sh"

cd "$REPO_DIR"
mkdir -p logs output/pcadapt

log_header "Pairwise pcadapt  |  Runtime: $RUNTIME"

require_file "output/quality_control/file7.bed" \
    "Run 01_qc first: bash scripts/cli/01_qc/run_qc.sh"

log_step 1 "Pairwise pcadapt (all population pairs)"
container_run "Rscript scripts/cli/02_selection_scans/pairwise_pcadapt.R \
    output/quality_control \
    output/pcadapt"
log_ok "Pairwise pcadapt complete → output/pcadapt/"

log_done "Pairwise pcadapt Complete"

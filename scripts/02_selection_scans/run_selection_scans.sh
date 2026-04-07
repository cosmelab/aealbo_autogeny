#!/bin/bash
# ============================================================
# Selection Scans (OutFLANK + pcadapt) → Files S2, S3
# ============================================================
# Usage: bash scripts/cli/02_selection_scans/run_selection_scans.sh [repo_dir]
# CONTAINER_RUNTIME=podman bash scripts/cli/02_selection_scans/run_selection_scans.sh
#
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
mkdir -p logs output/selection_scans output/outflank output/pcadapt

log_header "Selection Scans  |  Runtime: $RUNTIME"

require_file "output/quality_control/file7.bed" \
    "Run 01_qc first: bash scripts/cli/01_qc/run_qc.sh"
require_file "data/files/intergenic_SNPs.txt"

log_step 1 "OutFLANK + 3-way pcadapt"
container_run "Rscript scripts/cli/02_selection_scans/selection_scans.R \
    output/quality_control/file7 \
    output/selection_scans \
    data/files/intergenic_SNPs.txt"
log_ok "OutFLANK complete → output/selection_scans/"

log_step 2 "Pairwise pcadapt (AUTO vs NON-AUTO vs NON-AUTO-FIELD)"
container_run "Rscript scripts/cli/02_selection_scans/pairwise_pcadapt.R \
    output/quality_control \
    output/pcadapt"
log_ok "Pairwise pcadapt complete → output/pcadapt/"

log_done "Selection Scans Complete"

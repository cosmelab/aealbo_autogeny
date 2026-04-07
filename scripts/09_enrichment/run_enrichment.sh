#!/bin/bash
# ============================================================
# Functional Enrichment Analysis
# ============================================================
# Usage: bash scripts/cli/09_enrichment/run_enrichment.sh [repo_dir]
# CONTAINER_RUNTIME=podman bash scripts/cli/09_enrichment/run_enrichment.sh
#
# Note: Functional enrichment analysis was not included in the
#       final manuscript.
# Prerequisites: output/ldna/  (from 05_ldna)
# ============================================================
set -euo pipefail

REPO_DIR="${1:-$HOME/projects/albopictus-autogeny}"
CONTAINER="ghcr.io/cosmelab/albopictus-autogeny:latest"
RUNTIME="${CONTAINER_RUNTIME:-docker}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/cli/lib.sh
source "$SCRIPT_DIR/../lib.sh"

cd "$REPO_DIR"
mkdir -p logs output/enrichment

log_header "Functional Enrichment Analysis  |  Runtime: $RUNTIME"
log_warn "Note: This analysis was not included in the final manuscript."

require_dir "output/ldna" \
    "Run 05_ldna first: bash scripts/cli/05_ldna/run_ldna.sh"

log_step 1 "Enrichment analysis"
container_run "Rscript scripts/cli/09_enrichment/enrichment_analysis.R /project"
log_ok "Enrichment results → output/enrichment/"

log_done "Enrichment Analysis Complete"

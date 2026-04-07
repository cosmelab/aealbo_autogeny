#!/bin/bash
# ============================================================
# LDna Network Analysis → Files S4a–S4e
# ============================================================
# Usage: bash scripts/cli/05_ldna/run_ldna.sh [repo_dir] [cores]
# CONTAINER_RUNTIME=podman bash scripts/cli/05_ldna/run_ldna.sh
#
# Note: Requires ~32GB RAM. Use cores=2 on 32GB machines.
# Prerequisites: output/quality_control/file7  (from 01_qc)
# ============================================================
set -euo pipefail

REPO_DIR="${1:-$HOME/projects/albopictus-autogeny}"
CORES="${2:-2}"
CONTAINER="ghcr.io/cosmelab/albopictus-autogeny:latest"
RUNTIME="${CONTAINER_RUNTIME:-docker}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/cli/lib.sh
source "$SCRIPT_DIR/../lib.sh"

cd "$REPO_DIR"
mkdir -p logs output/ldna

log_header "LDna Network Analysis  |  Runtime: $RUNTIME  |  Cores: $CORES"
log_warn "Requires ~32GB RAM. This step takes ~45 minutes."

require_file "output/quality_control/file7.bed" \
    "Run 01_qc first: bash scripts/cli/01_qc/run_qc.sh"

log_step 1 "LD matrix computation + network analysis (all chromosomes)"
container_run "Rscript scripts/cli/05_ldna/ldna_analysis.R /project $CORES"
log_ok "LDna analysis complete → output/ldna/"

log_done "LDna Analysis Complete"

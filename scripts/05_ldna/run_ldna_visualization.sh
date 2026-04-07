#!/bin/bash
# ============================================================
# LDna Visualization → Files S4a–S4e
# ============================================================
# Usage: bash scripts/cli/05_ldna/run_ldna_visualization.sh [repo_dir]
# CONTAINER_RUNTIME=podman bash scripts/cli/05_ldna/run_ldna_visualization.sh
#
# Produces:
#   output/ldna/soc_networks/   — native plotLDnetwork SOC plots (per pop/chr)
#   output/ldna/                — cluster position PDFs (ggplot2 rectangles)
#
# Prerequisites: output/ldna/  (from run_ldna.sh)
# ============================================================
set -euo pipefail

REPO_DIR="${1:-$HOME/projects/albopictus-autogeny}"
CONTAINER="ghcr.io/cosmelab/albopictus-autogeny:latest"
RUNTIME="${CONTAINER_RUNTIME:-docker}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/cli/lib.sh
source "$SCRIPT_DIR/../lib.sh"

cd "$REPO_DIR"
mkdir -p logs output/ldna/soc_networks

log_header "LDna Visualization  |  Runtime: $RUNTIME"

require_dir "output/ldna" \
    "Run LDna analysis first: bash scripts/cli/05_ldna/run_ldna.sh"
require_file "output/ldna/pop/chr2/AUT.rds" \
    "Run LDna analysis first: bash scripts/cli/05_ldna/run_ldna.sh"

log_step 1 "SOC network plots (native plotLDnetwork — per population × chromosome)"
log_info "   Plots each SOC as an LD network: SNPs=nodes, LD=edges, position=color"
container_run "Rscript scripts/cli/05_ldna/ldna_soc_networks.R /project 50 100"
log_ok "SOC networks → output/ldna/soc_networks/"

log_step 2 "Cluster position plots (ggplot2 rectangles per chromosome)"
container_run "Rscript scripts/cli/05_ldna/ldna_visualization.R /project 100 0.1"
log_ok "Cluster plots → output/ldna/"

log_done "LDna Visualization Complete"

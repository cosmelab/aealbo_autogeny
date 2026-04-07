#!/bin/bash
# ============================================================
# LDna Network Plots → Files S4a–S4e
# ============================================================
# Usage: bash scripts/cli/05_ldna/run_network_plots.sh [repo_dir]
# CONTAINER_RUNTIME=podman bash scripts/cli/05_ldna/run_network_plots.sh
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
mkdir -p logs output/ldna/network_plots

log_header "LDna Network Plots  |  Runtime: $RUNTIME"

require_dir "output/ldna" \
    "Run LDna analysis first: bash scripts/cli/05_ldna/run_ldna.sh"

log_step 1 "Cluster comparison and network plots"
container_run "Rscript scripts/cli/05_ldna/ldna_cluster_comparison.R \
    --ldna-dir /project/output/ldna --outdir /project/output/ldna/network_plots"
log_ok "Network plots → output/ldna/network_plots/"

log_done "Network Plots Complete"

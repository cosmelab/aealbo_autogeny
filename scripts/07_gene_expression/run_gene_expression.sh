#!/bin/bash
# ============================================================
# Gene Expression Intersection → File S5
# ============================================================
# Usage: bash scripts/cli/07_gene_expression/run_gene_expression.sh [repo_dir]
# CONTAINER_RUNTIME=podman bash scripts/cli/07_gene_expression/run_gene_expression.sh
#
# Note: Full SNP-to-gene coordinate mapping (scaffold → chromosome)
#       is implemented in notebooks/File_S5.Intersection_gene_expression_SNPs.Rmd.
#       This script provides a summary only.
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
mkdir -p logs output/gene_expression

log_header "Gene Expression Intersection  |  Runtime: $RUNTIME"

require_dir "output/ldna" \
    "Run 05_ldna first: bash scripts/cli/05_ldna/run_ldna.sh"

log_step 1 "Intersect candidate SNPs with gene expression data"
container_run "python3 scripts/cli/07_gene_expression/gene_expression_intersection.py \
    --project-dir /project"
log_ok "Intersection summary → output/gene_expression/"

log_done "Gene Expression Intersection Complete"

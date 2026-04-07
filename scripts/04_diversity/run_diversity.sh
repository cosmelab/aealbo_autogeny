#!/bin/bash
# ============================================================
# Diversity, Fst, and Frequency Analysis → Files S7, S8
# ============================================================
# Usage: bash scripts/cli/04_diversity/run_diversity.sh [repo_dir]
# CONTAINER_RUNTIME=podman bash scripts/cli/04_diversity/run_diversity.sh
#
# Prerequisites:
#   output/quality_control/file7  (from 01_qc)
#   output/snpeff/SNPs_158.txt    (from 03_annotation)
#   output/ldna/files/file1       (from 05_ldna)
#
# Note: Sex-stratified Fst (S8 Sections 2 and 4) was not used
#       in the final manuscript and is not implemented here.
# ============================================================
set -euo pipefail

REPO_DIR="${1:-$HOME/projects/albopictus-autogeny}"
CONTAINER="ghcr.io/cosmelab/albopictus-autogeny:latest"
RUNTIME="${CONTAINER_RUNTIME:-docker}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/cli/lib.sh
source "$SCRIPT_DIR/../lib.sh"

cd "$REPO_DIR"
mkdir -p logs output/diversity output/fst/figures output/frequencies/figures

log_header "Diversity, Fst, and Frequency Analysis  |  Runtime: $RUNTIME"

require_file "output/quality_control/file7.bed" \
    "Run 01_qc first: bash scripts/cli/01_qc/run_qc.sh"
require_file "output/snpeff/SNPs_158.txt" \
    "Run 03_annotation first: bash scripts/cli/03_annotation/run_snpeff.sh"
require_file "output/ldna/files/file1.bed" \
    "Run 05_ldna first: bash scripts/cli/05_ldna/run_ldna.sh"

log_step 1 "He/Ho genetic diversity"
container_run "Rscript scripts/cli/04_diversity/genetic_diversity.R \
    output/quality_control/file7 output/diversity"
log_ok "Diversity → output/diversity/"

log_step 2 "Sliding window Fst"
container_run "bash scripts/cli/04_diversity/sliding_window_fst.sh"
log_ok "Sliding window Fst → output/fst/"

log_step 3 "Pairwise Fst (StAMPP)"
container_run "Rscript scripts/cli/04_diversity/pairwise_fst.R \
    output/ldna/files/file1 output/fst"
log_ok "Pairwise Fst → output/fst/figures/pairwise_estimates.pdf"

log_step 4 "Allele and genotype frequency plots (17 candidate SNPs)"
container_run "Rscript scripts/cli/04_diversity/allele_freq_plots.R \
    output/frequencies"
log_ok "Frequency plots → output/frequencies/figures/"

log_done "Diversity Analysis Complete"

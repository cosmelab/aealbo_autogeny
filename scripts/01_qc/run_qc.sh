#!/bin/bash
# ============================================================
# Quality Control Pipeline → File S1
# ============================================================
# Usage: bash scripts/cli/01_qc/run_qc.sh [repo_dir]
# CONTAINER_RUNTIME=podman bash scripts/cli/01_qc/run_qc.sh
#
# Output: output/quality_control/file7  (60 samples, ~33,836 SNPs)
# ============================================================
set -euo pipefail

REPO_DIR="${1:-$HOME/projects/albopictus-autogeny}"
CONTAINER="ghcr.io/cosmelab/albopictus-autogeny:latest"
RUNTIME="${CONTAINER_RUNTIME:-docker}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/cli/lib.sh
source "$SCRIPT_DIR/../lib.sh"

cd "$REPO_DIR"
mkdir -p logs output/quality_control output/qc data/qc

log_header "Quality Control Pipeline  |  Runtime: $RUNTIME"

require_file "data/genotype_calls/autogenous.vcf" \
    "Place raw VCF at data/genotype_calls/autogenous.vcf"

log_step 1 "Convert VCF to PLINK"
container_run "python3 scripts/cli/01_qc/convert_vcf_to_plink.py"
log_ok "VCF converted to PLINK format"

log_step 2 "Create chromosomal scale"
container_run "python3 scripts/cli/01_qc/create_chromosomal_scale.py"
log_ok "Chromosomal scale ready"

log_step 3 "Run QC pipeline (7-step filter)"
container_run "python3 scripts/cli/01_qc/quality_control.py"
log_ok "QC filters applied"

log_step 4 "Generate QC report"
container_run "python3 scripts/cli/01_qc/visualize_qc_summary.py"
log_ok "Report → output/qc/qc_report.html"

if [ -f "output/quality_control/file7.fam" ]; then
    log_ok "Samples : $(wc -l < output/quality_control/file7.fam)"
    log_ok "SNPs    : $(wc -l < output/quality_control/file7.bim)"
fi

log_done "Quality Control Complete"

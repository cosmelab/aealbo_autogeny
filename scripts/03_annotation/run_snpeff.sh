#!/bin/bash
# ============================================================
# SnpEff Functional Annotation → File S6
# ============================================================
# Usage: bash scripts/cli/03_annotation/run_snpeff.sh [repo_dir]
# CONTAINER_RUNTIME=podman bash scripts/cli/03_annotation/run_snpeff.sh
#
# Prerequisites:
#   output/quality_control/file7  (from 01_qc)
#   output/selection_scans/       (from 02_selection_scans)
# ============================================================
set -euo pipefail

REPO_DIR="${1:-$HOME/projects/albopictus-autogeny}"
CONTAINER="ghcr.io/cosmelab/albopictus-autogeny:latest"
RUNTIME="${CONTAINER_RUNTIME:-docker}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/cli/lib.sh
source "$SCRIPT_DIR/../lib.sh"

cd "$REPO_DIR"
mkdir -p logs output/snpeff output/snpeff/database/data/albo_v3

log_header "SnpEff Annotation  |  Runtime: $RUNTIME"

require_file "output/quality_control/file7.bed" \
    "Run 01_qc first: bash scripts/cli/01_qc/run_qc.sh"
require_file "data/files/genes.gff"
require_file "data/genome/albo.fasta.gz"

log_step 1 "Prepare SnpEff database files"
cp data/files/genes.gff output/snpeff/database/data/albo_v3/genes.gff
zcat data/genome/albo.fasta.gz > output/snpeff/database/data/albo_v3/sequences.fa
cat > output/snpeff/database/snpEff.config << 'EOF'
data.dir = ./data/
albo_v3.genome : Aedes_albopictus_AalbF3
EOF
log_ok "Database files prepared"

log_step 2 "Build SnpEff database"
container_run_in "/project/output/snpeff/database" \
    "snpEff build -gff3 -v albo_v3 -c snpEff.config"
log_ok "SnpEff database built"

log_step 3 "Export QC data to VCF"
container_run "plink2 \
    --bfile output/quality_control/file7 --export vcf \
    --out output/snpeff/all_snps --silent"
log_ok "VCF exported → output/snpeff/all_snps.vcf"

log_step 4 "Run SnpEff annotation"
container_run_in "/project/output/snpeff/database" \
    "snpEff ann -v albo_v3 -c snpEff.config \
        /project/output/snpeff/all_snps.vcf \
        > /project/output/snpeff/all_snps_annotated.vcf"
mv output/snpeff/database/snpEff_genes.txt output/snpeff/ 2>/dev/null || true
mv output/snpeff/database/snpEff_summary.html output/snpeff/ 2>/dev/null || true
log_ok "Annotation complete → output/snpeff/all_snps_annotated.vcf"

log_step 5 "Extract outlier SNP annotations"
container_run "python3 scripts/cli/03_annotation/extract_annotations.py \
    --project-dir /project"
log_ok "Outlier annotations → output/snpeff/SNPs_158.txt"

log_done "SnpEff Annotation Complete"

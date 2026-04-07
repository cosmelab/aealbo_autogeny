#!/bin/bash
# ============================================================
# Albopictus Autogeny — Full CLI Pipeline
# ============================================================
# Usage:
#   bash scripts/cli/run_pipeline.sh [repo_dir]           # run all
#   bash scripts/cli/run_pipeline.sh [repo_dir] --dry-run # check only
#
# CONTAINER_RUNTIME=podman bash scripts/cli/run_pipeline.sh
#
# Pipeline order:
#   01_qc → 02_selection_scans → 05_ldna → 03_annotation
#        → 07_gene_expression → 04_diversity
#
# > Population label key: AUT = AUTO, MAN = NON-AUTO, NEW = NON-AUTO-FIELD
# ============================================================
set -euo pipefail

REPO_DIR="${1:-$HOME/projects/albopictus-autogeny}"
DRY_RUN=false
for arg in "$@"; do [[ "$arg" == "--dry-run" ]] && DRY_RUN=true; done

CONTAINER="ghcr.io/cosmelab/albopictus-autogeny:latest"
RUNTIME="${CONTAINER_RUNTIME:-docker}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/cli/lib.sh
source "$SCRIPT_DIR/lib.sh"

cd "$REPO_DIR"

# ---- Helpers ---------------------------------------------------

_ok_file()  { [ -f "$1" ] && echo "done" || echo "pending"; }
_ok_dir()   { [ -d "$1" ] && ls -A "$1" &>/dev/null && [ "$(ls -A "$1")" ] && echo "done" || echo "pending"; }

_status_icon() {
    if [ "$1" = "done" ]; then
        echo -e "${GREEN}✓ done   ${NC}"
    else
        echo -e "${YELLOW}○ pending${NC}"
    fi
}

run_step() {
    local label="$1" script="$2"; shift 2
    if $DRY_RUN; then
        return
    fi
    log_step "" "$label"
    bash "$SCRIPT_DIR/$script" "$REPO_DIR" "$@"
}

# ---- Pre-flight checks -----------------------------------------

log_header "Albopictus Autogeny Pipeline  |  Runtime: $RUNTIME"

echo -e "\n${BOLD}  Checking prerequisites...${NC}\n"

# Container
if "$RUNTIME" image inspect "$CONTAINER" &>/dev/null; then
    log_ok "Container image available: $CONTAINER"
else
    log_err "Container image not found: $CONTAINER"
    log_info "  → Pull with: $RUNTIME pull $CONTAINER"
    exit 1
fi

# Raw data
require_file "$REPO_DIR/data/genotype_calls/autogenous.vcf" \
    "Place raw VCF at data/genotype_calls/autogenous.vcf"
log_ok "Raw VCF present"

require_file "$REPO_DIR/data/files/intergenic_SNPs.txt"
require_file "$REPO_DIR/data/files/genes.gff"
require_file "$REPO_DIR/data/genome/albo.fasta.gz"
log_ok "Reference data files present"

# Disk space (warn if <20 GB free)
FREE_GB=$(df -BG "$REPO_DIR" | awk 'NR==2 {gsub("G",""); print $4}')
if [ "$FREE_GB" -lt 20 ]; then
    log_warn "Low disk space: ${FREE_GB}GB free (LDna generates ~5-10 GB of LD files)"
else
    log_ok "Disk space: ${FREE_GB}GB free"
fi

# ---- Pipeline status -------------------------------------------

s1=$(_ok_file "$REPO_DIR/output/quality_control/file7.bed")
s2=$(_ok_dir  "$REPO_DIR/output/selection_scans")
s3=$(_ok_file "$REPO_DIR/output/snpeff/SNPs_158.txt")
s4=$(_ok_file "$REPO_DIR/output/ldna/files/file1.bed")
s5=$(_ok_dir  "$REPO_DIR/output/ldna/clusters")
s6=$(_ok_dir  "$REPO_DIR/output/gene_expression")
s7=$(_ok_dir  "$REPO_DIR/output/diversity")
s8=$(_ok_file "$REPO_DIR/output/fst/pop_pairwise2_df.csv")

echo -e "\n${BOLD}  Pipeline steps:${NC}\n"
echo -e "  $(_status_icon "$s1")  Step 1 · Quality Control          → output/quality_control/file7"
echo -e "  $(_status_icon "$s2")  Step 2 · Selection Scans           → output/selection_scans/"
echo -e "  $(_status_icon "$s3")  Step 3 · Functional Annotation     → output/snpeff/SNPs_158.txt"
echo -e "  $(_status_icon "$s4")  Step 4 · LDna Network Analysis     → output/ldna/files/file1"
echo -e "  $(_status_icon "$s5")  Step 5 · LDna Visualization        → output/ldna/clusters/"
echo -e "  $(_status_icon "$s6")  Step 6 · Gene Expression           → output/gene_expression/"
echo -e "  $(_status_icon "$s7")  Step 7 · He/Ho Genetic Diversity   → output/diversity/"
echo -e "  $(_status_icon "$s8")  Step 8 · Fst + Frequency Plots     → output/fst/, output/frequencies/"

echo ""
if $DRY_RUN; then
    log_warn "Dry run — nothing was executed."
    echo -e "  Run without --dry-run to execute pending steps.\n"
    exit 0
fi

# ---- Execution -------------------------------------------------

log_header "Starting Pipeline Execution  |  Runtime: $RUNTIME"

[ "$s1" = "pending" ] && run_step "Quality Control"        "01_qc/run_qc.sh"
[ "$s2" = "pending" ] && run_step "Selection Scans"        "02_selection_scans/run_selection_scans.sh"
[ "$s3" = "pending" ] && run_step "Functional Annotation"  "03_annotation/run_snpeff.sh"
[ "$s4" = "pending" ] && run_step "LDna Analysis"          "05_ldna/run_ldna.sh" "${LDNA_CORES:-2}"
[ "$s5" = "pending" ] && run_step "LDna Visualization"     "05_ldna/run_ldna_visualization.sh"
[ "$s6" = "pending" ] && run_step "Gene Expression"        "07_gene_expression/run_gene_expression.sh"
[ "$s7" = "pending" ] || [ "$s8" = "pending" ] && run_step "Diversity + Fst" "04_diversity/run_diversity.sh"

log_done "All Pipeline Steps Complete"

#!/usr/bin/env bash
# Download input data for aealbo_autogeny analysis
# Data archived at Zenodo: https://zenodo.org/records/19451616
# DOI: 10.5281/zenodo.19451616
set -euo pipefail

ZENODO_RECORD="19451616"
DATA_DIR="${1:-data}"      # optional: pass target directory as first argument

BASE_URL="https://zenodo.org/record/${ZENODO_RECORD}/files"

# ---- helpers ----------------------------------------------------------------

download_file() {
    local name="$1"
    local dest="$2"
    if [[ -f "$dest" ]]; then
        echo "  [skip] $dest already exists"
        return
    fi
    echo "  [download] $name → $dest"
    mkdir -p "$(dirname "$dest")"
    wget -q --show-progress -O "$dest" "${BASE_URL}/${name}"
}

# ---- VCF files --------------------------------------------------------------

echo "=== Downloading genotype calls ==="
download_file "autogenous.vcf" "${DATA_DIR}/genotype_calls/autogenous.vcf"
download_file "aut.vcf"        "${DATA_DIR}/genotype_calls/aut.vcf"
download_file "man.vcf"        "${DATA_DIR}/genotype_calls/man.vcf"
download_file "new.vcf"        "${DATA_DIR}/genotype_calls/new.vcf"

# ---- Annotation and support files -------------------------------------------

echo "=== Downloading annotation files ==="
download_file "MANvsAUTO_sig_mRNAs.csv" \
    "${DATA_DIR}/files/MANvsAUTO_sig_mRNAs.csv"
download_file "chip_ann.txt" \
    "${DATA_DIR}/files/chip_ann.txt"
download_file "genes.gff" \
    "${DATA_DIR}/files/genes.gff"
download_file "albopictus_SNPs_fail_segregation.txt" \
    "${DATA_DIR}/files/albopictus_SNPs_fail_segregation.txt"
download_file "intergenic_SNPs.txt" \
    "${DATA_DIR}/files/intergenic_SNPs.txt"

echo "=== Downloading genome support files ==="
download_file "scaffold_sizes.txt" \
    "${DATA_DIR}/genome/scaffold_sizes.txt"

# ---- Reference genome (external) --------------------------------------------

echo ""
echo "=== Reference genome (AalbF3) ==="
echo "The reference genome is not archived here — download from NCBI:"
echo ""
echo "  wget -P ${DATA_DIR}/genome/ \\"
echo "    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/006/496/715/GCA_006496715.1_AalbF3/GCA_006496715.1_AalbF3_genomic.fna.gz"
echo "  mv ${DATA_DIR}/genome/GCA_006496715.1_AalbF3_genomic.fna.gz ${DATA_DIR}/genome/albo.fasta.gz"
echo ""
echo "  # OR use the scripts/cli/00_setup pipeline which handles this automatically."

# ---- Done -------------------------------------------------------------------

echo ""
echo "=== Data download complete ==="
echo "  Genotype VCFs : ${DATA_DIR}/genotype_calls/"
echo "  Annotation    : ${DATA_DIR}/files/"
echo "  Genome support: ${DATA_DIR}/genome/"
echo ""
echo "Next step: bash scripts/render_notebooks.sh"

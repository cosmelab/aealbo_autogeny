#!/usr/bin/env bash
# Upload input data to Zenodo via REST API (no tarball — individual files)
#
# Usage:
#   On Spark (token from vault):
#     eval $(python3 ~/projects/tracking-hub/scripts/05_load_vault.py --export)
#     bash scripts/01_upload_zenodo.sh
#
#   Locally (set token manually):
#     export ZENODO_TOKEN="your_token_here"
#     bash scripts/01_upload_zenodo.sh
#
#   Sandbox test (safe, creates real deposit on sandbox.zenodo.org):
#     export ZENODO_SANDBOX=1
#     bash scripts/01_upload_zenodo.sh
#
# Token scopes required: deposit:write
# Get a token at: https://zenodo.org/account/settings/applications/tokens/new/
# On Spark the token is stored in the Supabase vault as ZENODO_UCR_API_KEY

set -euo pipefail

# ---- Config -----------------------------------------------------------------

# Allow sandbox mode for testing
if [[ "${ZENODO_SANDBOX:-0}" == "1" ]]; then
    ZENODO_API="https://sandbox.zenodo.org/api"
    echo "[INFO] SANDBOX MODE — using sandbox.zenodo.org (separate token required)"
else
    ZENODO_API="https://zenodo.org/api"
fi

# Token: prefer ZENODO_TOKEN; fall back to ZENODO_UCR_API_KEY (vault name on Spark)
TOKEN="${ZENODO_TOKEN:-${ZENODO_UCR_API_KEY:-}}"

if [[ -z "$TOKEN" ]]; then
    echo "ERROR: No Zenodo token found."
    echo ""
    echo "On Spark (from vault):"
    echo "  eval \$(python3 ~/projects/tracking-hub/scripts/05_load_vault.py --export)"
    echo "  bash scripts/01_upload_zenodo.sh"
    echo ""
    echo "Locally:"
    echo "  export ZENODO_TOKEN='your_token_here'"
    echo "  bash scripts/01_upload_zenodo.sh"
    exit 1
fi

# Files to upload: "source_path|zenodo_filename|description"
# Edit paths to match where the data lives on your machine
DATA_ROOT="${DATA_ROOT:-data}"
FILES=(
    "${DATA_ROOT}/genotype_calls/autogenous.vcf|autogenous.vcf|All populations combined VCF"
    "${DATA_ROOT}/genotype_calls/aut.vcf|aut.vcf|AUTO population VCF"
    "${DATA_ROOT}/genotype_calls/man.vcf|man.vcf|NON-AUTO population VCF"
    "${DATA_ROOT}/genotype_calls/new.vcf|new.vcf|NON-AUTO-FIELD population VCF"
    "${DATA_ROOT}/files/MANvsAUTO_sig_mRNAs.csv|MANvsAUTO_sig_mRNAs.csv|Differentially expressed genes (File S5 input)"
    "${DATA_ROOT}/files/chip_ann.txt|chip_ann.txt|SNP chip annotations (File S6 input)"
    "${DATA_ROOT}/files/genes.gff|genes.gff|Gene annotations for SnpEff"
    "${DATA_ROOT}/files/albopictus_SNPs_fail_segregation.txt|albopictus_SNPs_fail_segregation.txt|QC exclude list"
    "${DATA_ROOT}/files/intergenic_SNPs.txt|intergenic_SNPs.txt|Intergenic SNP filter list"
    "${DATA_ROOT}/genome/scaffold_sizes.txt|scaffold_sizes.txt|Scaffold lengths"
)

# ---- Metadata ---------------------------------------------------------------

METADATA='{
  "metadata": {
    "title": "Data for: Genomic basis of autogeny in Aedes albopictus (Sturiale et al., BMC Biology)",
    "upload_type": "dataset",
    "description": "Input data for GWAS and linkage disequilibrium network analysis of autogeny in Aedes albopictus. Includes VCF files for three populations (AUTO, NON-AUTO, NON-AUTO-FIELD), SNP chip annotations, gene annotations, and QC filter lists. Analysis code at https://github.com/cosmelab/aealbo_autogeny.",
    "creators": [
      {"name": "Sturiale, Stephanie L.", "affiliation": "Georgetown University"},
      {"name": "Heilig, Martha C.", "affiliation": "Georgetown University"},
      {"name": "Aardema, Matthew L.", "affiliation": "Montclair State University"},
      {"name": "Cosme, Luciano V.", "affiliation": "Georgetown University"},
      {"name": "Armbruster, Peter A.", "affiliation": "Georgetown University"}
    ],
    "keywords": ["Aedes albopictus", "autogeny", "GWAS", "SNP chip", "linkage disequilibrium", "selection scan", "mosquito"],
    "license": "cc-by",
    "access_right": "open",
    "related_identifiers": [
      {
        "relation": "isSupplementTo",
        "identifier": "https://github.com/cosmelab/aealbo_autogeny",
        "scheme": "url",
        "resource_type": "software"
      }
    ]
  }
}'

# ---- Functions --------------------------------------------------------------

check_file() {
    local path="$1"
    if [[ ! -f "$path" ]]; then
        echo "  [MISSING] $path — skipping"
        return 1
    fi
    return 0
}

upload_file() {
    local path="$1"
    local name="$2"
    local bucket_url="$3"
    local size
    size=$(du -h "$path" | cut -f1)
    echo "  Uploading $name ($size)..."
    curl --progress-bar \
        -X PUT "${bucket_url}/${name}" \
        -H "Authorization: Bearer ${TOKEN}" \
        --upload-file "$path" | cat
    echo ""
}

# ---- Main -------------------------------------------------------------------

echo "=== Zenodo Upload — aealbo_autogeny data ==="
echo "API: ${ZENODO_API}"
echo ""

# Step 1: Create deposition
echo "Step 1: Creating deposition..."
RESPONSE=$(curl -s -X POST "${ZENODO_API}/deposit/depositions" \
    -H "Authorization: Bearer ${TOKEN}" \
    -H "Content-Type: application/json" \
    -d '{}')

DEPOSITION_ID=$(echo "$RESPONSE" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('id','ERROR'))")
BUCKET_URL=$(echo "$RESPONSE"   | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('links',{}).get('bucket','ERROR'))")

if [[ "$DEPOSITION_ID" == "ERROR" ]] || [[ "$BUCKET_URL" == "ERROR" ]]; then
    echo "ERROR: Could not create deposition."
    echo "Response: $RESPONSE"
    exit 1
fi

echo "  Deposition ID: ${DEPOSITION_ID}"
echo "  Bucket URL:    ${BUCKET_URL}"
echo ""

# Step 2: Upload each file
echo "Step 2: Uploading files..."
UPLOADED=0
SKIPPED=0
for entry in "${FILES[@]}"; do
    IFS='|' read -r src name desc <<< "$entry"
    if check_file "$src"; then
        upload_file "$src" "$name" "$BUCKET_URL"
        (( UPLOADED++ )) || true
    else
        (( SKIPPED++ )) || true
    fi
done
echo "  Uploaded: ${UPLOADED}  Skipped (missing): ${SKIPPED}"
echo ""

# Step 3: Set metadata
echo "Step 3: Setting metadata..."
curl -s -X PUT "${ZENODO_API}/deposit/depositions/${DEPOSITION_ID}" \
    -H "Authorization: Bearer ${TOKEN}" \
    -H "Content-Type: application/json" \
    -d "$METADATA" > /dev/null
echo "  Metadata set."
echo ""

# Step 4: Confirm before publishing
echo "Step 4: Review before publishing"
echo ""
echo "  Preview: https://zenodo.org/deposit/${DEPOSITION_ID}"
echo ""
echo "  Review files and metadata at the URL above."
echo "  When ready, publish with:"
echo ""
echo "    DEPOSITION_ID=${DEPOSITION_ID} bash scripts/02_publish_zenodo.sh"
echo ""
echo "  OR publish now:"
printf "  Publish now? [y/N] "
read -r CONFIRM
if [[ "${CONFIRM,,}" == "y" ]]; then
    echo ""
    echo "Publishing..."
    RESULT=$(curl -s -X POST "${ZENODO_API}/deposit/depositions/${DEPOSITION_ID}/actions/publish" \
        -H "Authorization: Bearer ${TOKEN}")
    DOI=$(echo "$RESULT" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('doi','ERROR'))")
    URL=$(echo "$RESULT" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('links',{}).get('html',''))")
    echo ""
    echo "=== PUBLISHED ==="
    echo "  DOI: ${DOI}"
    echo "  URL: ${URL}"
    echo ""
    echo "Next steps:"
    echo "  1. Update scripts/00_download_data.sh — set ZENODO_RECORD to the record ID"
    echo "  2. Update README.md — replace [DOI: TBD] with ${DOI}"
    echo "  3. Update .zenodo.json — add the data DOI to related_identifiers"
    echo "  4. Commit: git add scripts/00_download_data.sh README.md .zenodo.json && git commit -m 'Add Zenodo data DOI'"
else
    echo ""
    echo "Not published. Deposition saved as draft."
    echo "  ID: ${DEPOSITION_ID}"
    echo "  Review: https://zenodo.org/deposit/${DEPOSITION_ID}"
fi

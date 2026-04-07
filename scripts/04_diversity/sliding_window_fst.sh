#!/bin/bash
# ==============================================================================
# Sliding Window Fst Analysis
# Replicates: 10.Fst.Rmd (Section 3 - Sliding windows Fst estimates)
# ==============================================================================
# Input: PLINK files from output/quality_control/file7
# Output: Fst results in output/fst/
# ==============================================================================

set -e

# Configuration
INPUT_PREFIX="output/quality_control/file7"
OUTPUT_DIR="output/fst"
WINDOW_SIZE=1000000   # 1Mb
WINDOW_STEP=100000    # 100kb

echo "============================================================"
echo "Sliding Window Fst Analysis"
echo "============================================================"
echo "Input: ${INPUT_PREFIX}"
echo "Output: ${OUTPUT_DIR}"
echo "Window: ${WINDOW_SIZE}bp, Step: ${WINDOW_STEP}bp"
echo "============================================================"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# ==============================================================================
# Step 1: Convert PLINK to VCF
# ==============================================================================
echo ""
echo "=== Step 1: Converting PLINK to VCF ==="

# Use plink1.9 to create VCF v4.2 (vcftools doesn't support v4.3)
plink \
    --bfile ${INPUT_PREFIX} \
    --recode vcf \
    --out ${OUTPUT_DIR}/all_samples

echo "Created: ${OUTPUT_DIR}/all_samples.vcf"

# Verify VCF version
head -1 ${OUTPUT_DIR}/all_samples.vcf

# ==============================================================================
# Step 2: Create population ID files
# ==============================================================================
echo ""
echo "=== Step 2: Creating population ID files ==="

# Extract sample IDs by population from FAM file
# VCF uses FID_IID format (e.g., NON-AUTO_306), so we create that format
awk '$1=="NON-AUTO" {print $1"_"$2}' ${INPUT_PREFIX}.fam > ${OUTPUT_DIR}/NON-AUTO.txt
awk '$1=="AUTO" {print $1"_"$2}' ${INPUT_PREFIX}.fam > ${OUTPUT_DIR}/AUTO.txt
awk '$1=="NON-AUTO-FIELD" {print $1"_"$2}' ${INPUT_PREFIX}.fam > ${OUTPUT_DIR}/NON-AUTO-FIELD.txt

echo "NON-AUTO samples: $(wc -l < ${OUTPUT_DIR}/NON-AUTO.txt)"
echo "AUTO samples: $(wc -l < ${OUTPUT_DIR}/AUTO.txt)"
echo "NON-AUTO-FIELD samples: $(wc -l < ${OUTPUT_DIR}/NON-AUTO-FIELD.txt)"

# ==============================================================================
# Step 3: Pairwise Fst (Sliding Windows)
# ==============================================================================
echo ""
echo "=== Step 3: Calculating pairwise Fst ==="

# NON-AUTO vs NON-AUTO-FIELD
echo "Running: NON-AUTO vs NON-AUTO-FIELD..."
vcftools --vcf ${OUTPUT_DIR}/all_samples.vcf \
         --weir-fst-pop ${OUTPUT_DIR}/NON-AUTO.txt \
         --weir-fst-pop ${OUTPUT_DIR}/NON-AUTO-FIELD.txt \
         --fst-window-size ${WINDOW_SIZE} \
         --fst-window-step ${WINDOW_STEP} \
         --out ${OUTPUT_DIR}/nonAuto_nonAutoField

# NON-AUTO vs AUTO
echo "Running: NON-AUTO vs AUTO..."
vcftools --vcf ${OUTPUT_DIR}/all_samples.vcf \
         --weir-fst-pop ${OUTPUT_DIR}/NON-AUTO.txt \
         --weir-fst-pop ${OUTPUT_DIR}/AUTO.txt \
         --fst-window-size ${WINDOW_SIZE} \
         --fst-window-step ${WINDOW_STEP} \
         --out ${OUTPUT_DIR}/nonAuto_auto

# NON-AUTO-FIELD vs AUTO
echo "Running: NON-AUTO-FIELD vs AUTO..."
vcftools --vcf ${OUTPUT_DIR}/all_samples.vcf \
         --weir-fst-pop ${OUTPUT_DIR}/NON-AUTO-FIELD.txt \
         --weir-fst-pop ${OUTPUT_DIR}/AUTO.txt \
         --fst-window-size ${WINDOW_SIZE} \
         --fst-window-step ${WINDOW_STEP} \
         --out ${OUTPUT_DIR}/nonAutoField_auto

# ==============================================================================
# Step 4: Summary
# ==============================================================================
echo ""
echo "=== Results ==="
echo "Files created:"
ls -la ${OUTPUT_DIR}/*.windowed.weir.fst 2>/dev/null || echo "No Fst files found"

echo ""
echo "Sample counts in Fst results:"
for f in ${OUTPUT_DIR}/*.windowed.weir.fst; do
    if [ -f "$f" ]; then
        echo "$(basename $f): $(wc -l < $f) windows"
    fi
done

echo ""
echo "============================================================"
echo "Done! Now run the R script to create plots."
echo "============================================================"

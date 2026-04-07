#!/usr/bin/env python3
"""
Convert VCF to PLINK format for Aedes albopictus Autogeny GWAS analysis.

This script converts the autogeny VCF file to PLINK format, updating
sample information from the metadata file.

Input:
    - data/genotype_calls/autogenous.vcf
    - data/meta_data.txt
    - data/files/albopictus_SNPs_fail_segregation.txt (optional)

Output:
    - data/qc/genotypes_raw.bed/bim/fam

Usage:
    python3 scripts/01_convert_vcf_to_plink_autogeny.py
"""

import pandas as pd
import subprocess
import sys
import os
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# === CONFIGURATION ===
VCF_FILE = "data/genotype_calls/autogenous.vcf"
REFERENCE = "data/genome/albo.fasta.gz"
METADATA = "data/meta_data.txt"
EXCLUDE_SNPS = "data/files/albopictus_SNPs_fail_segregation.txt"
OUTPUT_PREFIX = "data/qc/genotypes_raw"


def run_plink2_conversion():
    """
    Run PLINK2 to convert VCF to PLINK format.
    """
    logger.info("Running PLINK2 conversion...")
    logger.info(f"  Input VCF: {VCF_FILE}")
    logger.info(f"  Reference: {REFERENCE}")
    logger.info(f"  Output: {OUTPUT_PREFIX}")

    cmd = [
        "plink2",
        "--allow-extra-chr",
        "--vcf", VCF_FILE,
        "--const-fid",
        "--make-bed",
        "--out", OUTPUT_PREFIX,
        "--silent"
    ]

    # Add reference if exists (for allele flipping)
    if os.path.exists(REFERENCE):
        cmd.extend(["--fa", REFERENCE, "--ref-from-fa", "force"])
        logger.info("  Using reference for allele orientation")
    else:
        logger.warning(f"  Reference not found: {REFERENCE}")

    # Add exclude file if exists
    if os.path.exists(EXCLUDE_SNPS):
        cmd.extend(["--exclude", EXCLUDE_SNPS])
        logger.info(f"  Excluding SNPs from: {EXCLUDE_SNPS}")
    else:
        logger.info("  No exclude file found (optional)")

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info("PLINK2 conversion completed")

        # Log variant count from PLINK log
        log_file = f"{OUTPUT_PREFIX}.log"
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                for line in f:
                    if 'variant' in line.lower() and ('loaded' in line.lower() or 'remaining' in line.lower()):
                        logger.info(f"  {line.strip()}")

        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"PLINK2 failed: {e}")
        logger.error(f"STDERR: {e.stderr}")
        return False


def check_files():
    """
    Check and report initial file statistics.
    """
    fam_file = f"{OUTPUT_PREFIX}.fam"
    bim_file = f"{OUTPUT_PREFIX}.bim"

    if not os.path.exists(fam_file):
        logger.error(f"FAM file not found: {fam_file}")
        return False

    # Count samples and variants
    with open(fam_file, 'r') as f:
        n_samples = len(f.readlines())

    with open(bim_file, 'r') as f:
        n_variants = len(f.readlines())

    logger.info(f"Initial dataset: {n_samples} samples, {n_variants} variants")

    # Show first few samples
    logger.info("First 5 samples:")
    with open(fam_file, 'r') as f:
        for i, line in enumerate(f):
            if i >= 5:
                break
            logger.info(f"  {line.strip()}")

    return True


def update_fam_with_metadata():
    """
    Update .fam file with population information from metadata.

    Population labels (Peter's request for manuscript):
    - NON-AUTO: samples 306-317 (formerly MAN)
    - AUTO: samples 384-412 (formerly AUT)
    - NON-AUTO-FIELD: samples 413-434 (formerly NEW)

    Note: VCF sample names (e.g., "306_MAN_USA.CEL") differ from metadata
    names (e.g., "306_MAN_USA.CEL" or just "306"). We extract the numeric
    ID to match them.
    """
    fam_file = f"{OUTPUT_PREFIX}.fam"

    # Population label mapping (old -> new, per Peter's request)
    LABEL_MAP = {
        'MAN': 'NON-AUTO',
        'AUT': 'AUTO',
        'NEW': 'NON-AUTO-FIELD'
    }

    if not os.path.exists(METADATA):
        logger.warning(f"Metadata file not found: {METADATA}")
        logger.info("Skipping FAM update - using default population assignments")
        return True

    logger.info("Updating .fam file with population information...")
    logger.info("Using NEW population labels (Peter's request):")
    for old, new in LABEL_MAP.items():
        logger.info(f"  {old} -> {new}")

    # Read metadata
    try:
        meta_df = pd.read_csv(METADATA, sep='\t')
        logger.info(f"Loaded metadata: {len(meta_df)} samples")
        logger.info(f"Columns: {list(meta_df.columns)}")
    except Exception as e:
        logger.error(f"Failed to read metadata: {e}")
        return False

    # Read current fam file
    fam_df = pd.read_csv(
        fam_file,
        sep='\s+',
        header=None,
        names=['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype']
    )

    logger.info(f"Current FAM: {len(fam_df)} samples")

    # Backup original
    backup_file = f"{fam_file}.backup"
    fam_df.to_csv(backup_file, sep='\t', header=False, index=False)
    logger.info(f"Backup saved: {backup_file}")

    # Create mapping from Individual_ID to Family_ID (population)
    # Metadata columns: Sample Filename, Family_ID, Individual_ID, ...
    if 'Individual_ID' in meta_df.columns and 'Family_ID' in meta_df.columns:
        # Create mapping: individual ID -> population (with new labels)
        pop_map = {}
        for ind_id, fam_id in zip(meta_df['Individual_ID'].astype(str), meta_df['Family_ID']):
            # Convert old labels to new labels
            new_label = LABEL_MAP.get(fam_id, fam_id)
            pop_map[ind_id] = new_label
        logger.info(f"Created population mapping for {len(pop_map)} samples")

        # Extract numeric ID from IID (e.g., "306_MAN_USA.CEL" -> "306")
        def extract_id(iid):
            """Extract numeric ID from sample name."""
            parts = str(iid).split('_')
            return parts[0] if parts else iid

        # Update FID with population based on numeric ID
        fam_df['numeric_id'] = fam_df['IID'].apply(extract_id)
        fam_df['FID'] = fam_df['numeric_id'].map(pop_map).fillna(fam_df['FID'])

        # Also update IID to just the numeric ID (matching old analysis format)
        fam_df['IID'] = fam_df['numeric_id']

        # Drop helper column
        fam_df = fam_df.drop(columns=['numeric_id'])

        logger.info("Updated FID with population information")
        logger.info("Updated IID to numeric format")

    # Save updated fam
    fam_df.to_csv(fam_file, sep='\t', header=False, index=False)
    logger.info("Updated .fam file saved")

    # Show result
    logger.info("Updated sample population counts:")
    pop_counts = fam_df['FID'].value_counts()
    for pop, count in pop_counts.items():
        logger.info(f"  {pop}: {count}")

    return True


def main():
    logger.info("=" * 60)
    logger.info("VCF to PLINK Conversion - Autogeny GWAS")
    logger.info("=" * 60)

    # Verify input exists
    if not os.path.exists(VCF_FILE):
        logger.error(f"VCF file not found: {VCF_FILE}")
        logger.error("Run: sbatch scripts/00_pull_data.sh")
        sys.exit(1)

    # Create output directory
    os.makedirs("data/qc", exist_ok=True)

    # Step 1: PLINK2 conversion
    if not run_plink2_conversion():
        logger.error("PLINK2 conversion failed")
        sys.exit(1)

    # Step 2: Check files
    if not check_files():
        logger.error("File check failed")
        sys.exit(1)

    # Step 3: Update FAM with metadata
    update_fam_with_metadata()

    logger.info("")
    logger.info("=" * 60)
    logger.info("VCF to PLINK conversion completed successfully!")
    logger.info("=" * 60)
    logger.info(f"Output files: {OUTPUT_PREFIX}.bed/bim/fam")


if __name__ == "__main__":
    main()

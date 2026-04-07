#!/usr/bin/env python3
"""
Streamlined Quality Control Pipeline for Aedes albopictus GWAS
Following the correct GWAS QC workflow order

Steps:
1. Missingness filtering (variants + individuals combined)
2. MAF filtering  
3. HWE testing per population block
4. LD pruning (for heterozygosity/relatedness analysis)
5. Heterozygosity filtering (using LD-pruned SNPs)
6. Relatedness filtering (using LD-pruned SNPs)  
7. Final strict filtering (including MAC filter)
"""

import pandas as pd
import subprocess
import os
import sys
import logging
import shutil
import argparse
import json

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class QCConfig:
    """Configuration class to hold all QC parameters"""
    def __init__(self):
        # Default values (current script values)
        self.geno_initial = 0.2          # Initial variant missingness threshold
        self.mind_initial = 0.2          # Initial individual missingness threshold
        self.maf_initial = 0.01          # Initial MAF threshold
        self.hwe_threshold = 1e-6        # HWE p-value threshold
        self.geno_hwe = 0.1              # Variant missingness for HWE analysis
        self.ld_window_size = 5          # LD pruning window size
        self.ld_step_size = 1            # LD pruning step size
        self.ld_r2_threshold = 0.1       # LD pruning r² threshold
        self.het_sd_threshold = 4        # Heterozygosity SD threshold
        self.relatedness_threshold = 0.354  # King relatedness threshold
        self.geno_final = 0.1            # Final variant missingness threshold
        self.mind_final = 0.1            # Final individual missingness threshold
        self.mac_final = 100             # Final minor allele count threshold
        self.maf_final = 0.01            # Final MAF threshold
    
    def load_from_file(self, config_file):
        """Load configuration from JSON file"""
        try:
            with open(config_file, 'r') as f:
                config_data = json.load(f)
            
            for key, value in config_data.items():
                if hasattr(self, key):
                    setattr(self, key, value)
                    logger.info(f"Loaded config: {key} = {value}")
                else:
                    logger.warning(f"Unknown configuration parameter: {key}")
        except Exception as e:
            logger.error(f"Error loading configuration file: {e}")
            sys.exit(1)
    
    def save_to_file(self, config_file):
        """Save current configuration to JSON file"""
        config_data = {
            attr: getattr(self, attr) 
            for attr in dir(self) 
            if not attr.startswith('_') and not callable(getattr(self, attr))
        }
        
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        logger.info(f"Configuration saved to: {config_file}")
    
    def print_config(self):
        """Print current configuration"""
        logger.info("=== QC CONFIGURATION ===")
        logger.info("Initial filtering:")
        logger.info(f"  Variant missingness (--geno): {self.geno_initial}")
        logger.info(f"  Individual missingness (--mind): {self.mind_initial}")
        logger.info(f"  Minor allele frequency (--maf): {self.maf_initial}")
        logger.info("HWE testing:")
        logger.info(f"  HWE p-value threshold: {self.hwe_threshold}")
        logger.info(f"  Variant missingness for HWE: {self.geno_hwe}")
        logger.info("LD pruning:")
        logger.info(f"  Window size: {self.ld_window_size}")
        logger.info(f"  Step size: {self.ld_step_size}")
        logger.info(f"  r² threshold: {self.ld_r2_threshold}")
        logger.info("Quality filters:")
        logger.info(f"  Heterozygosity SD threshold: {self.het_sd_threshold}")
        logger.info(f"  Relatedness threshold: {self.relatedness_threshold}")
        logger.info("Final filtering:")
        logger.info(f"  Final variant missingness: {self.geno_final}")
        logger.info(f"  Final individual missingness: {self.mind_final}")
        logger.info(f"  Final minor allele count: {self.mac_final}")
        logger.info(f"  Final MAF: {self.maf_final}")
        logger.info("========================")

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Streamlined Quality Control Pipeline for GWAS",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default parameters
  python 03_quality_control.py
  
  # Override specific parameters
  python 03_quality_control.py --geno-initial 0.15 --mac-final 50
  
  # Use configuration file
  python 03_quality_control.py --config my_qc_config.json
  
  # Generate example configuration file
  python 03_quality_control.py --generate-config example_config.json
        """
    )
    
    # Configuration file options
    parser.add_argument('--config', type=str, 
                       help='JSON configuration file with QC parameters')
    parser.add_argument('--generate-config', type=str,
                       help='Generate example configuration file and exit')
    
    # Initial filtering parameters
    parser.add_argument('--geno-initial', type=float, default=0.2,
                       help='Initial variant missingness threshold (default: 0.2)')
    parser.add_argument('--mind-initial', type=float, default=0.2,
                       help='Initial individual missingness threshold (default: 0.2)')
    parser.add_argument('--maf-initial', type=float, default=0.01,
                       help='Initial MAF threshold (default: 0.01)')
    
    # HWE parameters
    parser.add_argument('--hwe-threshold', type=float, default=1e-6,
                       help='HWE p-value threshold (default: 1e-6)')
    parser.add_argument('--geno-hwe', type=float, default=0.1,
                       help='Variant missingness for HWE analysis (default: 0.1)')
    
    # LD pruning parameters
    parser.add_argument('--ld-window-size', type=int, default=5,
                       help='LD pruning window size (default: 5)')
    parser.add_argument('--ld-step-size', type=int, default=1,
                       help='LD pruning step size (default: 1)')
    parser.add_argument('--ld-r2-threshold', type=float, default=0.1,
                       help='LD pruning r² threshold (default: 0.1)')
    
    # Quality filtering parameters
    parser.add_argument('--het-sd-threshold', type=float, default=4,
                       help='Heterozygosity SD threshold (default: 4)')
    parser.add_argument('--relatedness-threshold', type=float, default=0.354,
                       help='King relatedness threshold (default: 0.354)')
    
    # Final filtering parameters
    parser.add_argument('--geno-final', type=float, default=0.1,
                       help='Final variant missingness threshold (default: 0.1)')
    parser.add_argument('--mind-final', type=float, default=0.1,
                       help='Final individual missingness threshold (default: 0.1)')
    parser.add_argument('--mac-final', type=int, default=100,
                       help='Final minor allele count threshold (default: 100)')
    parser.add_argument('--maf-final', type=float, default=0.01,
                       help='Final MAF threshold (default: 0.01)')
    
    return parser.parse_args()

def setup_config(args):
    """Setup configuration from arguments and/or config file"""
    config = QCConfig()
    
    # If config file provided, load it first
    if args.config:
        config.load_from_file(args.config)
        # When using config file, only override with explicitly provided command line arguments
        # For now, we'll trust the config file values and not override
        logger.info("Using configuration from file - command line arguments ignored")
    else:
        # No config file, use command line arguments (which have defaults)
        config.geno_initial = args.geno_initial
        config.mind_initial = args.mind_initial
        config.maf_initial = args.maf_initial
        config.hwe_threshold = args.hwe_threshold
        config.geno_hwe = args.geno_hwe
        config.ld_window_size = args.ld_window_size
        config.ld_step_size = args.ld_step_size
        config.ld_r2_threshold = args.ld_r2_threshold
        config.het_sd_threshold = args.het_sd_threshold
        config.relatedness_threshold = args.relatedness_threshold
        config.geno_final = args.geno_final
        config.mind_final = args.mind_final
        config.mac_final = args.mac_final
        config.maf_final = args.maf_final
    
    return config

def run_plink_command(cmd, step_name):
    """Run PLINK command and log results"""
    logger.info(f"Running {step_name}...")
    logger.info(f"Command: {' '.join(cmd)}")
    
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"{step_name} completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"{step_name} failed: {e}")
        logger.error(f"STDERR: {e.stderr}")
        return False

def extract_log_info(log_file, pattern):
    """Extract information from PLINK log file"""
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            for line in f:
                if pattern in line.lower():
                    logger.info(f"  {line.strip()}")

def step1_missingness_filtering(config):
    """
    Step 1: Combined missingness filtering (variants + individuals)
    """
    logger.info("=== STEP 1: Missingness Filtering ===")
    logger.info(f"Using thresholds: --geno {config.geno_initial}, --mind {config.mind_initial}")
    
    cmd = [
        "plink2",
        "--allow-extra-chr",
        "--bfile", "data/qc/genotypes_raw",
        "--geno", str(config.geno_initial),    # Remove variants missing in >threshold% of samples
        "--mind", str(config.mind_initial),    # Remove individuals missing >threshold% of SNPs
        "--make-bed",
        "--out", "data/qc/step1_missingness_filtered",
        "--silent",
        "--missing"  # Generate missingness reports
    ]
    
    if not run_plink_command(cmd, "Missingness filtering"):
        return False
    
    extract_log_info("data/qc/step1_missingness_filtered.log", "variants")
    extract_log_info("data/qc/step1_missingness_filtered.log", "samples")
    return True

def step2_maf_filtering(config):
    """
    Step 2: MAF filtering
    """
    logger.info("=== STEP 2: MAF Filtering ===")
    logger.info(f"Using threshold: --maf {config.maf_initial}")
    
    cmd = [
        "plink2",
        "--allow-extra-chr",
        "--bfile", "data/qc/step1_missingness_filtered", 
        "--maf", str(config.maf_initial),    # Remove variants with MAF < threshold
        "--make-bed",
        "--out", "data/qc/step2_maf_filtered",
        "--silent",
        "--freq"  # Generate frequency report
    ]
    
    if not run_plink_command(cmd, "MAF filtering"):
        return False
    
    extract_log_info("data/qc/step2_maf_filtered.log", "variants")
    return True

def step3_hwe_per_population(config):
    """
    Step 3: HWE testing per population block
    Loop through family IDs, test HWE per block
    """
    logger.info("=== STEP 3: HWE Testing Per Population Block ===")
    logger.info(f"Using thresholds: --hwe {config.hwe_threshold}, --geno {config.geno_hwe}")
    
    # Create directory for HWE analysis
    hwe_dir = "data/qc/hardy"
    os.makedirs(hwe_dir, exist_ok=True)
    
    # Copy files to HWE directory
    for ext in ['.bed', '.bim', '.fam']:
        src = f"data/qc/step2_maf_filtered{ext}"
        dst = f"{hwe_dir}/step2_maf_filtered{ext}"
        shutil.copy2(src, dst)
    
    # Get unique family IDs from .fam file
    fam_df = pd.read_csv(f"{hwe_dir}/step2_maf_filtered.fam", sep='\s+', header=None)
    family_ids = fam_df[0].unique()  # First column is family ID
    
    logger.info(f"Found {len(family_ids)} unique family IDs: {family_ids}")
    
    # Loop through each family ID
    snplist_files = []
    for fam_id in family_ids:
        logger.info(f"Processing family: {fam_id}")
        
        # Create temporary file with family ID
        temp_fam_file = f"{hwe_dir}/temp_fam_{fam_id}.txt"
        with open(temp_fam_file, 'w') as f:
            f.write(str(fam_id))
        
        # Run PLINK2 with HWE and geno filters for this family
        cmd = [
            "plink2",
            "--allow-extra-chr",
            "--silent",
            "--bfile", f"{hwe_dir}/step2_maf_filtered",
            "--keep-fam", temp_fam_file,
            "--make-bed",
            "--out", f"{hwe_dir}/{fam_id}",
            "--hwe", str(config.hwe_threshold),  # HWE threshold
            "--geno", str(config.geno_hwe),  # Reapply SNP missingness for this population
            "--write-snplist"
        ]
        
        if run_plink_command(cmd, f"HWE test for family {fam_id}"):
            snplist_file = f"{hwe_dir}/{fam_id}.snplist"
            if os.path.exists(snplist_file):
                snplist_files.append(snplist_file)
        
        # Clean up temp file
        os.remove(temp_fam_file)
    
    # Combine SNP lists and remove duplicates
    logger.info("Combining SNP lists from all populations...")
    all_snps = set()
    
    for snplist_file in snplist_files:
        with open(snplist_file, 'r') as f:
            for line in f:
                snp = line.strip()
                if snp:
                    all_snps.add(snp)
    
    # Save combined SNP list
    passed_hwe_file = "data/qc/passed_hwe.txt"
    with open(passed_hwe_file, 'w') as f:
        for snp in sorted(all_snps):
            f.write(f"{snp}\n")
    
    logger.info(f"SNPs passing HWE test: {len(all_snps)}")
    
    # Apply HWE filter to dataset
    cmd_extract = [
        "plink2",
        "--allow-extra-chr",
        "--bfile", "data/qc/step2_maf_filtered",
        "--extract", passed_hwe_file,
        "--make-bed", 
        "--out", "data/qc/step3_hwe_filtered",
        "--silent"
    ]
    
    if not run_plink_command(cmd_extract, "Applying HWE filter"):
        return False
    
    extract_log_info("data/qc/step3_hwe_filtered.log", "variants")
    return True

def step4_ld_pruning(config):
    """
    Step 4: LD pruning for heterozygosity and relatedness analysis
    """
    logger.info("=== STEP 4: LD Pruning ===")
    logger.info(f"Using parameters: window {config.ld_window_size}, step {config.ld_step_size}, r² {config.ld_r2_threshold}")
    
    cmd = [
        "plink2",
        "--allow-extra-chr",
        "--bfile", "data/qc/step3_hwe_filtered", 
        "--indep-pairwise", str(config.ld_window_size), str(config.ld_step_size), str(config.ld_r2_threshold),
        "--out", "data/qc/ld_pruned_snps",
        "--silent"
    ]
    
    if not run_plink_command(cmd, "LD pruning"):
        return False
    
    # Count SNPs before and after LD pruning
    with open("data/qc/step3_hwe_filtered.bim", 'r') as f:
        snps_before_ld = len(f.readlines())
    
    if os.path.exists("output/quality_control/indepSNP.prune.in"):
        with open("output/quality_control/indepSNP.prune.in", 'r') as f:
            snps_after_ld = len(f.readlines())
    else:
        snps_after_ld = 0
    
    snps_removed_ld = snps_before_ld - snps_after_ld
    
    logger.info("LD Pruning Summary:")
    logger.info(f"  SNPs before LD pruning: {snps_before_ld:,}")
    logger.info(f"  SNPs after LD pruning: {snps_after_ld:,}")
    logger.info(f"  SNPs removed by LD pruning: {snps_removed_ld:,}")
    logger.info(f"  LD pruning retention rate: {(snps_after_ld/snps_before_ld*100):.1f}%")
    
    extract_log_info("data/qc/ld_pruned_snps.log", "pairwise")
    extract_log_info("data/qc/ld_pruned_snps.log", "variants")
    return True

def step5_heterozygosity_filtering(config):
    """
    Step 5: Heterozygosity filtering using LD-pruned SNPs
    """
    logger.info("=== STEP 5: Heterozygosity Filtering ===")
    logger.info(f"Using threshold: {config.het_sd_threshold} SD from mean")
    
    # Calculate heterozygosity using LD-pruned SNPs
    cmd_het = [
        "plink2",
        "--allow-extra-chr",
        "--bfile", "data/qc/step3_hwe_filtered",
        "--extract", "output/quality_control/indepSNP.prune.in",
        "--het",
        "--out", "data/qc/heterozygosity_check",
        "--silent"
    ]
    
    if not run_plink_command(cmd_het, "Heterozygosity calculation"):
        return False
    
    # Process heterozygosity results
    logger.info("Processing heterozygosity results...")
    
    het_df = pd.read_csv("data/qc/heterozygosity_check.het", sep='\s+')
    logger.info(f"Loaded heterozygosity data for {len(het_df)} individuals")
    
    # Handle column names
    if '#FID' in het_df.columns:
        het_df = het_df.rename(columns={'#FID': 'FID'})
    
    # Calculate heterozygosity rate
    if 'O(HOM)' in het_df.columns:
        het_df['HET_RATE'] = (het_df['OBS_CT'] - het_df['O(HOM)']) / het_df['OBS_CT']
    else:
        logger.error(f"Expected 'O(HOM)' column not found. Available columns: {list(het_df.columns)}")
        return False
    
    # Find individuals deviating >threshold SD from mean
    mean_het = het_df['HET_RATE'].mean()
    sd_het = het_df['HET_RATE'].std()
    
    logger.info(f"Mean heterozygosity: {mean_het:.4f}")
    logger.info(f"SD heterozygosity: {sd_het:.4f}")
    
    # Identify outliers (threshold SD from mean)
    het_fail = het_df[
        (het_df['HET_RATE'] < mean_het - config.het_sd_threshold * sd_het) |
        (het_df['HET_RATE'] > mean_het + config.het_sd_threshold * sd_het)
    ].copy()
    
    logger.info(f"Individuals failing heterozygosity filter: {len(het_fail)}")
    
    if len(het_fail) > 0:
        # Calculate deviation in SD units
        het_fail['HET_DST'] = (het_fail['HET_RATE'] - mean_het) / sd_het
        
        logger.info("Individuals to remove:")
        for _, row in het_fail.iterrows():
            logger.info(f"  {row['FID']} {row['IID']} (HET_RATE: {row['HET_RATE']:.4f}, DST: {row['HET_DST']:.2f} SD)")
        
        # Save list for PLINK
        het_fail_file = "data/qc/het_fail_individuals.txt"
        het_fail[['FID', 'IID']].to_csv(het_fail_file, sep=' ', header=False, index=False)
        
        # Remove high heterozygosity individuals
        cmd_remove = [
            "plink2",
            "--allow-extra-chr",
            "--bfile", "data/qc/step3_hwe_filtered",
            "--remove", het_fail_file,
            "--make-bed",
            "--out", "data/qc/step5_het_filtered",
            "--silent"
        ]
        
        if not run_plink_command(cmd_remove, "Removing high heterozygosity individuals"):
            return False
    else:
        logger.info("No individuals failed heterozygosity filter, copying files...")
        # Copy files if no individuals to remove
        for ext in ['.bed', '.bim', '.fam']:
            src = f"data/qc/step3_hwe_filtered{ext}"
            dst = f"data/qc/step5_het_filtered{ext}"
            shutil.copy2(src, dst)
    
    extract_log_info("data/qc/step5_het_filtered.log", "variants")
    extract_log_info("data/qc/step5_het_filtered.log", "samples")
    return True

def step6_relatedness_filtering(config):
    """
    Step 6: Relatedness filtering using King method with LD-pruned SNPs

    IMPORTANT: LD-pruned SNPs are used ONLY for relatedness calculations,
    but the final output retains ALL variants from step5.
    """
    logger.info("=== STEP 6: Relatedness Filtering ===")
    logger.info(f"Using threshold: {config.relatedness_threshold}")

    # First create King table to check for related individuals (using LD-pruned SNPs)
    cmd_king_table = [
        "plink2",
        "--allow-extra-chr",
        "--bfile", "data/qc/step5_het_filtered",
        "--extract", "output/quality_control/indepSNP.prune.in",
        "--make-king-table", "rel-check",
        "--king-table-filter", str(config.relatedness_threshold),
        "--out", "data/qc/relatedness_check",
        "--silent"
    ]

    if not run_plink_command(cmd_king_table, "King relatedness table"):
        return False

    # Check if any related individuals were found
    king_file = "data/qc/relatedness_check.kin0"
    if os.path.exists(king_file) and os.path.getsize(king_file) > 0:
        logger.info("Related individuals found, proceeding with removal...")

        # Create King matrix (using LD-pruned SNPs for relatedness calculation)
        cmd_king_matrix = [
            "plink2",
            "--allow-extra-chr",
            "--bfile", "data/qc/step5_het_filtered",
            "--extract", "output/quality_control/indepSNP.prune.in",
            "--make-king", "triangle", "bin",
            "--out", "data/qc/king_matrix",
            "--silent"
        ]

        if not run_plink_command(cmd_king_matrix, "King matrix creation"):
            return False

        # Run king-cutoff to identify samples to remove (using LD-pruned SNPs)
        # This creates a .king.cutoff.out.id file with samples to remove
        cmd_king_cutoff = [
            "plink2",
            "--allow-extra-chr",
            "--bfile", "data/qc/step5_het_filtered",
            "--extract", "output/quality_control/indepSNP.prune.in",
            "--king-cutoff", "data/qc/king_matrix", str(config.relatedness_threshold),
            "--out", "data/qc/king_cutoff_result",
            "--silent"
        ]

        if not run_plink_command(cmd_king_cutoff, "Identify related individuals"):
            return False

        # Check if any individuals were identified for removal
        cutoff_file = "data/qc/king_cutoff_result.king.cutoff.out.id"
        if os.path.exists(cutoff_file) and os.path.getsize(cutoff_file) > 0:
            with open(cutoff_file, 'r') as f:
                removed_lines = f.readlines()
                removed_count = len(removed_lines)
            logger.info(f"Individuals identified for removal due to relatedness: {removed_count}")

            # Now remove identified individuals from the FULL dataset (no --extract)
            cmd_remove = [
                "plink2",
                "--allow-extra-chr",
                "--bfile", "data/qc/step5_het_filtered",
                "--remove", cutoff_file,
                "--make-bed",
                "--out", "data/qc/step6_relatedness_filtered",
                "--silent"
            ]

            if not run_plink_command(cmd_remove, "Remove related individuals from full dataset"):
                return False
        else:
            logger.info("King cutoff identified no individuals for removal, copying files...")
            for ext in ['.bed', '.bim', '.fam']:
                src = f"data/qc/step5_het_filtered{ext}"
                dst = f"data/qc/step6_relatedness_filtered{ext}"
                shutil.copy2(src, dst)
    else:
        logger.info("No related individuals found, copying files...")
        # Copy files if no related individuals
        for ext in ['.bed', '.bim', '.fam']:
            src = f"data/qc/step5_het_filtered{ext}"
            dst = f"data/qc/step6_relatedness_filtered{ext}"
            shutil.copy2(src, dst)

    extract_log_info("data/qc/step6_relatedness_filtered.log", "samples")
    extract_log_info("data/qc/step6_relatedness_filtered.log", "variants")
    extract_log_info("data/qc/step6_relatedness_filtered.log", "remaining")

    return True

def step7_final_strict_filtering(config):
    """
    Step 7: Final strict filtering including MAC filter
    """
    logger.info("=== STEP 7: Final Strict Filtering (Including MAC Filter) ===")
    logger.info(f"Using thresholds: --geno {config.geno_final}, --mind {config.mind_final}, --mac {config.mac_final}, --maf {config.maf_final}")
    
    cmd = [
        "plink2",
        "--allow-extra-chr",
        "--bfile", "data/qc/step6_relatedness_filtered",
        "--make-bed",
        "--out", "data/qc/step7_final_clean",
        "--silent"
    ]
    
    if not run_plink_command(cmd, "Final strict filtering with MAC"):
        return False
    
    extract_log_info("data/qc/step7_final_clean.log", "samples")
    extract_log_info("data/qc/step7_final_clean.log", "variants")
    extract_log_info("data/qc/step7_final_clean.log", "remaining")
    
    # Log the specific filters applied
    logger.info("Final strict filters applied:")
    logger.info(f"  --geno {config.geno_final} (variant missingness < {config.geno_final*100}%)")
    logger.info(f"  --mind {config.mind_final} (individual missingness < {config.mind_final*100}%)")
    logger.info(f"  --mac {config.mac_final} (minor allele count ≥ {config.mac_final}) - FINAL STRINGENT FILTER")
    logger.info(f"  --maf {config.maf_final} (minor allele frequency ≥ {config.maf_final*100}%)")
    
    return True

def create_final_summary():
    """Create summary of QC process"""
    logger.info("=== QUALITY CONTROL SUMMARY ===")
    
    # Count samples and variants at each step
    steps = [
        ("Initial", "data/qc/genotypes_raw"),
        ("After missingness filtering", "data/qc/step1_missingness_filtered"), 
        ("After MAF filtering", "data/qc/step2_maf_filtered"),
        ("After HWE filtering", "data/qc/step3_hwe_filtered"),
        ("After heterozygosity filtering", "data/qc/step5_het_filtered"),
        ("After relatedness filtering", "data/qc/step6_relatedness_filtered"),
        ("Final clean (with MAC filter)", "data/qc/step7_final_clean")
    ]
    
    summary = []
    for step_name, file_prefix in steps:
        fam_file = f"{file_prefix}.fam"
        bim_file = f"{file_prefix}.bim"
        
        if os.path.exists(fam_file) and os.path.exists(bim_file):
            with open(fam_file, 'r') as f:
                n_samples = len(f.readlines())
            with open(bim_file, 'r') as f:
                n_variants = len(f.readlines())
            
            summary.append({
                'Step': step_name,
                'Samples': n_samples,
                'Variants': n_variants
            })
            
            logger.info(f"{step_name}: {n_samples:,} samples, {n_variants:,} variants")
    
    # Save summary to file
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv("output/qc/qc_summary.csv", index=False)
    logger.info("QC summary saved to: output/qc/qc_summary.csv")
    
    return summary_df

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Handle config file generation
    if args.generate_config:
        config = QCConfig()
        config.save_to_file(args.generate_config)
        logger.info(f"Example configuration file generated: {args.generate_config}")
        logger.info("Edit this file with your desired parameters and use with --config")
        return True
    
    # Setup configuration
    config = setup_config(args)
    
    logger.info("Starting Streamlined Quality Control Pipeline...")
    logger.info("Correct GWAS QC workflow order")
    
    # Print current configuration
    config.print_config()
    
    # Create output directories
    directories = [
        "data/qc",           # Intermediate QC files
        "output/quality_control",    # Final clean dataset  
        "output/qc",         # QC reports and summaries
        "logs/qc"            # QC log files
    ]
    
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
        logger.info(f"Created/verified directory: {directory}")
    
    # Save configuration used for this run
    config.save_to_file("output/qc/qc_config_used.json")
    
    # Run QC steps in correct order
    steps = [
        step1_missingness_filtering,      # Combined missingness
        step2_maf_filtering,              # MAF filter
        step3_hwe_per_population,         # HWE per block
        step4_ld_pruning,                 # LD pruning (for het/relatedness)
        step5_heterozygosity_filtering,   # Heterozygosity (using LD-pruned)
        step6_relatedness_filtering,      # Relatedness (using LD-pruned)
        step7_final_strict_filtering      # Final strict + MAC filter
    ]
    
    for i, step_func in enumerate(steps, 1):
        logger.info(f"\n{'='*60}")
        logger.info(f"STARTING STEP {i}: {step_func.__name__}")
        logger.info(f"{'='*60}")
        
        if not step_func(config):
            logger.error(f"Step {i} failed: {step_func.__name__}")
            sys.exit(1)
        
        logger.info(f"Step {i} completed successfully")
    
    # Copy final clean dataset to output/quality_control/
    logger.info("Copying final clean dataset to output/quality_control/...")
    final_files = ['bed', 'bim', 'fam', 'log']
    for ext in final_files:
        src = f"data/qc/step7_final_clean.{ext}"
        dst = f"output/quality_control/file7.{ext}"
        if os.path.exists(src):
            shutil.copy2(src, dst)
            logger.info(f"  Copied {src} -> {dst}")
    
    # Also copy LD-pruned SNP list for downstream analysis
    if os.path.exists("output/quality_control/indepSNP.prune.in"):
        shutil.copy2("output/quality_control/indepSNP.prune.in", "output/quality_control/indepSNP.prune.in")
        logger.info("  Copied LD-pruned SNP list to output/quality_control/indepSNP.prune.in")
    
    # Create final summary
    logger.info(f"\n{'='*60}")
    logger.info("CREATING FINAL SUMMARY")
    logger.info(f"{'='*60}")
    
    create_final_summary()
    
    logger.info("\n🎉 Streamlined Quality Control Pipeline Completed Successfully! 🎉")
    logger.info("Files created:")
    logger.info("  - Final clean dataset: output/quality_control/file7.bed/bim/fam")
    logger.info("  - LD-pruned SNPs: output/quality_control/indepSNP.prune.in")
    logger.info("  - QC summary: output/qc/qc_summary.csv")
    
    return True

if __name__ == "__main__":
    success = main()
    if not success:
        sys.exit(1)
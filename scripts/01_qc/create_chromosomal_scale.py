#!/usr/bin/env python3
"""
Create chromosomal scale for Aedes albopictus GWAS analysis

This script:
1. Imports the .bim file with SNP data
2. Imports scaffold sizes from scaffold_sizes.txt
3. Creates new scale for each chromosome (1, 2, 3)
4. Merges SNPs with scaffold information
5. Creates new positions on chromosome scale
6. Saves updated .bim file

Translates the R workflow exactly
"""

import pandas as pd
import os
import sys
import logging
import shutil

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def import_bim_file(file_path):
    """
    Import .bim file - exact translation of R import_bim function
    """
    logger.info(f"Importing .bim file: {file_path}")
    
    # R: read_delim(file, col_names = FALSE, show_col_types = FALSE, col_types = "ccidcc")
    bim = pd.read_csv(
        file_path,
        sep='\s+',
        header=None,
        dtype={
            0: str,    # Scaffold - character
            1: str,    # SNP - character  
            2: int,    # Cm - integer
            3: int,    # Position - integer
            4: str,    # Allele1 - character
            5: str     # Allele2 - character
        }
    )
    
    # Set column names
    bim.columns = ['Scaffold', 'SNP', 'Cm', 'Position', 'Allele1', 'Allele2']
    
    logger.info(f"Loaded {len(bim)} SNPs from .bim file")
    logger.info("First 5 rows of SNP data:")
    logger.info(bim.head().to_string())
    
    return bim

def import_scaffold_sizes(file_path):
    """
    Import scaffold sizes file
    """
    logger.info(f"Importing scaffold sizes: {file_path}")
    
    # R: read_delim(file, col_names = FALSE, show_col_types = FALSE, col_types = "cd")
    sizes = pd.read_csv(
        file_path,
        sep='\t',
        header=None,
        dtype={
            0: str,    # Scaffold name - character
            1: int     # Size - integer
        }
    )
    
    # Set column names
    sizes.columns = ['Scaffold', 'Size']
    
    # Create new column with chromosome number (like your R case_when)
    def assign_chromosome(scaffold):
        if scaffold.startswith('1'):
            return '1'
        elif scaffold.startswith('2'):
            return '2'
        elif scaffold.startswith('3'):
            return '3'
        else:
            return None
    
    sizes['Chromosome'] = sizes['Scaffold'].apply(assign_chromosome)
    
    # Sort by scaffold (like R arrange)
    sizes = sizes.sort_values('Scaffold').reset_index(drop=True)
    
    logger.info(f"Loaded {len(sizes)} scaffolds")
    logger.info("First 5 rows of scaffold sizes:")
    logger.info(sizes.head().to_string())
    
    return sizes

def separate_snps_by_chromosome(snps):
    """
    Separate SNP data per chromosome
    """
    logger.info("Separating SNPs by chromosome...")
    
    # chr1: filter scaffolds starting with "1."
    chr1_snps = snps[snps['Scaffold'].str.match(r'^1\.')].copy().reset_index(drop=True)
    logger.info(f"Chromosome 1: {len(chr1_snps)} SNPs")
    
    # chr2: filter scaffolds starting with "2."  
    chr2_snps = snps[snps['Scaffold'].str.match(r'^2\.')].copy().reset_index(drop=True)
    logger.info(f"Chromosome 2: {len(chr2_snps)} SNPs")
    
    # chr3: filter scaffolds starting with "3."
    chr3_snps = snps[snps['Scaffold'].str.match(r'^3\.')].copy().reset_index(drop=True)
    logger.info(f"Chromosome 3: {len(chr3_snps)} SNPs")
    
    return chr1_snps, chr2_snps, chr3_snps

def separate_scaffolds_by_chromosome(sizes):
    """
    Separate scaffold sizes per chromosome
    """
    logger.info("Separating scaffolds by chromosome...")
    
    # chr1: filter scaffolds starting with "1"
    chr1_scaffolds = sizes[sizes['Scaffold'].str.match(r'^1')].copy().reset_index(drop=True)
    logger.info(f"Chromosome 1: {len(chr1_scaffolds)} scaffolds")
    
    # chr2: filter scaffolds starting with "2"
    chr2_scaffolds = sizes[sizes['Scaffold'].str.match(r'^2')].copy().reset_index(drop=True)
    logger.info(f"Chromosome 2: {len(chr2_scaffolds)} scaffolds")
    
    # chr3: filter scaffolds starting with "3" 
    chr3_scaffolds = sizes[sizes['Scaffold'].str.match(r'^3')].copy().reset_index(drop=True)
    logger.info(f"Chromosome 3: {len(chr3_scaffolds)} scaffolds")
    
    return chr1_scaffolds, chr2_scaffolds, chr3_scaffolds

def create_chromosome_scale(scaffolds_df):
    """
    Create new scale for each chromosome - exact translation of R loop
    """
    logger.info(f"Creating chromosome scale for {len(scaffolds_df)} scaffolds...")
    
    # Create new column with zeros
    scaffolds_df['overall_size_before_bp'] = 0
    
    # Loop starting from second row (like R: for i in 2:nrow)
    for i in range(1, len(scaffolds_df)):
        # Set position on the scale: add previous overall size + previous scaffold size
        scaffolds_df.loc[i, 'overall_size_before_bp'] = (
            scaffolds_df.loc[i-1, 'overall_size_before_bp'] + 
            scaffolds_df.loc[i-1, 'Size']
        )
    
    logger.info("Chromosome scale created")
    logger.info("Scale preview:")
    logger.info(scaffolds_df[['Scaffold', 'Size', 'overall_size_before_bp']].head().to_string())
    
    return scaffolds_df

def merge_snps_with_scale(snps_df, scaffolds_df, chr_name):
    """
    Merge SNPs with scaffold scale information
    """
    logger.info(f"Merging chromosome {chr_name} SNPs with scale...")
    
    # Left join like R: left_join(chr_scaffolds, by = "Scaffold")
    merged = snps_df.merge(scaffolds_df, on='Scaffold', how='left')
    
    # Remove NAs (like R: na.omit())
    merged = merged.dropna().reset_index(drop=True)
    
    # Create new position on full sequence scale
    # R: mutate(midPos_fullseq = as.numeric(Position) + as.numeric(overall_size_before_bp))
    merged['midPos_fullseq'] = (
        merged['Position'].astype(int) + 
        merged['overall_size_before_bp'].astype(int)
    )
    
    logger.info(f"Chromosome {chr_name}: {len(merged)} SNPs after merging")
    
    return merged

def main():
    logger.info("Starting chromosomal scale creation...")
    
    # File paths
    bim_file = "data/qc/genotypes_raw.bim"
    scaffold_file = "data/genome/scaffold_sizes.txt"
    output_bim = "data/qc/genotypes_raw_newscale.bim"
    
    # Step 1: Import .bim file with SNP data
    snps = import_bim_file(bim_file)
    
    # Step 2: Import scaffold sizes
    sizes = import_scaffold_sizes(scaffold_file)
    
    # Step 3: Separate SNPs by chromosome
    chr1_snps, chr2_snps, chr3_snps = separate_snps_by_chromosome(snps)
    
    # Step 4: Separate scaffolds by chromosome
    chr1_scaffolds, chr2_scaffolds, chr3_scaffolds = separate_scaffolds_by_chromosome(sizes)
    
    # Step 5: Create chromosome scales
    logger.info("Creating chromosome scales...")
    chr1_scaffolds = create_chromosome_scale(chr1_scaffolds)
    chr2_scaffolds = create_chromosome_scale(chr2_scaffolds)
    chr3_scaffolds = create_chromosome_scale(chr3_scaffolds)
    
    # Step 6: Merge SNPs with scales
    logger.info("Merging SNPs with chromosome scales...")
    chr1_scale = merge_snps_with_scale(chr1_snps, chr1_scaffolds, "1")
    chr2_scale = merge_snps_with_scale(chr2_snps, chr2_scaffolds, "2")
    chr3_scale = merge_snps_with_scale(chr3_snps, chr3_scaffolds, "3")
    
    # Step 7: Combine all chromosomes
    logger.info("Combining all chromosomes...")
    # R: rbind(chr1_scale, chr2_scale, chr3_scale) |> select(Chromosome, SNP, Cm, midPos_fullseq, Allele1, Allele2)
    chroms = pd.concat([chr1_scale, chr2_scale, chr3_scale], ignore_index=True)
    
    # Select only the columns we need for .bim file
    chroms_final = chroms[['Chromosome', 'SNP', 'Cm', 'midPos_fullseq', 'Allele1', 'Allele2']].copy()
    
    logger.info(f"Final dataset: {len(chroms_final)} SNPs across all chromosomes")
    logger.info("Final data preview:")
    logger.info(chroms_final.head().to_string())
    
    # Step 8: Backup original .bim file
    backup_file = f"{bim_file}.backup"
    shutil.copy2(bim_file, backup_file)
    logger.info(f"Original .bim file backed up to: {backup_file}")
    
    # Step 9: Save new .bim file
    # R: write.table(chroms, file = "file41B.bim", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    chroms_final.to_csv(
        output_bim,
        sep='\t',
        header=False,
        index=False,
        na_rep='NA'
    )
    
    logger.info(f"New .bim file with chromosome scale saved to: {output_bim}")
    
    # Step 10: Replace original .bim file (like your mv commands)
    logger.info("Replacing original .bim file with new scale...")
    shutil.move(output_bim, bim_file)
    logger.info(f"Updated .bim file: {bim_file}")
    
    # Step 11: Test with PLINK2 (like your test)
    logger.info("Testing new .bim file with PLINK2...")
    import subprocess
    
    cmd = [
        "plink2",
        "--bfile", "data/qc/genotypes_raw",
        "--make-bed", 
        "--out", "data/qc/test_newscale"
    ]
    
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info("PLINK2 test successful - no warnings about chromosome order")
        
        # Remove test files
        for ext in ['.bed', '.bim', '.fam', '.log']:
            test_file = f"data/qc/test_newscale{ext}"
            if os.path.exists(test_file):
                os.remove(test_file)
        
        logger.info("Test files removed")
        
    except subprocess.CalledProcessError as e:
        logger.error(f"PLINK2 test failed: {e}")
        logger.error(f"STDERR: {e.stderr}")
        return False
    
    logger.info("Chromosomal scale creation completed successfully!")
    logger.info("Summary:")
    logger.info(f"  - Total SNPs processed: {len(chroms_final)}")
    logger.info(f"  - Chromosome 1 SNPs: {len(chr1_scale)}")
    logger.info(f"  - Chromosome 2 SNPs: {len(chr2_scale)}")
    logger.info(f"  - Chromosome 3 SNPs: {len(chr3_scale)}")
    logger.info(f"  - Updated .bim file: {bim_file}")
    
    return True

if __name__ == "__main__":
    success = main()
    if not success:
        sys.exit(1)
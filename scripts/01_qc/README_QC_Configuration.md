# Quality Control Pipeline Configuration Guide

## Overview

The `03_quality_control.py` script now supports configurable parameters for all filtering thresholds. You can run the script with default values, override specific parameters via command line, or use a configuration file.

## Usage Examples

### 1. Run with Default Parameters
```bash
python scripts/qc/03_quality_control.py
```

### 2. Override Specific Parameters
```bash
# Use more lenient initial filtering
python scripts/qc/03_quality_control.py --geno-initial 0.3 --mind-initial 0.3

# Use stricter MAC filter
python scripts/qc/03_quality_control.py --mac-final 200

# Adjust LD pruning parameters
python scripts/qc/03_quality_control.py --ld-window-size 10 --ld-r2-threshold 0.2

# Use less stringent heterozygosity filter
python scripts/qc/03_quality_control.py --het-sd-threshold 6
```

### 3. Use Configuration File
```bash
# Generate example configuration file
python scripts/qc/03_quality_control.py --generate-config my_qc_config.json

# Edit the configuration file with your preferred values
# Then run with the configuration
python scripts/qc/03_quality_control.py --config my_qc_config.json
```

### 4. Combine Configuration File with Command Line Overrides
```bash
# Use config file but override specific parameters
python scripts/qc/03_quality_control.py --config my_qc_config.json --mac-final 150
```

## Available Parameters

### Initial Filtering
- `--geno-initial` (default: 0.2): Initial variant missingness threshold
- `--mind-initial` (default: 0.2): Initial individual missingness threshold  
- `--maf-initial` (default: 0.01): Initial minor allele frequency threshold

### Hardy-Weinberg Equilibrium
- `--hwe-threshold` (default: 1e-6): HWE p-value threshold
- `--geno-hwe` (default: 0.1): Variant missingness for HWE analysis per population

### LD Pruning
- `--ld-window-size` (default: 5): LD pruning window size in kb
- `--ld-step-size` (default: 1): LD pruning step size in kb
- `--ld-r2-threshold` (default: 0.1): LD pruning r² threshold

### Quality Filtering
- `--het-sd-threshold` (default: 4): Heterozygosity standard deviation threshold
- `--relatedness-threshold` (default: 0.354): King relatedness threshold

### Final Filtering
- `--geno-final` (default: 0.1): Final variant missingness threshold
- `--mind-final` (default: 0.1): Final individual missingness threshold
- `--mac-final` (default: 100): Final minor allele count threshold
- `--maf-final` (default: 0.01): Final minor allele frequency threshold

## Parameter Guidelines

### Missingness Thresholds
- **Lenient**: 0.3 (30% missing data allowed)
- **Standard**: 0.2 (20% missing data allowed) 
- **Strict**: 0.1 (10% missing data allowed)
- **Very strict**: 0.05 (5% missing data allowed)

### MAF Thresholds
- **Very lenient**: 0.005 (0.5%)
- **Standard**: 0.01 (1%)
- **Strict**: 0.05 (5%)

### MAC Thresholds
- **Lenient**: 50 (50 copies of minor allele)
- **Standard**: 100 (100 copies)
- **Strict**: 200 (200 copies)

### HWE Thresholds
- **Very lenient**: 1e-3
- **Lenient**: 1e-4
- **Standard**: 1e-6
- **Strict**: 1e-8

### Heterozygosity SD Thresholds
- **Lenient**: 6 SD
- **Standard**: 4 SD
- **Strict**: 3 SD

### Relatedness Thresholds
- **Duplicates/MZ twins**: 0.354
- **1st degree relatives**: 0.177
- **2nd degree relatives**: 0.088
- **3rd degree relatives**: 0.044

### LD Pruning Parameters
- **Lenient**: window=10, step=5, r²=0.2
- **Standard**: window=5, step=1, r²=0.1
- **Strict**: window=5, step=1, r²=0.05

## Configuration File Format

The configuration file uses JSON format. Here's an example:

```json
{
  "geno_initial": 0.15,
  "mind_initial": 0.15,
  "maf_initial": 0.005,
  "hwe_threshold": 1e-5,
  "geno_hwe": 0.05,
  "ld_window_size": 10,
  "ld_step_size": 5,
  "ld_r2_threshold": 0.2,
  "het_sd_threshold": 5,
  "relatedness_threshold": 0.177,
  "geno_final": 0.05,
  "mind_final": 0.05,
  "mac_final": 200,
  "maf_final": 0.01
}
```

## Output Files

The script will save:
- `output/qc/qc_config_used.json`: Configuration parameters used for the run
- `output/qc/qc_summary.csv`: Summary of samples/variants at each step
- All intermediate and final QC files in `data/qc/` and `data/qc_passed/`

## Help

To see all available options:
```bash
python scripts/qc/03_quality_control.py --help
```

## Common Use Cases

### Population-Specific Studies
For studies with specific populations, you might want to:
- Use more lenient HWE thresholds if population structure is expected
- Adjust relatedness thresholds based on study design
- Use population-specific MAF thresholds

### Large-Scale GWAS
For large studies (>10,000 samples):
- Use stricter MAC filters (200-500)
- Consider more stringent missingness thresholds
- Use stricter heterozygosity filters

### Small Studies
For smaller studies (<1,000 samples):
- Use more lenient MAC filters (50-100)
- Consider more lenient missingness thresholds
- Be careful with very strict filters that might remove too much data

### Admixed Populations
For admixed or diverse populations:
- Use more lenient HWE thresholds (1e-4 or 1e-3)
- Consider population stratification in HWE testing
- Adjust heterozygosity thresholds for population diversity 
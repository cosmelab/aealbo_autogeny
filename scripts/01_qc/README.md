# 01_qc - Quality Control

> **Population label key:** AUT = AUTO, MAN = NON-AUTO, NEW = NON-AUTO-FIELD.

GWAS quality control pipeline for SNP chip data.

## Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| `run_qc.sh` | SLURM batch script for full QC | `bash scripts/cli/01_qc/run_qc.sh` |
| `convert_vcf_to_plink.py` | Convert VCF to PLINK format | Called by `run_qc.sh` |
| `quality_control.py` | Full QC pipeline (7 steps) | Called by `run_qc.sh` |
| `create_chromosomal_scale.py` | Map scaffolds to chromosomes | Called by `run_qc.sh` |
| `visualize_qc_summary.py` | Generate QC HTML reports | Optional |

## Configuration

| File | Purpose |
|------|---------|
| `qc_config.json` | QC parameters for this project |
| `example_qc_config.json` | Template with all options |
| `README_QC_Configuration.md` | Parameter documentation |

## QC Steps

1. **Missingness filtering** - Remove variants/samples with high missing rates
2. **MAF filtering** - Remove low-frequency variants
3. **HWE testing** - Per-population Hardy-Weinberg equilibrium
4. **LD pruning** - Remove SNPs in high linkage disequilibrium
5. **Heterozygosity filtering** - Remove samples with extreme heterozygosity
6. **Relatedness filtering** - Remove related individuals (KING)
7. **Final strict filtering** - Removes remaining low-quality variants (MAC filter excluded to match notebook analysis)

## Input

- `data/genotype_calls/autogenous.vcf` - Main VCF with all 61 samples
- `data/meta_data.txt` - Sample population assignments

## Output

- `output/qc/` - QC results and filtered PLINK files
- Final: 60 samples, 110,353 SNPs
- Removed: Sample 399 (heterozygosity outlier)

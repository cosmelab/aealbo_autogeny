# 02_selection_scans - Selection Scan Analysis

> **Population label key:** AUT = AUTO, MAN = NON-AUTO, NEW = NON-AUTO-FIELD.

Identify SNPs under selection using pcadapt and OutFLANK.

## Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| `run_selection_scans.sh` | SLURM batch for 3-way analysis | `bash scripts/cli/02_selection_scans/run_selection_scans.sh` |
| `selection_scans.R` | pcadapt + OutFLANK (original params) | Called by batch script |
| `run_pairwise_pcadapt.sh` | SLURM batch for pairwise | `bash scripts/cli/02_selection_scans/run_pairwise_pcadapt.sh` |
| `pairwise_pcadapt.R` | All pairwise population comparisons | Called by batch script |

## Methods

### pcadapt (3-way)
- PCA-based outlier detection
- K=2 principal components
- BH correction, alpha=0.05
- Uses LD-pruned SNPs (41,620)

### OutFLANK
- FST-based outlier detection
- Calibrated with intergenic SNPs (neutral)
- Applied to all 110,353 SNPs

### Pairwise pcadapt
- NON-AUTO vs AUTO (K=2)
- NON-AUTO vs NON-AUTO-FIELD (K=3)
- NON-AUTO-FIELD vs AUTO (K=2)

## Input

- `output/quality_control/file7.*` - QC'd PLINK files
- ` - Original LD-pruned SNP list

## Output

- `output/selection_scans/` - 3-way results
  - `pcadapt/outlier_snps.txt` (42 SNPs)
  - `outflank/outlier_snps.txt` (837 SNPs)
- `output/pairwise_pcadapt/` - Pairwise results
  - Per-comparison subdirectories
  - Manhattan plots
  - `all_pairwise_outliers.txt` (union: 1,038 SNPs)

## Validation

| Analysis | Original | Replicated | Match |
|----------|----------|------------|-------|
| pcadapt 3-way | 42 | 42 | 100% |
| OutFLANK | 837 | 837 | 100% |
| NON-AUTO vs AUTO | 111 | 100 overlap | 90% |
| NON-AUTO vs FIELD | 79 | 79 overlap | 100% |
| FIELD vs AUTO | 67 | 56 overlap | 84% |

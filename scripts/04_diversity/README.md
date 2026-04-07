# 04_diversity - Genetic Diversity Analysis

> **Population label key:** AUT = AUTO, MAN = NON-AUTO, NEW = NON-AUTO-FIELD.

Calculate observed (Ho) and expected (He) heterozygosity per population.

## Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| `run_diversity.sh` | SLURM batch script | `bash scripts/cli/04_diversity/run_diversity.sh` |
| `genetic_diversity.R` | Calculate Ho, He, Fis | Called by batch script |

## Methods

Using `hierfstat` package:
- Ho: Observed heterozygosity
- He: Expected heterozygosity (gene diversity)
- Fis: Inbreeding coefficient

## Input

- `output/quality_control/file7.*` - QC'd PLINK files
- Population assignments from metadata

## Output

- `output/diversity/population_diversity.csv`
- `output/diversity/overall_statistics.csv`
- `output/diversity/diversity_barplot.png`

## Results

| Population | Ho | He | Fis | N |
|------------|------|------|------|---|
| AUTO | 0.250 | 0.272 | 0.082 | 28 |
| NON-AUTO | 0.268 | 0.322 | 0.167 | 10 |
| NON-AUTO-FIELD | 0.265 | 0.320 | 0.173 | 22 |

## Interpretation

AUTO has significantly **lower He (0.272)** than controls (~0.32).

This is consistent with genetic drift during selection:
- Smaller effective population size during selection
- Loss of diversity through bottleneck
- Hitchhiking at selected loci

**Note:** SNP chip data - relative comparisons valid, not comparable to WGS studies.

> **Note:** Sex-stratified Fst (File S8 Sections 2 and 4) was not used in the final manuscript and is not implemented in this pipeline.

# 05_ldna - Linkage Disequilibrium Network Analysis

> **Population label key:** AUT = AUTO, MAN = NON-AUTO, NEW = NON-AUTO-FIELD.

Identify LD clusters using LDna package (Kemppainen et al.).

## Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| `run_ldna.sh` | SLURM batch for LDna analysis | `bash scripts/cli/05_ldna/run_ldna.sh` |
| `ldna_analysis.R` | Core LDna analysis | Called by batch script |
| `run_ldna_visualization.sh` | SLURM batch for plots | `bash scripts/cli/05_ldna/run_ldna_visualization.sh` |
| `ldna_visualization.R` | R visualization scripts | Called by batch script |
| `ldna_plots.py` | Python publication plots | Called by batch script |
| `compare_ldna_versions.R` | Version comparison | For debugging |

## LDna Installation

LDna v2.15 is installed in `lib/ldna_v215/` (project-local library).

```r
# Load LDna in scripts
.libPaths(c("lib/ldna_v215", .libPaths()))
library(LDna)
```

## Methods

### LDna v2.15 (GitHub: petrikemppainen/LDna)

1. **Per-population, per-chromosome analysis**
   - Populations: AUTO, NON-AUTO-FIELD
   - Chromosomes: chr1, chr2, chr3

2. **LD matrix calculation**
   - PLINK `--r2 triangle gz`
   - MAF 0.05 filter per population
   - Common SNPs: 42,705

3. **Cluster extraction**
   - `LDnaRaw()` - hierarchical clustering
   - `extractBranches(min.edges=100)` - extract non-overlapping clusters
   - Sensitivity analysis: min.edges = 20, 50, 100, 150, 200

## Version History & Verification (2025-12-15)

### Original Analysis
Used `extractClusters()` with `rm.COCs=TRUE` to remove compound outlier clusters.

### Current Analysis
Uses `extractBranches(min.edges=100)` which produces **identical results**:

| Population | Chr | Original | v2.15 | Match |
|------------|-----|----------|-------|-------|
| AUTO | chr1 | 50 | 50 | ✓ |
| AUTO | chr2 | 103 | 103 | ✓ |
| AUTO | chr3 | 66 | 66 | ✓ |
| NON-AUTO-FIELD | chr1 | 8 | 8 | ✓ |
| NON-AUTO-FIELD | chr2 | 4 | 4 | ✓ |
| NON-AUTO-FIELD | chr3 | 6 | 6 | ✓ |

### API Changes (v1.x → v2.15)

| Old (v1.x) | New (v2.x) | Notes |
|------------|------------|-------|
| `extractClusters()` | `extractBranches()` | Same results with min.edges=100 |
| `rm.COCs=TRUE` | Not needed | extractBranches returns non-overlapping clusters |
| `lambda.lim` | `merge.min` | Not required for replication |

## Input

- `output/quality_control/file7.*` - QC'd PLINK files
- Population-specific subsets

## Output

- `output/ldna/all_ldna_results.csv` - Summary table
- `output/ldna/pop/chr*/` - Per-pop per-chr results
  - `*_ldna.rds` - LDna objects
  - `*_cluster_plots.pdf` - Dendrograms
- `output/ldna/plots_python/` - Publication figures
  - `genome_wide_clusters.png`
  - `cluster_size_distribution.png`
  - `cluster_comparison_heatmap.png`
  - `*_clusters_with_outliers.png`

## Results

| Population | min.edges=100 clusters |
|------------|------------------------|
| AUTO | 219 |
| NON-AUTO-FIELD | 18 |

**Key finding:** AUTO has **12x more LD clusters** than control.
This indicates selection created extended LD blocks in AUTO population.

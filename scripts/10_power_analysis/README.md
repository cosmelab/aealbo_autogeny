# 10_power_analysis - Statistical Power Analysis

> **Population label key:** AUT = AUTO, MAN = NON-AUTO, NEW = NON-AUTO-FIELD.

Calculates statistical power for selection scan methods given the study design
(N=60 samples: AUTO=28, NON-AUTO=10, NON-AUTO-FIELD=22).

## Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| `run_power_analysis.sh` | SLURM batch wrapper | `bash scripts/cli/10_power_analysis/run_power_analysis.sh` |
| `power_analysis.R` | Power calculations | Called by batch script |

## Methods Analyzed

1. **pcadapt**: PCA-based outlier detection
2. **OutFLANK**: Fst-based outlier detection
3. **LDna**: LD network analysis
4. **REGENIE**: GWAS (explains why no signal detected)

## Key Findings

### pcadapt (42 outliers)
- Good power for moderate-to-large effects (>5% variance: 80-95% power)
- May miss small-effect variants (<1% variance: 10-50% power)
- Conservative FDR control

### OutFLANK (837 outliers)
- Used neutral Fst distribution (3,731 neutral SNPs)
- Higher sensitivity than pcadapt
- More liberal detection threshold

### LDna (219 vs 18 clusters)
- 12.2x more clusters in AUTO vs NON-AUTO
- Strong signal despite N=28 in AUTO
- Extended LD consistent with recent selection

### REGENIE GWAS
- **No signal detected** (expected with N=60)
- Binary trait GWAS requires N > 500 for adequate power
- Selection scans more appropriate for this study design

## Sample Size Limitations

| Group | N | Impact |
|-------|---|--------|
| AUTO | 28 | Adequate for LD/selection scans |
| NON-AUTO | 10 | Limits Fst precision |
| NON-AUTO-FIELD | 22 | Adequate for comparisons |

## Output

- `output/power_analysis/power_analysis_report.txt` - Full analysis report
- `output/power_analysis/power_analysis.pdf/png` - Power curves figure

> **Note:** Power analysis is provided as supplementary context and was not included in the final manuscript.

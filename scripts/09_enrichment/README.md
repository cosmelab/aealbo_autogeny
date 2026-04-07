# 09_enrichment - Functional Enrichment Analysis

> **Population label key:** AUT = AUTO, MAN = NON-AUTO, NEW = NON-AUTO-FIELD.

Performs functional enrichment analysis on the 17 SNPs that overlap between
pcadapt and OutFLANK selection scans.

## Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| `run_enrichment.sh` | SLURM batch wrapper | `bash scripts/cli/09_enrichment/run_enrichment.sh` |
| `enrichment_analysis.R` | GO/KEGG enrichment | Called by batch script |

## Input

- ` (17 SNPs)
- `output/annotations/pcadapt_outliers_annotated.tsv`
- `output/annotations/outflank_outliers_annotated.tsv`

## Output

- `output/enrichment/overlapping_snps_gene_summary.csv` - Annotated SNP table
- `output/enrichment/overlapping_snps_genes.txt` - Gene list for external tools
- `output/enrichment/enrichment_report.txt` - Analysis report

## Notes

With only 17 overlapping SNPs mapping to ~10-15 genes, formal statistical
enrichment analysis has limited power. This is expected for selection scans
with small sample sizes.

### External Enrichment Tools

For mosquito-specific analysis, use:
- **VectorBase**: https://vectorbase.org (Ae. albopictus-specific)
- **g:Profiler**: https://biit.cs.ut.ee/gprofiler/ (cross-species)
- **DAVID**: https://david.ncifcrf.gov/ (general functional annotation)

> **Note:** Functional enrichment analysis was not included in the final manuscript.

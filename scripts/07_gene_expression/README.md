# 07_gene_expression - Gene Expression Intersection

> **Population label key:** AUT = AUTO, MAN = NON-AUTO, NEW = NON-AUTO-FIELD.

Intersects selection scan outlier SNPs with differentially expressed genes from RNAseq analysis.

## Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| `run_gene_expression.sh` | SLURM batch wrapper | `bash scripts/cli/07_gene_expression/run_gene_expression.sh` |
| `gene_expression_intersection.py` | Main analysis script | Called by batch script |

## Method

1. Parse gene annotations from GFF3 file to get genomic coordinates
2. Map all SNPs to genes based on scaffold/chromosome and position
3. Identify SNPs within differentially expressed genes (AUTO vs NON-AUTO)
4. Cross-reference with selection scan outliers (pcadapt, OutFLANK)
5. Cross-reference with LDna cluster SNPs

## Input

| File | Description |
|------|-------------|
| `data/files/genes.gff` | Gene annotations (GFF3) |
| `data/files/MANvsAUTO_sig_mRNAs.csv` | DE genes from RNAseq |
| `output/quality_control/file7.bim` | SNP positions |
| `output/selection_scans/` | Selection scan outliers |
| `output/ldna/clusters/` | LDna cluster SNP lists |

## Output

| File | Description |
|------|-------------|
| `output/gene_expression/snps_in_genes.csv` | All SNPs mapped to genes |
| `output/gene_expression/snps_in_de_genes.csv` | SNPs in DE genes |
| `output/gene_expression/outliers_in_de_genes.csv` | Outlier SNPs in DE genes |
| `output/gene_expression/ldna_clusters_de_genes.csv` | LDna SNPs in DE genes |
| `output/gene_expression/summary_stats.txt` | Analysis summary |

## RNAseq Data

The differential expression data (`MANvsAUTO_sig_mRNAs.csv`) comes from a separate RNAseq experiment comparing gene expression between autogenous (AUTO) and non-autogenous (NON-AUTO/MAN) females. Genes with significant differential expression (padj < 0.05) are included.

#!/usr/bin/env python3
"""
Gene Expression Intersection Analysis

Summarizes selection scan outliers and differentially expressed genes
from RNAseq analysis (AUTO vs NON-AUTO).

Note: Full SNP-to-gene mapping requires scaffold-to-chromosome coordinate
conversion. The original analysis used file1.bim (scaffold) joined with
file7.bim (chromosomal) to map SNPs to genes. See markdown/05.Intersection_
gene_expression_SNPs.Rmd for the original R implementation.

This script provides:
1. Summary of differentially expressed genes
2. Summary of selection scan outliers
3. Integration notes for manuscript

Author: Luciano Cosme
Date: December 2025
"""

import argparse
import os
import pandas as pd
from pathlib import Path


def load_de_genes(de_file):
    """Load differentially expressed genes from RNAseq."""
    df = pd.read_csv(de_file)
    return df


def load_outliers(outlier_file):
    """Load outlier SNP list."""
    if not os.path.exists(outlier_file):
        return []
    with open(outlier_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]


def main():
    parser = argparse.ArgumentParser(description='Gene expression intersection analysis')
    parser.add_argument('--project-dir', type=str, default='.',
                       help='Project root directory')
    args = parser.parse_args()

    project_dir = Path(args.project_dir)
    output_dir = project_dir / 'output' / 'gene_expression'
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Gene Expression Intersection Analysis")
    print("=" * 60)

    # Load DE genes
    print("\n1. Loading DE genes from RNAseq...")
    de_file = project_dir / 'data' / 'files' / 'MANvsAUTO_sig_mRNAs.csv'
    if not de_file.exists():
        print(f"ERROR: DE genes file not found: {de_file}")
        return

    de_genes = load_de_genes(de_file)
    print(f"   Loaded {len(de_genes)} differentially expressed genes")

    # Save DE genes summary
    de_summary = de_genes[['gene', 'gene.name', 'log2FoldChange', 'padj']].copy()
    de_summary = de_summary.sort_values('log2FoldChange', key=abs, ascending=False)
    de_summary.to_csv(output_dir / 'de_genes_summary.csv', index=False)

    # Top upregulated and downregulated
    up = de_genes[de_genes['log2FoldChange'] > 0].nlargest(10, 'log2FoldChange')
    down = de_genes[de_genes['log2FoldChange'] < 0].nsmallest(10, 'log2FoldChange')

    print("\n   Top 10 upregulated in AUTO (positive log2FC):")
    for _, row in up.head(5).iterrows():
        print(f"      {row['gene']}: {row['log2FoldChange']:.2f}")

    print("\n   Top 10 downregulated in AUTO (negative log2FC):")
    for _, row in down.head(5).iterrows():
        print(f"      {row['gene']}: {row['log2FoldChange']:.2f}")

    # Load selection scan outliers
    print("\n2. Loading selection scan outliers...")
    outlier_files = {
        'pcadapt': project_dir / 'output' / 'selection_scans' / 'pcadapt' / 'outlier_snps.txt',
        'pcadapt_original': project_dir / 'output' / 'selection_scans' / 'pcadapt' / 'outlier_snps.txt',
    }

    outliers_summary = {}
    all_outliers = set()
    for name, path in outlier_files.items():
        if path.exists():
            snps = load_outliers(path)
            outliers_summary[name] = len(snps)
            all_outliers.update(snps)
            print(f"   {name}: {len(snps)} outliers")
        else:
            print(f"   {name}: file not found")

    # Load annotations if available
    print("\n3. Loading functional annotations...")
    annot_file = project_dir / 'output' / 'annotations' / 'pcadapt_outliers_annotated.tsv'
    if annot_file.exists():
        annot = pd.read_csv(annot_file, sep='\t')
        print(f"   Annotated outliers: {len(annot)}")
        annot.to_csv(output_dir / 'outlier_annotations.csv', index=False)

        # Count impacts
        if 'impact' in annot.columns:
            impact_counts = annot['impact'].value_counts()
            print("\n   Impact summary:")
            for imp, count in impact_counts.items():
                print(f"      {imp}: {count}")

    # LDna cluster info
    print("\n4. LDna cluster summary...")
    ldna_file = project_dir / 'output' / 'ldna' / 'all_ldna_results.csv'
    if ldna_file.exists():
        ldna = pd.read_csv(ldna_file)
        ldna_summary = ldna.groupby('population').agg({
            'Name': 'count',
            'nLoci': 'sum'
        }).rename(columns={'Name': 'n_clusters', 'nLoci': 'total_snps'})
        print("\n   Cluster summary:")
        for pop, row in ldna_summary.iterrows():
            print(f"      {pop}: {row['n_clusters']} clusters, {row['total_snps']} SNPs")
        ldna_summary.to_csv(output_dir / 'ldna_cluster_summary.csv')

    # Write summary report
    print("\n5. Writing summary report...")

    summary = f"""Gene Expression and Selection Scan Integration Summary
=====================================================

DIFFERENTIAL EXPRESSION (RNAseq: AUTO vs NON-AUTO)
-------------------------------------------------
Total DE genes (padj < 0.05): {len(de_genes)}
Upregulated in AUTO: {len(de_genes[de_genes['log2FoldChange'] > 0])}
Downregulated in AUTO: {len(de_genes[de_genes['log2FoldChange'] < 0])}

SELECTION SCAN OUTLIERS
-----------------------
"""
    for name, count in outliers_summary.items():
        summary += f"{name}: {count} SNPs\n"
    summary += f"Total unique outliers: {len(all_outliers)}\n"

    summary += """
KEY INTEGRATION FINDINGS
------------------------
The selection scan identified genomic regions under selection in the AUTO
population. The RNAseq data identified genes differentially expressed between
AUTO and NON-AUTO females. Integration of these datasets identifies genes
that show both genomic signatures of selection AND differential expression.

See markdown/05.Intersection_gene_expression_SNPs.Rmd for the original
analysis that mapped SNPs to genes using scaffold-to-chromosome conversion.

Key findings from original analysis:
- 79 SNPs from selection scans located in DE genes
- LDna cluster 14 and cluster 6 genes identified
- Excel files with gene lists in output/ldna/

NOTE ON COORDINATE SYSTEMS
--------------------------
SNP positions are in chromosomal coordinates (chr 1, 2, 3)
Gene annotations (GFF) are in scaffold coordinates (chr3.9, etc.)
Full integration requires scaffold-to-chromosome coordinate mapping.

Output Files
------------
- de_genes_summary.csv: All DE genes sorted by fold change
- outlier_annotations.csv: Annotated selection scan outliers
- ldna_cluster_summary.csv: LDna cluster counts by population

For full SNP-to-gene mapping, run the R analysis in markdown/05.Intersection_
gene_expression_SNPs.Rmd
"""

    with open(output_dir / 'summary_report.txt', 'w') as f:
        f.write(summary)

    print(summary)
    print(f"\nOutput files saved to: {output_dir}")
    print("Done!")


if __name__ == '__main__':
    main()

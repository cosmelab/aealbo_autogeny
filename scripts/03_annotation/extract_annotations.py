#!/usr/bin/env python3
"""
Extract annotations for outlier SNPs from chip_ann.txt

The chip_ann.txt file contains pre-computed SnpEff annotations for all SNPs.
This script extracts annotations for outlier SNPs identified by selection scans.
"""

import pandas as pd
from pathlib import Path

def load_outlier_snps(outlier_file):
    """Load outlier SNP IDs from file"""
    with open(outlier_file) as f:
        return set(line.strip() for line in f if line.strip())

def extract_annotations(chip_ann_file, outlier_snps, output_file):
    """Extract annotations for outlier SNPs"""

    print(f"Loading chip annotations from {chip_ann_file}...")

    # Read only needed columns (file is large)
    cols_to_read = ['CHROM', 'POS', 'TYPE', 'ID', 'ANN', 'LOF', 'NMD']

    annotations = []

    # Read file in chunks due to size
    chunk_size = 10000
    for chunk in pd.read_csv(chip_ann_file, sep='\t', usecols=cols_to_read,
                             chunksize=chunk_size, low_memory=False):
        # Filter for outlier SNPs
        matches = chunk[chunk['ID'].isin(outlier_snps)]
        if len(matches) > 0:
            annotations.append(matches)

    if annotations:
        result = pd.concat(annotations, ignore_index=True)

        # Parse ANN field to extract key information
        def parse_ann(ann_str):
            if pd.isna(ann_str) or ann_str == '.':
                return {
                    'effect': 'NA',
                    'impact': 'NA',
                    'gene': 'NA',
                    'feature_type': 'NA',
                    'transcript': 'NA',
                    'biotype': 'NA',
                    'hgvs_c': 'NA',
                    'hgvs_p': 'NA'
                }

            # Take first annotation (most severe)
            first_ann = ann_str.split(',')[0]
            fields = first_ann.split('|')

            return {
                'effect': fields[1] if len(fields) > 1 else 'NA',
                'impact': fields[2] if len(fields) > 2 else 'NA',
                'gene': fields[3] if len(fields) > 3 else 'NA',
                'feature_type': fields[5] if len(fields) > 5 else 'NA',
                'transcript': fields[6] if len(fields) > 6 else 'NA',
                'biotype': fields[7] if len(fields) > 7 else 'NA',
                'hgvs_c': fields[9] if len(fields) > 9 else 'NA',
                'hgvs_p': fields[10] if len(fields) > 10 else 'NA'
            }

        # Parse annotations
        parsed = result['ANN'].apply(parse_ann).apply(pd.Series)
        result = pd.concat([result[['CHROM', 'POS', 'ID', 'TYPE', 'LOF', 'NMD']], parsed], axis=1)

        # Save results
        result.to_csv(output_file, sep='\t', index=False)
        print(f"Saved {len(result)} annotated outliers to {output_file}")

        return result
    else:
        print("No outlier SNPs found in chip annotations")
        return None

def summarize_annotations(df, output_prefix):
    """Create summary statistics of annotations"""

    if df is None or len(df) == 0:
        return

    # Effect summary
    effect_counts = df['effect'].value_counts()
    effect_counts.to_csv(f"{output_prefix}_effect_summary.txt", sep='\t', header=['count'])

    # Impact summary
    impact_counts = df['impact'].value_counts()
    impact_counts.to_csv(f"{output_prefix}_impact_summary.txt", sep='\t', header=['count'])

    # Gene list
    genes = df[df['gene'] != 'NA']['gene'].unique()
    with open(f"{output_prefix}_genes.txt", 'w') as f:
        for gene in sorted(genes):
            f.write(f"{gene}\n")

    print("\nSummary:")
    print(f"  Total annotated SNPs: {len(df)}")
    print(f"  Unique genes: {len(genes)}")
    print("\n  Impact distribution:")
    for impact, count in impact_counts.items():
        print(f"    {impact}: {count}")
    print("\n  Top effects:")
    for effect, count in effect_counts.head(10).items():
        print(f"    {effect}: {count}")

def main():
    project_dir = Path("$HOME/projects/albopictus-autogeny")
    chip_ann = project_dir / "data/files/chip_ann.txt"
    output_dir = project_dir / "output/annotations"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Process pcadapt outliers (original workflow)
    pcadapt_outliers = project_dir / "output/selection_scans/pcadapt/outlier_snps.txt"
    if pcadapt_outliers.exists():
        print("\n=== Processing pcadapt outliers ===")
        outliers = load_outlier_snps(pcadapt_outliers)
        print(f"Loaded {len(outliers)} outlier SNPs")

        result = extract_annotations(chip_ann, outliers,
                                    output_dir / "pcadapt_outliers_annotated.tsv")
        summarize_annotations(result, str(output_dir / "pcadapt_outliers"))

    # Process OutFLANK outliers
    outflank_outliers = project_dir / "output/selection_scans/outflank/outlier_snps.txt"
    if outflank_outliers.exists():
        print("\n=== Processing OutFLANK outliers ===")
        outliers = load_outlier_snps(outflank_outliers)
        print(f"Loaded {len(outliers)} outlier SNPs")

        result = extract_annotations(chip_ann, outliers,
                                    output_dir / "outflank_outliers_annotated.tsv")
        summarize_annotations(result, str(output_dir / "outflank_outliers"))

    print("\n=== Annotation extraction complete ===")

if __name__ == "__main__":
    main()

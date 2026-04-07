# 03_annotation - Functional Annotation

> **Population label key:** AUT = AUTO, MAN = NON-AUTO, NEW = NON-AUTO-FIELD.

Annotate outlier SNPs with gene information.

## Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| `run_snpeff.sh` | SLURM batch for SnpEff annotation | `bash scripts/cli/03_annotation/run_snpeff.sh` |
| `extract_annotations.py` | Extract annotations from chip_ann.txt | `python extract_annotations.py` |

## Methods

### Primary: chip_ann.txt extraction
The SNP chip annotation file (`data/files/chip_ann.txt`) contains pre-computed SnpEff annotations for all SNPs. We extract annotations for outlier SNPs directly from this file.

### Alternative: Direct SnpEff (partial)
Direct SnpEff database build had GFF structural issues. Use chip_ann.txt instead.

## Input

- `data/files/chip_ann.txt` - SNP chip annotation (379MB)
- `data/files/genes.gff` - Gene annotation
- Outlier SNP lists from selection scans

## Output

- `output/annotations/` - Annotated outlier SNPs
  - `pcadapt_annotated.txt`
  - `outflank_annotated.txt`

## Annotation Results

### pcadapt outliers (41 SNPs)
| Impact | Count |
|--------|-------|
| MODERATE | 1 (missense) |
| LOW | 5 (synonymous) |
| MODIFIER | 12 (introns, UTRs) |
| Unannotated | 23 |

Key variant: AX-581302901 (chr3.17:12036013) - Leu156Phe missense

### OutFLANK outliers (63 SNPs from intergenic calibration)
- All MODIFIER impact (intergenic by design)

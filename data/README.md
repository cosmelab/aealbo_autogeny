# Data

Input data for the analysis notebooks. Large files are not tracked
in git -- obtain them from the sources below.

## Files in this directory

- `meta_data.txt` -- sample metadata (populations, IDs)

## Required data (not in git)

Obtain from Dropbox or Dryad before running notebooks.

### genotype_calls/

| File | Size | Description |
|------|------|-------------|
| `autogenous.vcf` | 42 MB | Full VCF (all populations) |
| `aut.vcf` | 23 MB | AUTO population VCF |
| `man.vcf` | 14 MB | NON-AUTO population VCF |
| `new.vcf` | 20 MB | NON-AUTO-FIELD population VCF |

### files/

| File | Size | Description |
|------|------|-------------|
| `chip_ann.txt` | 397 MB | SNP chip annotations |
| `genes.gff` | 108 MB | Gene annotations (GFF3) |
| `genes.txt` | 975 KB | Gene list |
| `MANvsAUTO_sig_mRNAs.csv` | 91 KB | DE gene expression results |
| `intergenic_SNPs.txt` | 179 KB | Intergenic SNP list |

### genome/

Reference genome files (from NCBI or Dryad).

### raw_data/

Raw genotyping data (from Dryad deposit `doi:10.5061/dryad.47v3c`).

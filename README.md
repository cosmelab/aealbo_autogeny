[![Container Build](https://github.com/cosmelab/aealbo_autogeny/actions/workflows/docker-build.yml/badge.svg)](https://github.com/cosmelab/aealbo_autogeny/actions/workflows/docker-build.yml) [![GitHub Pages](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://cosmelab.github.io/aealbo_autogeny/)

<!-- ![Banner](docs/images/banner.png) -->

# Autogeny in *Aedes albopictus*: Supplementary Analysis Notebooks

This repository contains the supplementary R Markdown notebooks for [Sturiale et al., BMC Biology, in revision]. We used a high-density SNP chip (131,048 SNPs) to analyze 60 *Aedes albopictus* mosquitoes across three populations: AUTO (autogenous selected line, n=28), NON-AUTO (Manassas VA lab colony, n=10), and NON-AUTO-FIELD (Montclair NJ field-collected, n=22). Analyses include quality control, genome-wide selection scans (OutFLANK, pcadapt), linkage disequilibrium network analysis (LDna), functional SNP annotation, and population genetic statistics.

## Supplementary Files

| File | Analysis | View |
|------|----------|------|
| File S1 | Quality control | [HTML](https://cosmelab.github.io/aealbo_autogeny/html/File_S1.Quality_control.html) |
| File S2 | Selection scan with OutFLANK | [HTML](https://cosmelab.github.io/aealbo_autogeny/html/File_S2_OutFlank.html) |
| File S3 | Selection scan with pcadapt | [HTML](https://cosmelab.github.io/aealbo_autogeny/html/File_S3.pcadapt.html) |
| File S4 | Linkage network analysis (LDna) | [HTML](https://cosmelab.github.io/aealbo_autogeny/html/File_S4.Linkage_network_analysis.html) |
| File S5 | Gene expression x selection scan intersection | [HTML](https://cosmelab.github.io/aealbo_autogeny/html/File_S5.Intersection_gene_expression_SNPs.html) |
| File S6 | Functional annotation of SNPs | [HTML](https://cosmelab.github.io/aealbo_autogeny/html/File_S6.SNPs_functional_annotation.html) |
| File S7 | Allele and genotype frequencies | [HTML](https://cosmelab.github.io/aealbo_autogeny/html/File_S7.Frequencies.html) |
| File S8 | FST estimations | [HTML](https://cosmelab.github.io/aealbo_autogeny/html/File_S8.Fst.html) |

## Quick Start

```bash
# 1. Clone
git clone https://github.com/cosmelab/aealbo_autogeny && cd aealbo_autogeny
# 2. Pull container (Singularity / HPC)
singularity pull docker://ghcr.io/cosmelab/aealbo_autogeny:latest
# OR Docker (local)
docker pull ghcr.io/cosmelab/aealbo_autogeny:latest
# 3. Render a notebook
singularity exec --bind $PWD:/workspace --pwd /workspace aealbo_autogeny_latest.sif \
  Rscript -e "rmarkdown::render('notebooks/File_S1.Quality_control.Rmd', output_dir='docs/html/')"
```

Note: See `scripts/render_notebooks.sh` for full pipeline and `docs/spark_instructions.md` for HPC.

## Data Availability

Raw genotype data (VCF files) and supporting files are available at Zenodo [DOI: TBD — will be updated upon acceptance]. Reference genome AalbF3 is available at NCBI (GCA_006496715.1). See `data/README.md` for complete data provenance.

## Citation

Sturiale SL, Heilig MC, Aardema ML, Cosme LV, Corley M, Marzec S, Hamilton M, Vizcarra D, Anderson L, Holzapfel CM, Bradshaw WE, Meuti ME, Caccone A, Armbruster PA. [Title]. BMC Biology [in revision]. GitHub: https://github.com/cosmelab/aealbo_autogeny

## License

MIT

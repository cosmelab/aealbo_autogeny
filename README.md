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

## Reproducing the Analysis

### Requirements

- Docker ≥ 20.x **or** Singularity ≥ 3.x (for HPC)
- ~10 GB disk space (container + data)
- `wget` (for data download)

### Step 1 — Clone the repo

```bash
git clone https://github.com/cosmelab/aealbo_autogeny && cd aealbo_autogeny
```

### Step 2 — Download data

```bash
bash scripts/00_download_data.sh
```

This downloads VCF files and annotation inputs from Zenodo [DOI: TBD]. The reference genome (AalbF3, ~1.9 GB) is fetched separately from NCBI — see the script's printed instructions.

### Step 3 — Pull the container

```bash
# Docker (local workstation)
docker pull ghcr.io/cosmelab/aealbo_autogeny:latest

# Singularity (HPC / Spark)
singularity pull docker://ghcr.io/cosmelab/aealbo_autogeny:latest
```

### Step 4 — Render notebooks

```bash
# Render all notebooks in dependency order (auto-detects Docker vs Singularity)
bash scripts/render_notebooks.sh

# Render a single notebook
docker run --rm -v $PWD:/workspace --workdir /workspace \
    ghcr.io/cosmelab/aealbo_autogeny:latest \
    bash -c 'eval "$(pixi shell-hook)" && Rscript -e \
    "rmarkdown::render(\"notebooks/File_S1.Quality_control.Rmd\", output_dir=\"docs/html/\")"'
```

See `docs/spark_instructions.md` for the full HPC workflow and `docs/spark_agent_task.json` for agent-driven rendering.

## Data Availability

Input data (VCF files, SNP chip annotations, gene annotations) are archived at Zenodo [DOI: TBD — updated upon acceptance]. Reference genome AalbF3 is available at NCBI (GCA_006496715.1). See `data/README.md` for complete provenance and `docs/zenodo_manifest.md` for the full archive contents.

## Citation

Sturiale SL, Heilig MC, Aardema ML, Cosme LV, Corley M, Marzec S, Hamilton M, Vizcarra D, Anderson L, Holzapfel CM, Bradshaw WE, Meuti ME, Caccone A, Armbruster PA. [Title]. BMC Biology [in revision]. GitHub: https://github.com/cosmelab/aealbo_autogeny

## License

MIT

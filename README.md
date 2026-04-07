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

Choose the method that fits your setup. All three produce identical results.

---

### Option A — pixi (local, no container)

[pixi](https://pixi.sh) manages the full R + Python + genomics tool environment from `pixi.toml`. No Docker or Singularity needed.

```bash
# 1. Clone
git clone https://github.com/cosmelab/aealbo_autogeny && cd aealbo_autogeny

# 2. Install pixi (if not already installed)
curl -fsSL https://pixi.sh/install.sh | bash

# 3. Install all dependencies (R, Python, plink2, etc.)
pixi install
pixi run install-r-extras   # installs OutFLANK, LDna, and other CRAN/GitHub packages

# 4. Download data
bash scripts/00_download_data.sh

# 5. Render a notebook
pixi run Rscript -e "rmarkdown::render('notebooks/File_S1.Quality_control.Rmd', output_dir='docs/html/')"

# Render all notebooks in dependency order
bash scripts/render_notebooks.sh
```

---

### Option B — Docker (local workstation)

```bash
# 1. Clone and download data
git clone https://github.com/cosmelab/aealbo_autogeny && cd aealbo_autogeny
bash scripts/00_download_data.sh

# 2. Pull container
docker pull ghcr.io/cosmelab/aealbo_autogeny:latest

# 3. Render a notebook
docker run --rm -v $PWD:/workspace --workdir /workspace \
    ghcr.io/cosmelab/aealbo_autogeny:latest \
    bash -c 'eval "$(pixi shell-hook)" && Rscript -e \
    "rmarkdown::render(\"notebooks/File_S1.Quality_control.Rmd\", output_dir=\"docs/html/\")"'

# Render all
bash scripts/render_notebooks.sh   # auto-detects Docker
```

---

### Option C — Singularity (HPC / Spark)

```bash
# 1. Clone and download data
git clone https://github.com/cosmelab/aealbo_autogeny && cd aealbo_autogeny
bash scripts/00_download_data.sh

# 2. Pull container
singularity pull docker://ghcr.io/cosmelab/aealbo_autogeny:latest

# 3. Render a notebook
singularity exec --bind $PWD:/workspace --pwd /workspace \
    aealbo_autogeny_latest.sif \
    bash -c 'eval "$(pixi shell-hook)" && Rscript -e \
    "rmarkdown::render(\"notebooks/File_S1.Quality_control.Rmd\", output_dir=\"docs/html/\")"'

# Render all
bash scripts/render_notebooks.sh   # auto-detects Singularity
```

See `docs/spark_instructions.md` for the full HPC workflow and `docs/spark_agent_task.json` for agent-driven rendering.

## CLI Analysis Pipeline

For programmatic or HPC use, run the full pipeline from the command line without opening notebooks. Each step maps directly to its supplementary notebook.

```bash
# Check what's pending (dry run — nothing executes)
bash scripts/run_pipeline.sh --dry-run

# Run the full pipeline (skips already-completed steps)
bash scripts/run_pipeline.sh
```

Requires the container and data (see Reproducing the Analysis above). The pipeline auto-detects Docker or Podman via `CONTAINER_RUNTIME`.

| Step | Script | Output | Notebook |
|------|--------|--------|----------|
| 01 | `scripts/01_qc/run_qc.sh` | `output/quality_control/file7.*` | File S1 |
| 02 | `scripts/02_selection_scans/run_selection_scans.sh` | `output/selection_scans/` | Files S2, S3 |
| 03 | `scripts/05_ldna/run_ldna.sh` | `output/ldna/` | Files S4a–S4e |
| 04 | `scripts/03_annotation/run_snpeff.sh` | `output/snpeff/` | File S6 |
| 05 | `scripts/07_gene_expression/run_gene_expression.sh` | `output/gene_expression/` | File S5 |
| 06 | `scripts/04_diversity/run_diversity.sh` | `output/diversity/`, `output/fst/` | Files S7, S8 |

> **Note:** Step 03 (LDna) requires ~45 min and 32 GB RAM. See `scripts/README.md` for per-step details.

## Data Availability

Input data (VCF files, SNP chip annotations, gene annotations) are archived at Zenodo [DOI: TBD — updated upon acceptance]. Reference genome AalbF3 is available at NCBI (GCA_006496715.1). See `data/README.md` for complete provenance and `docs/zenodo_manifest.md` for the full archive contents.

## Citation

Sturiale SL, Heilig MC, Aardema ML, Cosme LV, Corley M, Marzec S, Hamilton M, Vizcarra D, Anderson L, Holzapfel CM, Bradshaw WE, Meuti ME, Caccone A, Armbruster PA. [Title]. BMC Biology [in revision]. GitHub: https://github.com/cosmelab/aealbo_autogeny

## License

MIT

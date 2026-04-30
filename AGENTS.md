# AGENTS.md — aealbo_autogeny

Supplementary analysis repository for Sturiale et al. (*BMC Biology*, in revision): selection scans and LD network analysis in *Aedes albopictus* autogeny. GitHub Pages site: <https://cosmelab.github.io/aealbo_autogeny/>.

## Development workflow

Active development happens in the companion dev repo (`albopictus-autogeny`) on Spark. This repo holds the published notebooks and rendered HTML. Changes made here should be kept in sync with the dev repo.

## Repository structure

```text
notebooks/     R Markdown source (File_S1–S8; S4 split into S4a–S4e)
docs/html/     Rendered HTML for GitHub Pages
scripts/       CLI reproducibility pipeline (see scripts/README.md)
data/          Input data (downloaded via scripts/00_download_data.sh)
```

## Notebook index

| File | Analysis |
|------|----------|
| File_S1 | Quality control |
| File_S2 | Selection scan — OutFLANK |
| File_S3 | Selection scan — pcadapt |
| File_S4a | LDna data preparation |
| File_S4b | LDna chromosome 1 |
| File_S4c | LDna chromosome 2 |
| File_S4d | LDna chromosome 3 |
| File_S4e | LDna combined / cluster mapping |
| File_S5 | Gene expression × selection scan intersection |
| File_S6 | Functional annotation of SNPs |
| File_S7 | Allele and genotype frequencies |
| File_S8 | FST estimations |

File_S4 (`docs/html/File_S4.Linkage_network_analysis.html`) is a landing page linking to S4a–S4e; there is no corresponding monolithic Rmd.

## Population labels

| Code in data | Manuscript label | Description |
|---|---|---|
| AUT | AUTO | Autogenous selected line (n=28) |
| MAN | NON-AUTO | Non-autogenous lab colony, Manassas VA (n=10) |
| NEW | NON-AUTO-FIELD | Non-autogenous field sample, Montclair NJ (n=22) |

MAN (n=10) is excluded from LD network analysis due to insufficient sample size for reliable LD estimation.

## Rendering notebooks

### Option A — pixi (local, no container)

```bash
pixi install
pixi run install-r-extras
bash scripts/00_download_data.sh
pixi run Rscript -e "rmarkdown::render('notebooks/File_S1.Quality_control.Rmd', output_dir='docs/html/')"
bash scripts/render_notebooks.sh   # renders all in dependency order
```

### Option B — container (podman or docker)

```bash
# Podman
CONTAINER_RUNTIME=podman bash scripts/render_notebooks.sh

# Docker (default if CONTAINER_RUNTIME not set)
bash scripts/render_notebooks.sh
```

Container image: `ghcr.io/cosmelab/albopictus-autogeny:latest`

LDna rendering (~45 min, requires 32 GB RAM):

```bash
CONTAINER_RUNTIME=podman bash scripts/05_ldna/run_ldna.sh
```

## CLI pipeline (full reproducibility)

The `scripts/` directory contains a command-line pipeline that reproduces all analyses without the notebooks. Run order:

```text
01_qc → 02_selection_scans → 05_ldna → 07_gene_expression → 03_annotation → 04_diversity
```

See `scripts/README.md` for full details and per-step output descriptions.

## Key results (for sanity checks)

| Check | Expected |
|---|---|
| Samples after QC | 60 |
| SNPs after QC | ~33,836 |
| LDna cluster AUTO chr2 | `3215_0.78` (324 SNPs) |
| Pairwise Fst AUT vs MAN | 0.13 |
| Pairwise Fst MAN vs NEW | 0.04 |

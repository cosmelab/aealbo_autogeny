# CLI Analysis Pipeline

Reproducible command-line pipeline for autogeny GWAS in *Aedes albopictus*.
Runs the same analyses as the R Markdown supplementary notebooks using Docker or Podman.

> **Population label key:** The code labels used in this pipeline correspond to the following manuscript labels: AUT = AUTO, MAN = NON-AUTO, NEW = NON-AUTO-FIELD.

## Container

```bash
ghcr.io/cosmelab/albopictus-autogeny:latest
```

Works with Docker (default) or Podman:

```bash
# Docker (default)
bash scripts/cli/01_qc/run_qc.sh

# Podman
CONTAINER_RUNTIME=podman bash scripts/cli/01_qc/run_qc.sh
```

## Pipeline — Run Order

```
01_qc → 02_selection_scans → 05_ldna → 07_gene_expression → 03_annotation → 04_diversity
```

All run scripts accept an optional `repo_dir` argument (defaults to `$HOME/projects/albopictus-autogeny`).

### 1. Quality Control → File S1

```bash
bash scripts/cli/01_qc/run_qc.sh
```

Output: `output/quality_control/file7` — 60 samples, ~33,836 common SNPs

### 2. Selection Scans → Files S2, S3

```bash
bash scripts/cli/02_selection_scans/run_selection_scans.sh
bash scripts/cli/02_selection_scans/run_pairwise_pcadapt.sh
```

Output: `output/selection_scans/`, `output/pcadapt/`

> **Note:** Pairwise pcadapt between NON-AUTO and NON-AUTO-FIELD was not used in the final manuscript.

### 3. LD Network Analysis → Files S4a–S4e

```bash
bash scripts/cli/05_ldna/run_ldna.sh           # ~45 min, 32GB RAM required
bash scripts/cli/05_ldna/run_ldna_visualization.sh
```

Output: `output/ldna/`

### 4. Gene Expression Intersection → File S5

```bash
bash scripts/cli/07_gene_expression/run_gene_expression.sh
```

Output: `output/gene_expression/`

> **Note:** Full SNP-to-gene coordinate mapping (scaffold → chromosome) is implemented in `notebooks/File_S5.Intersection_gene_expression_SNPs.Rmd`. This CLI script provides a summary only.

### 5. Functional Annotation → File S6

```bash
bash scripts/cli/03_annotation/run_snpeff.sh
```

Output: `output/snpeff/`

### 6. Diversity and Fst → Files S7, S8

```bash
bash scripts/cli/04_diversity/run_diversity.sh
```

Output: `output/diversity/`, `output/fst/`, `output/frequencies/`

> **Note:** Sex-stratified Fst (S8 Sections 2 and 4) was not used in the final manuscript and is not implemented in this pipeline.

### 7. Enrichment Analysis

```bash
bash scripts/cli/09_enrichment/run_enrichment.sh
```

> **Note:** Functional enrichment analysis was not included in the final manuscript.

### 8. Power Analysis

```bash
bash scripts/cli/10_power_analysis/run_power_analysis.sh
```

> **Note:** Power analysis is provided as supplementary context and was not included in the final manuscript.

## Directory Structure

| Directory | Notebook(s) | Description |
|-----------|------------|-------------|
| `01_qc/` | File S1 | Quality control (7 steps) |
| `02_selection_scans/` | Files S2, S3 | OutFLANK + pcadapt |
| `03_annotation/` | File S6 | SnpEff functional annotation |
| `04_diversity/` | Files S7, S8 | He/Ho, pairwise Fst, frequency plots |
| `05_ldna/` | Files S4a–S4e | LDna network analysis |
| `07_gene_expression/` | File S5 | Gene expression intersection (summary) |
| `09_enrichment/` | — | Enrichment analysis (not in manuscript) |
| `10_power_analysis/` | — | Power analysis (not in manuscript) |
| `12_ldna_figures/` | Files S4e, S7, S8 | Publication figure scripts |

## Population Labels

| Code in data | Manuscript label | Description |
|-------------|-----------------|-------------|
| AUT | AUTO | Autogenous (selected line) |
| MAN | NON-AUTO | Non-autogenous (Manassas, VA colony) |
| NEW | NON-AUTO-FIELD | Non-autogenous (Montclair, NJ field) |

## Key Results

| Analysis | Result |
|----------|--------|
| QC final | 60 samples, ~33,836 common SNPs |
| LDna clusters AUTO chr2 | Cluster 14 (`3215_0.78`, 324 SNPs) |
| LDna clusters AUTO chr1 | Cluster 6 (`1971_0.64`, 90 SNPs) |
| Pairwise Fst AUT vs MAN | 0.13 |
| Pairwise Fst AUT vs NEW | 0.13 |
| Pairwise Fst MAN vs NEW | 0.04 |

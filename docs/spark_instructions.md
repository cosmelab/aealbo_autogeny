# Rendering Notebooks on Spark HPC

Step-by-step instructions for rendering all supplementary notebooks on the
Spark high-performance computing cluster using Singularity.

---

## Prerequisites

### 1. Clone the repository

```bash
git clone https://github.com/cosmelab/aealbo_autogeny
cd aealbo_autogeny
```

### 2. Pull the container image

```bash
singularity pull docker://ghcr.io/cosmelab/aealbo_autogeny:latest
```

This produces `aealbo_autogeny_latest.sif` in the current directory.

### 3. Sync input data

Option A — sync from Dropbox:

```bash
rsync -avz user@dropbox-host:/path/to/data/ data/
```

Option B — download from Zenodo (use when DOI is published):

```bash
# Replace DOI_HERE with the actual Zenodo DOI when available
wget https://zenodo.org/record/DOI_HERE/files/aealbo_autogeny_data.tar.gz
tar -xzvf aealbo_autogeny_data.tar.gz
```

---

## Render Command Template

All notebooks are rendered with Singularity using the pattern below.
Replace `NOTEBOOK.Rmd` with the target notebook filename.

```bash
singularity exec --bind $PWD:/workspace --pwd /workspace \
  aealbo_autogeny_latest.sif \
  Rscript -e "rmarkdown::render('notebooks/NOTEBOOK.Rmd', output_dir='docs/html/')"
```

The `--bind $PWD:/workspace --pwd /workspace` flags are required so that `here()`
resolves paths correctly inside the container.

---

## CRITICAL BLOCKER — File_S3.pcadapt.Rmd

`File_S3.pcadapt.Rmd` contains Python chunks that call `matplotlib_venn`,
which is **not installed in the container**. Rendering will fail without
the following fix.

Before rendering S3, open the notebook and set `eval=FALSE` on the Python
chunks at **lines 335, 820, and 1300**:

```
# Change each affected chunk header from:
{python ...}
# to:
{python eval=FALSE, ...}
```

Do this once before submitting the render job. The rest of the notebook
renders normally.

---

## Render Order

Notebooks have data dependencies. Follow this order. Steps at the same
numbered level may be run in parallel.

```
Step 1:  File_S1.Quality_control.Rmd
         # Generates output/quality_control/file7.*
         # Required by steps 2, 3, and 4

Step 2:  File_S2_OutFlank.Rmd
         # Needs: file7.* from Step 1

Step 3:  File_S3.pcadapt.Rmd
         # Needs: file7.* from Step 1
         # BLOCKER: set eval=FALSE on Python chunks first (lines 335, 820, 1300)

Step 4:  File_S4a.LDna_data_prep.Rmd
         # Needs: file7.* from Step 1

Step 5:  File_S4b.LDna_chr1.Rmd       (parallel)
         File_S4c.LDna_chr2.Rmd       (parallel)
         File_S4d.LDna_chr3.Rmd       (parallel)
         # All need: output from Step 4 (S4a)

Step 6:  File_S4e.LDna_combined.Rmd
         # Needs: output from Step 5 (S4b, S4c, S4d)

Step 7:  File_S4.Linkage_network_analysis.Rmd
         # Needs: output from Step 6 (S4e)
         # Note: this notebook takes approximately 2 hours

Step 8a: Run SnpEff annotation first:
         bash scripts/03_annotation/run_snpeff.sh

Step 8:  File_S5.Intersection_gene_expression_SNPs.Rmd
         # Needs: SnpEff output from Step 8a

Step 9:  File_S6.SNPs_functional_annotation.Rmd
         # Needs: SnpEff output from Step 8a

Step 10: File_S7.Frequencies.Rmd

Step 11: File_S8.Fst.Rmd
```

Alternatively, run all notebooks in the correct order automatically:

```bash
bash scripts/render_notebooks.sh
```

---

## After Rendering

Commit the rendered HTML files and push to GitHub:

```bash
git add docs/html/
git commit -m "Add rendered HTML notebooks"
git push origin main
```

---

## Enable GitHub Pages

After pushing, enable GitHub Pages to serve the HTML notebooks publicly:

1. Go to: https://github.com/cosmelab/aealbo_autogeny/settings/pages
2. Under **Source**, select: branch `main`, folder `/docs` → click **Save**
3. Wait approximately 2 minutes for deployment to complete
4. Verify: https://cosmelab.github.io/aealbo_autogeny/

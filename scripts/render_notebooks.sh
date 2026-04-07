#!/usr/bin/env bash
set -euo pipefail

CONTAINER="ghcr.io/cosmelab/aealbo_autogeny:latest"
REPO_DIR=$(realpath "$(dirname "$0")/..")
SIF="${REPO_DIR}/aealbo_autogeny_latest.sif"
FORCE=0

# Parse flags
for arg in "$@"; do
  case "$arg" in
    --force) FORCE=1 ;;
  esac
done

# Auto-detect runtime
if command -v singularity &>/dev/null; then
  RUNTIME="singularity"
else
  RUNTIME="docker"
fi

echo "[$(date +%H:%M:%S)] Runtime detected: $RUNTIME"

render_notebook() {
  local NB="$1"
  local HTML_OUT="${REPO_DIR}/docs/html/${NB%.Rmd}.html"

  if [[ "$FORCE" -eq 0 && -f "$HTML_OUT" ]]; then
    echo "[$(date +%H:%M:%S)] Skipping $NB (HTML already exists; use --force to re-render)"
    return 0
  fi

  echo "[$(date +%H:%M:%S)] Rendering $NB..."

  if [[ "$RUNTIME" == "singularity" ]]; then
    singularity exec \
      --bind "${REPO_DIR}":/workspace \
      --pwd /workspace \
      "$SIF" \
      Rscript -e "rmarkdown::render('notebooks/${NB}', output_dir='docs/html/')"
  else
    docker run --rm \
      -v "${REPO_DIR}":/workspace \
      --workdir /workspace \
      --entrypoint bash \
      "$CONTAINER" \
      --norc -c "eval \"\$(pixi shell-hook)\" && Rscript -e \"rmarkdown::render('notebooks/${NB}', output_dir='docs/html/')\""
  fi

  echo "[$(date +%H:%M:%S)] Done: $NB"
}

warn_s3_blocker() {
  echo ""
  echo "WARNING: File_S3.pcadapt.Rmd contains Python chunks that use matplotlib_venn,"
  echo "         which is NOT installed in the container."
  echo "         Before rendering, ensure eval=FALSE is set on the python chunk headers"
  echo "         at lines 335, 820, and 1300 in notebooks/File_S3.pcadapt.Rmd."
  echo ""
  echo "Continuing in 5 seconds (Ctrl+C to abort)..."
  sleep 5
}

main() {
  mkdir -p "${REPO_DIR}/docs/html"

  render_notebook "File_S1.Quality_control.Rmd"
  render_notebook "File_S2_OutFlank.Rmd"

  warn_s3_blocker
  render_notebook "File_S3.pcadapt.Rmd"

  render_notebook "File_S4a.LDna_data_prep.Rmd"
  render_notebook "File_S4b.LDna_chr1.Rmd"
  render_notebook "File_S4c.LDna_chr2.Rmd"
  render_notebook "File_S4d.LDna_chr3.Rmd"
  render_notebook "File_S4e.LDna_combined.Rmd"
  render_notebook "File_S4.Linkage_network_analysis.Rmd"
  render_notebook "File_S5.Intersection_gene_expression_SNPs.Rmd"
  render_notebook "File_S6.SNPs_functional_annotation.Rmd"
  render_notebook "File_S7.Frequencies.Rmd"
  render_notebook "File_S8.Fst.Rmd"

  echo ""
  echo "[$(date +%H:%M:%S)] All notebooks rendered."
  echo "Run: git add docs/html/ && git commit -m 'Render notebooks' && git push"
}

main

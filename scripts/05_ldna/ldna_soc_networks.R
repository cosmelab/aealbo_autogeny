#!/usr/bin/env Rscript
# ============================================================
# LDna SOC Network Visualization
# ============================================================
# Uses the native plotLDnetwork(option=2) from LDna v2.15 to
# visualize each SOC (Single Outlier Cluster) as an actual LD
# network — SNPs as nodes, LD edges, colored by chromosome
# position (red → green via pos= argument).
#
# Outputs one PDF per population × chromosome with multi-panel
# network plots for the top SOCs by cluster size.
#
# Usage (inside container):
#   Rscript scripts/cli/05_ldna/ldna_soc_networks.R [project_dir] [min_loci] [min_edges]
#
# Defaults: min_loci=20, min_edges=100
# ============================================================

suppressPackageStartupMessages({
  .libPaths(c("/workspace/lib/ldna_v215", .libPaths()))
  library(LDna)
  library(data.table)
})

# ---- Arguments -----------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[1] else "/project"
min_loci    <- if (length(args) >= 2) as.integer(args[2]) else 20
min_edges   <- if (length(args) >= 3) as.integer(args[3]) else 100

setwd(project_dir)

outdir <- file.path("output", "ldna", "soc_networks")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("LDna SOC Network Visualization\n")
cat("  Project dir :", project_dir, "\n")
cat("  Min loci    :", min_loci, "\n")
cat("  Min edges   :", min_edges, "(extractBranches parameter)\n")
cat("  Output dir  :", outdir, "\n")
cat("============================================================\n\n")

# ---- Helpers -------------------------------------------------

import_bim <- function(path) {
  fread(path, header = FALSE,
        col.names = c("CHR", "SNP", "cM", "Position", "A1", "A2"))
}

# Build a named numeric vector of positions for plotLDnetwork pos=
make_pos_vector <- function(snp_names, bim) {
  pos <- bim$Position[match(snp_names, bim$SNP)]
  names(pos) <- snp_names
  pos
}

# ---- Main loop -----------------------------------------------

populations  <- c("AUT", "NEW")
chromosomes  <- c("chr1", "chr2", "chr3")
chr_numbers  <- c(1, 2, 3)

for (pop in populations) {
  for (ci in seq_along(chromosomes)) {
    chr    <- chromosomes[ci]
    chrnum <- chr_numbers[ci]

    rds_ldna <- file.path("output", "ldna", "pop", chr,
                          paste0(pop, ".rds"))
    rds_ld2  <- file.path("output", "ldna", "pop", chr,
                          paste0(pop, "_ld2.rds"))
    bim_file <- file.path("output", "ldna", "files", "file1.bim")

    if (!file.exists(rds_ldna) || !file.exists(rds_ld2)) {
      cat("SKIP:", pop, chr, "— rds files not found\n")
      next
    }

    cat("------------------------------------------------------------\n")
    cat("Processing:", pop, chr, "\n")

    ldna <- readRDS(rds_ldna)
    ld2  <- readRDS(rds_ld2)
    bim  <- import_bim(bim_file)

    # plotLDnetwork cannot handle NAs in the LD matrix
    na_count <- sum(is.na(ld2))
    if (na_count > 0) {
      cat("  Replacing", na_count, "NAs in LD matrix with 0\n")
      ld2[is.na(ld2)] <- 0
    }

    # Extract SOCs
    clusters <- extractBranches(ldna, min.edges = min_edges,
                                plot.tree = FALSE, cores = 1)
    summary  <- summaryLDna(ldna, clusters, ld2)

    cat("  Total SOCs     :", length(clusters), "\n")

    # Filter to SOCs with >= min_loci SNPs
    keep <- summary$nLoci >= min_loci
    clusters_filt <- clusters[keep]
    summary_filt  <- summary[keep, ]

    cat("  SOCs >= ", min_loci, " loci: ", sum(keep), "\n", sep = "")

    if (length(clusters_filt) == 0) {
      cat("  No clusters to plot — skipping\n")
      next
    }

    # Sort by nLoci descending
    ord           <- order(summary_filt$nLoci, decreasing = TRUE)
    clusters_filt <- clusters_filt[ord]
    summary_filt  <- summary_filt[ord, ]

    # Cap at 16 panels per PDF page
    n_plot <- min(length(clusters_filt), 16)

    # ---- PDF output ------------------------------------------
    pdf_file <- file.path(outdir,
                          paste0(pop, "_", chr, "_soc_networks.pdf"))

    # Determine layout: up to 4×4
    n_cols <- min(4, ceiling(sqrt(n_plot)))
    n_rows <- ceiling(n_plot / n_cols)

    pdf(pdf_file, width = n_cols * 4, height = n_rows * 4 + 1)
    par(mfrow = c(n_rows, n_cols),
        mar = c(1, 1, 3.5, 1),
        oma = c(0, 0, 2, 0))

    for (i in seq_len(n_plot)) {
      cl_name  <- names(clusters_filt)[i]
      cl_snps  <- clusters_filt[[i]]
      cl_nloci <- summary_filt$nLoci[i]
      cl_ne    <- summary_filt$nE[i]
      cl_r2    <- round(summary_filt$Median.LD[i], 2)

      # Position vector for color gradient (red=low pos, green=high pos)
      pos_vec <- make_pos_vector(cl_snps, bim)

      cat("  Plotting SOC", i, "/", n_plot, ":", cl_name,
          "(", cl_nloci, "loci,", cl_ne, "edges )\n")

      # plotLDnetwork adds its own title: "cluster_name @threshold"
      # We add a subtitle with loci/edges/r2 via mtext (line 3 = below main)
      plotLDnetwork(
        ldna         = ldna,
        LDmat        = ld2,
        option       = 2,
        clusters     = clusters_filt[i],
        summary      = summary_filt[i, ],
        full.network = FALSE,
        pos          = pos_vec
      )

      mtext(
        text = paste0(cl_nloci, " loci · ", cl_ne, " edges · r²=", cl_r2),
        side = 3, line = 0.2, cex = 0.65, col = "grey40"
      )
    }

    # Outer title: population + chromosome
    mtext(
      text = paste0(pop, "  ·  ", chr,
                    "  (SOCs ≥ ", min_loci, " loci, top ", n_plot, " shown)"),
      outer = TRUE, side = 3, line = 0.5, cex = 1, font = 2
    )

    dev.off()
    cat("  Saved →", pdf_file, "\n\n")

    # ---- Also save a summary table ---------------------------
    tsv_file <- file.path(outdir,
                          paste0(pop, "_", chr, "_soc_summary.tsv"))
    write.table(summary_filt, tsv_file,
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

cat("============================================================\n")
cat("Done. Output in:", outdir, "\n")
cat("============================================================\n")

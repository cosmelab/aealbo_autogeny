#!/usr/bin/env Rscript
# Fix LDna Network Plots - Separate plots for each population
# Uses LDna package plotLDnetwork and extractClusters functions

library(here)

# Add LDna library path
.libPaths(c(here("output", "ldna", "lib_v215"), .libPaths()))
library(LDna)

# Create output directory
dir.create(here("output", "ldna", "network_figures"), showWarnings = FALSE, recursive = TRUE)

# Parameters (matching original analysis)
min_edges <- 20
lambda_lim <- 1

# ===========================================
# Function to generate network plot for a population
# ===========================================
generate_ldna_plots <- function(pop_name, chr, output_dir) {
  cat("\n========================================\n")
  cat("Processing:", pop_name, "Chromosome", chr, "\n")
  cat("========================================\n")

  # Load data
  ldna_file <- here("output", "ldna", "pop", paste0("chr", chr), paste0(pop_name, "_ldna.rds"))
  ld2_file <- here("output", "ldna", "pop", paste0("chr", chr), paste0(pop_name, "_ld2.rds"))
  clusters_file <- here("output", "ldna", "pop", paste0("chr", chr), paste0(pop_name, "_clusters.rds"))

  if (!file.exists(ldna_file) || !file.exists(ld2_file)) {
    cat("Files not found for", pop_name, "chr", chr, "\n")
    return(NULL)
  }

  cat("Loading ldna object...\n")
  ldna <- readRDS(ldna_file)

  cat("Loading LD matrix...\n")
  ld2 <- readRDS(ld2_file)

  cat("LD matrix dimensions:", dim(ld2), "\n")

  # Generate tree plot using extractBranches
  pdf_file <- file.path(output_dir, paste0(pop_name, "_chr", chr, "_ldna_tree.pdf"))
  cat("Creating tree plot:", pdf_file, "\n")

  pdf(file = pdf_file, width = 14, height = 12)

  # Extract branches with tree plot
  branches <- extractBranches(
    ldna,
    min.edges = min_edges,
    plot.tree = TRUE
  )

  # Add title
  title(main = paste(pop_name, "- Chromosome", chr, "- LDna Tree"),
        sub = paste("min.edges =", min_edges))

  dev.off()
  cat("Saved:", pdf_file, "\n")

  # Also save as PNG
  png_file <- file.path(output_dir, paste0(pop_name, "_chr", chr, "_ldna_tree.png"))
  cat("Creating PNG:", png_file, "\n")

  png(file = png_file, width = 1400, height = 1200, res = 100)

  branches <- extractBranches(
    ldna,
    min.edges = min_edges,
    plot.tree = TRUE
  )

  title(main = paste(pop_name, "- Chromosome", chr, "- LDna Tree"),
        sub = paste("min.edges =", min_edges))

  dev.off()
  cat("Saved:", png_file, "\n")

  # Now create network plot using plotLDnetwork
  # This uses the extracted branches to show the network
  if (!is.null(branches) && length(branches) > 0) {
    net_pdf <- file.path(output_dir, paste0(pop_name, "_chr", chr, "_ldna_network.pdf"))
    cat("Creating network plot:", net_pdf, "\n")

    pdf(file = net_pdf, width = 14, height = 12)

    # Get summary for plotLDnetwork
    summ <- summaryLDna(ldna, branches, ld2)

    plotLDnetwork(
      ldna,
      LDmat = ld2,
      option = 2,  # option 2 shows clusters
      threshold = 0.8,
      clusters = branches,
      summary = summ,
      full.network = FALSE
    )

    title(main = paste(pop_name, "- Chromosome", chr, "- LDna Network"))

    dev.off()
    cat("Saved:", net_pdf, "\n")

    # PNG version
    net_png <- file.path(output_dir, paste0(pop_name, "_chr", chr, "_ldna_network.png"))
    png(file = net_png, width = 1400, height = 1200, res = 100)

    plotLDnetwork(
      ldna,
      LDmat = ld2,
      option = 2,
      threshold = 0.8,
      clusters = branches,
      summary = summ,
      full.network = FALSE
    )

    title(main = paste(pop_name, "- Chromosome", chr, "- LDna Network"))

    dev.off()
    cat("Saved:", net_png, "\n")
  }

  # Load and report clusters if available
  if (file.exists(clusters_file)) {
    clusters <- readRDS(clusters_file)
    cat("Number of clusters:", length(clusters), "\n")
  }

  TRUE
}

# ===========================================
# Generate plots for each population and chromosome
# ===========================================
output_dir <- here("output", "ldna", "network_figures")

populations <- c("AUTO", "NON-AUTO-FIELD")
chromosomes <- c(1, 2, 3)

for (pop in populations) {
  for (chr in chromosomes) {
    tryCatch({
      generate_ldna_plots(pop, chr, output_dir)
    }, error = function(e) {
      cat("Error processing", pop, "chr", chr, ":", e$message, "\n")
    })
  }
}

cat("\n=== Done! ===\n")
cat("Network plots saved to:", output_dir, "\n")

#!/usr/bin/env Rscript
# Compare LDna v.65 (extractClusters) vs v2.15 (extractBranches)
#
# Key difference:
# - v.65 extractClusters() has rm.COCs=TRUE to remove nested clusters (SOCs only)
# - v2.15 extractBranches() treats ALL branches as clusters (includes nested)
#
# Author: Luciano Cosme
# Date: 2025-12-15

suppressPackageStartupMessages({
  library(remotes)
})

# Create separate library paths for each version
lib_v65 <- "output/ldna/lib_v65"
lib_v215 <- "output/ldna/lib_v215"

dir.create(lib_v65, recursive = TRUE, showWarnings = FALSE)
dir.create(lib_v215, recursive = TRUE, showWarnings = FALSE)

cat("=== Installing LDna versions ===\n\n")

# Install v.63 (has extractClusters with rm.COCs) - with dependencies
if (!file.exists(file.path(lib_v65, "LDna"))) {
  cat("Installing LDna v.63 (has extractClusters)...\n")
  remotes::install_github("petrikemppainen/LDna@v.63",
                          lib = lib_v65,
                          upgrade = "never",
                          dependencies = TRUE,
                          build_vignettes = FALSE,
                          quiet = FALSE)
}

# Install v2.15 in separate path
if (!file.exists(file.path(lib_v215, "LDna"))) {
  cat("Installing LDna v2.15...\n")
  remotes::install_github("petrikemppainen/LDna@v.2.15",
                          lib = lib_v215,
                          upgrade = "never",
                          dependencies = TRUE,
                          build_vignettes = FALSE,
                          quiet = FALSE)
}

cat("Using both LDna versions from separate lib paths...\n")

# Load test data
cat("\n=== Loading test data ===\n")
ld2_file <- "output/ldna/pop/chr1/AUTO_ld2.rds"
ldna_file <- "output/ldna/pop/chr1/AUTO_ldna.rds"

if (!file.exists(ld2_file)) {
  stop("LD matrix file not found: ", ld2_file)
}

ld2 <- readRDS(ld2_file)
ldna_raw <- readRDS(ldna_file)

cat("LD matrix size:", nrow(ld2), "x", ncol(ld2), "\n")

# Test parameters
min_edges_test <- 100

cat("\n=== Testing v.65 (extractClusters) ===\n")

# Load v.65 in isolated environment
library(LDna, lib.loc = lib_v65)

# Check if extractClusters exists
if (exists("extractClusters")) {
  cat("extractClusters found!\n")

  # Run with original parameters
  tryCatch({
    clusters_v65 <- extractClusters(
      ldna_raw,
      LDmat = ld2,
      min.edges = min_edges_test,
      lambda.lim = 1,  # Original parameter
      rm.COCs = TRUE,   # Remove nested clusters
      extract = TRUE,
      plot.tree = FALSE,
      plot.graph = FALSE
    )

    n_clusters_v65 <- length(clusters_v65$clusters)
    cat("v.65 clusters (min.edges=", min_edges_test, ", lambda.lim=1, rm.COCs=TRUE):", n_clusters_v65, "\n")

    # Get cluster sizes
    cluster_sizes_v65 <- sapply(clusters_v65$clusters, length)
    cat("  Mean SNPs per cluster:", round(mean(cluster_sizes_v65), 1), "\n")
    cat("  Total SNPs:", sum(cluster_sizes_v65), "\n")

  }, error = function(e) {
    cat("Error with v.65:", e$message, "\n")
  })

} else {
  cat("extractClusters NOT found in v.65!\n")
}

# Detach v.65
detach("package:LDna", unload = TRUE)

cat("\n=== Testing v2.15 (extractBranches) ===\n")

# Load v2.15 from system
library(LDna)

# Check if extractBranches exists
if (exists("extractBranches")) {
  cat("extractBranches found!\n")

  # Test different merge.min values
  for (merge_min in c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) {
    tryCatch({
      clusters_v215 <- extractBranches(
        ldna_raw,
        min.edges = min_edges_test,
        merge.min = merge_min,
        plot.tree = FALSE,
        cores = 1
      )

      n_clusters_v215 <- length(clusters_v215)
      cluster_sizes_v215 <- sapply(clusters_v215, length)

      cat("v2.15 (min.edges=", min_edges_test, ", merge.min=", merge_min, "):",
          n_clusters_v215, "clusters,",
          round(mean(cluster_sizes_v215), 1), "mean SNPs\n")

    }, error = function(e) {
      cat("  Error with merge.min=", merge_min, ":", e$message, "\n")
    })
  }
}

cat("\n=== Summary ===\n")
cat("v.65 uses extractClusters() with rm.COCs=TRUE to remove nested clusters\n")
cat("v2.15 uses extractBranches() with merge.min to control merging\n")
cat("Try adjusting merge.min in v2.15 to match v.65 behavior\n")

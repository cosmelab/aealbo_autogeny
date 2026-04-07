#!/usr/bin/env Rscript
# LDna Network Visualization - igraph-based SOC network plots
# Creates actual network visualizations showing LD connections between SNPs
# Author: Luciano Cosme
# Date: 2025-12-15

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(dplyr)
})

# Simple argument parsing (no optparse dependency)
args <- commandArgs(trailingOnly = TRUE)

# Default values
opt <- list(
  `ldna-dir` = "output/ldna",
  `r2-threshold` = 0.5,
  `max-snps` = 500,
  outdir = "output/ldna/network_plots"
)

# Parse arguments
for (arg in args) {
  if (grepl("^--ldna-dir=", arg)) opt$`ldna-dir` <- sub("^--ldna-dir=", "", arg)
  if (grepl("^--r2-threshold=", arg)) opt$`r2-threshold` <- as.numeric(sub("^--r2-threshold=", "", arg))
  if (grepl("^--max-snps=", arg)) opt$`max-snps` <- as.integer(sub("^--max-snps=", "", arg))
  if (grepl("^--outdir=", arg)) opt$outdir <- sub("^--outdir=", "", arg)
}

# Create output directory
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

cat("=== LDna Network Visualization ===\n")
cat("LDna directory:", opt$`ldna-dir`, "\n")
cat("R2 threshold:", opt$`r2-threshold`, "\n")
cat("Output directory:", opt$outdir, "\n\n")

# Populations and chromosomes
populations <- c("AUTO", "NON-AUTO-FIELD")
chromosomes <- c("chr1", "chr2", "chr3")

# Function to create network from LD matrix
create_ld_network <- function(ld_file, snp_list, r2_threshold = 0.5, max_snps = 500) {
  cat("  Loading LD matrix...\n")

  # Read compressed LD matrix
  ld_df <- read.table(gzfile(ld_file), header = TRUE, row.names = 1,
                      check.names = FALSE, na.strings = c("NA", "nan", "NaN"))

  # Subset if too large
  if (nrow(ld_df) > max_snps) {
    cat("    Subsetting to", max_snps, "SNPs\n")
    idx <- sample(seq_len(nrow(ld_df)), max_snps)
    ld_df <- ld_df[idx, idx]
  }

  # Convert to matrix
  ld_mat <- as.matrix(ld_df)
  ld_mat[is.na(ld_mat)] <- 0

  # Create adjacency matrix (edges where r2 >= threshold)
  adj_mat <- ld_mat >= r2_threshold
  diag(adj_mat) <- FALSE

  # Create igraph object
  g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)

  # Add edge weights (r2 values)
  edges <- get.edgelist(g)
  if (nrow(edges) > 0) {
    weights <- sapply(seq_len(nrow(edges)), function(i) {
      ld_mat[edges[i, 1], edges[i, 2]]
    })
    E(g)$weight <- weights
    E(g)$r2 <- weights
  }

  g
}

# Function to identify clusters in network
identify_network_clusters <- function(g) {
  # Use community detection
  clusters <- cluster_louvain(g, weights = E(g)$weight)
  V(g)$cluster <- membership(clusters)
  V(g)$cluster_size <- sizes(clusters)[membership(clusters)]

  g
}

# Function to plot network using base igraph (no ggraph dependency)
plot_ld_network <- function(g, title, output_file, max_edges = 5000) {
  if (ecount(g) == 0) {
    cat("    No edges to plot\n")
    return(NULL)
  }

  # Subsample edges if too many
  if (ecount(g) > max_edges) {
    cat("    Subsampling to", max_edges, "edges\n")
    edge_weights <- E(g)$weight
    keep_edges <- order(edge_weights, decreasing = TRUE)[1:max_edges]
    g <- subgraph.edges(g, keep_edges, delete.vertices = FALSE)
    # Remove isolated vertices
    g <- delete.vertices(g, degree(g) == 0)
  }

  # Set up colors by cluster
  n_clusters <- max(V(g)$cluster, na.rm = TRUE)
  cluster_colors <- rainbow(n_clusters)
  V(g)$color <- cluster_colors[V(g)$cluster]

  # Node sizes based on cluster size
  V(g)$size <- log1p(V(g)$cluster_size) * 2 + 2

  # Edge widths based on r2
  E(g)$width <- E(g)$r2 * 2

  # Save as PNG
  png(output_file, width = 1200, height = 1000, res = 100)
  par(mar = c(1, 1, 3, 1))

  # Plot using Fruchterman-Reingold layout
  layout <- layout_with_fr(g)
  plot(g, layout = layout,
       vertex.label = NA,
       edge.color = rgb(0.5, 0.5, 0.5, 0.3),
       main = paste0(title, "\n(", vcount(g), " SNPs, ", ecount(g), " edges)"))

  # Add legend
  legend("topright", legend = paste("Cluster", 1:min(n_clusters, 5)),
         col = cluster_colors[1:min(n_clusters, 5)], pch = 19, cex = 0.8)

  dev.off()
  cat("    Saved:", output_file, "\n")

  g
}

# Function to create summary network plot (all chromosomes)
plot_summary_network <- function(networks, pop, output_file) {
  # Combine all networks
  all_vertices <- data.frame()
  all_edges <- data.frame()

  for (chrom in names(networks)) {
    g <- networks[[chrom]]
    if (is.null(g) || vcount(g) == 0) next

    # Get vertex data
    v_df <- data.frame(
      name = V(g)$name,
      cluster = V(g)$cluster,
      cluster_size = V(g)$cluster_size,
      chromosome = chrom
    )
    all_vertices <- rbind(all_vertices, v_df)
  }

  if (nrow(all_vertices) == 0) {
    cat("  No data for summary plot\n")
    return(NULL)
  }

  # Create summary barplot of cluster sizes
  cluster_summary <- all_vertices |>
    group_by(chromosome, cluster) |>
    summarise(n_snps = n(), .groups = "drop") |>
    arrange(chromosome, desc(n_snps))

  p <- ggplot(cluster_summary, aes(x = reorder(paste(chromosome, cluster, sep = "_"), -n_snps),
                                    y = n_snps, fill = chromosome)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste0(pop, " - LD Network Cluster Sizes"),
         x = "Cluster (Chromosome_ID)",
         y = "Number of SNPs",
         fill = "Chromosome") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 6))

  ggsave(output_file, p, width = 10, height = 12, dpi = 150)
  cat("  Saved summary:", output_file, "\n")

  p
}

# Main analysis loop
for (pop in populations) {
  cat("\n=== Processing", pop, "===\n")

  networks <- list()

  for (chrom in chromosomes) {
    cat("\n  Chromosome:", chrom, "\n")

    # Check for LD matrix file
    ld_file <- file.path(opt$`ldna-dir`, "pop", chrom,
                         paste0(pop, ".", chrom, ".txt.gz"))

    if (!file.exists(ld_file)) {
      cat("    LD file not found:", ld_file, "\n")
      next
    }

    tryCatch({
      # Create network
      g <- create_ld_network(ld_file, NULL, opt$`r2-threshold`, opt$`max-snps`)

      cat("    Network: ", vcount(g), " nodes, ", ecount(g), " edges\n", sep = "")

      if (vcount(g) > 0 && ecount(g) > 0) {
        # Identify clusters
        g <- identify_network_clusters(g)
        n_clusters <- max(V(g)$cluster)
        cat("    Detected", n_clusters, "clusters\n")

        # Plot network
        output_file <- file.path(opt$outdir,
                                 paste0(pop, "_", chrom, "_network.png"))
        plot_ld_network(g,
                       paste0(pop, " - ", chrom, " LD Network"),
                       output_file)

        networks[[chrom]] <- g
      }
    }, error = function(e) {
      cat("    Error:", e$message, "\n")
    })
  }

  # Create summary plot
  if (length(networks) > 0) {
    summary_file <- file.path(opt$outdir,
                              paste0(pop, "_cluster_summary.png"))
    plot_summary_network(networks, pop, summary_file)
  }
}

# Create comparison plot
cat("\n=== Creating Population Comparison ===\n")

# Read all cluster position files
all_clusters <- data.frame()
for (pop in populations) {
  for (chrom in chromosomes) {
    cluster_file <- file.path(opt$`ldna-dir`, "clusters",
                              paste0(pop, "_", chrom, "_cluster_positions.txt"))
    if (file.exists(cluster_file)) {
      df <- read.table(cluster_file, header = TRUE, sep = "\t")
      df$Population <- pop
      df$Chromosome <- chrom
      all_clusters <- rbind(all_clusters, df)
    }
  }
}

if (nrow(all_clusters) > 0) {
  # Summary statistics
  summary_stats <- all_clusters |>
    group_by(Population, Chromosome) |>
    summarise(
      n_clusters = n(),
      total_snps = sum(n_snps),
      mean_size = mean(Size_Mb),
      max_size = max(Size_Mb),
      .groups = "drop"
    )

  cat("\nCluster Summary:\n")
  print(summary_stats)

  # Save summary
  write.csv(summary_stats,
            file.path(opt$outdir, "network_summary_stats.csv"),
            row.names = FALSE)
}

cat("\n=== Done! ===\n")
cat("Output directory:", opt$outdir, "\n")

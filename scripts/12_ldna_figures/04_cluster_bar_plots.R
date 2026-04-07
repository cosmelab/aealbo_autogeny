#!/usr/bin/env Rscript
# ==============================================================================
# 04_cluster_bar_plots.R
# Figure 3: Linkage cluster size bar plots
#
# EXTRACTED FROM: markdown/08.Linkage_network_analysis.Rmd
# Lines: 4205-4245
# ==============================================================================

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(scales)
})

# ==============================================================================
# Configuration
# ==============================================================================

# Input: Original data from dropbox (already computed)
input_file <- here("dropbox_original", "autogenous", "output", "ldna", "results_with_counts.rds")

# Output: Comparison directory
output_dir <- here("output", "figures_comparison", "replicated")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# Helper function for Mb labels
# ==============================================================================

label_mb <- function(x) {
  paste0(x / 1e6, "Mb")
}

# ==============================================================================
# Load Data
# ==============================================================================

cat("Loading cluster data...\n")

if (!file.exists(input_file)) {
  stop("Data file not found: ", input_file)
}

results_with_counts <- readRDS(input_file)

cat("Data loaded:\n")
cat("  Rows:", nrow(results_with_counts), "\n")
cat("  Columns:", paste(names(results_with_counts), collapse = ", "), "\n")
cat("  Populations:", paste(unique(results_with_counts$Population), collapse = ", "), "\n")
cat("  Chromosomes:", paste(unique(results_with_counts$Chromosome), collapse = ", "), "\n")

# ==============================================================================
# Prepare Data (EXACTLY as in original Rmd lines 4206-4219)
# ==============================================================================

cat("\nPreparing data for plotting...\n")

# Convert population labels to manuscript format
pop_labels <- c(
  "AUT" = "AUTO",
  "NEW" = "NON-AUTO-FIELD",
  "MAN" = "NON-AUTO"
)
results_with_counts$Population <- pop_labels[results_with_counts$Population]
cat("  Population labels updated to manuscript format\n")

# Check if Cluster_ID exists, if not create it
if (!"Cluster_ID" %in% names(results_with_counts)) {
  cluster_ids <- data.frame(
    Cluster    = sort(unique(results_with_counts$Cluster)),
    Cluster_ID = seq_along(unique(results_with_counts$Cluster))
  )
  results_with_counts <- merge(results_with_counts, cluster_ids, by = "Cluster")
}

# Rebuild the fill label: add ID before Cluster and drop " SNPs"
results_with_counts$Cluster_SNP <- with(
  results_with_counts,
  paste0(Cluster_ID, ": ", Cluster, " (", SNP_Count, ")")
)

# Reorder Cluster_SNP factor by Cluster_ID so the legend follows numeric order
results_with_counts$Cluster_SNP <- factor(
  results_with_counts$Cluster_SNP,
  levels = results_with_counts |>
    distinct(Cluster_ID, Cluster_SNP) |>
    arrange(Cluster_ID) |>
    pull(Cluster_SNP)
)

cat("  Clusters prepared:", length(unique(results_with_counts$Cluster_ID)), "\n")

# ==============================================================================
# Create Plot (EXACTLY as in original Rmd lines 4222-4245)
# ==============================================================================

cat("\nCreating cluster bar plot...\n")

p <- ggplot(results_with_counts, aes(xmin = Start, xmax = End, ymin = 0, ymax = Size)) +
  geom_rect(aes(fill = Cluster_SNP), color = "black", linewidth = 0.2) +
  geom_text(
    aes(x = (Start + End) / 2, y = Size + 0.5e6, label = Cluster_ID),
    size = 3, vjust = -0.2
  ) +
  scale_fill_discrete(name = "Cluster and SNP Count") +
  scale_x_continuous(labels = label_mb, breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(labels = label_mb, breaks = pretty_breaks(n = 5)) +
  labs(x = "Position", y = "Cluster Size") +
  theme_minimal() +
  theme(
    axis.text.y        = element_text(margin = margin(t = 0, r = 5)),
    axis.ticks.y       = element_blank(),
    panel.grid.major.x = element_line(color = "gray", linetype = "dotted"),
    strip.background   = element_rect(fill = "#e8e8e8", colour = NA),
    panel.grid.minor.x = element_blank(),
    legend.position    = "right",
    panel.spacing.x    = unit(1, "lines"),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA)
  ) +
  facet_grid(Population ~ Chromosome, scales = "fixed", space = "fixed")

# ==============================================================================
# Save Plot
# ==============================================================================

cat("\nSaving plots...\n")

# Save PDF
ggsave(
  file.path(output_dir, "figure_3_cluster_bars.pdf"),
  p, width = 12, height = 5, bg = "white"
)
cat("  Saved: figure_3_cluster_bars.pdf\n")

# Save PNG
ggsave(
  file.path(output_dir, "figure_3_cluster_bars.png"),
  p, width = 12, height = 5, dpi = 300, bg = "white"
)
cat("  Saved: figure_3_cluster_bars.png\n")

# ==============================================================================
# Summary
# ==============================================================================

cat("\n========================================\n")
cat("FIGURE 3 REPLICATION COMPLETE\n")
cat("========================================\n")
cat("Output directory:", output_dir, "\n")
cat("Files created:\n")
cat("  - figure_3_cluster_bars.pdf\n")
cat("  - figure_3_cluster_bars.png\n")
cat("========================================\n")

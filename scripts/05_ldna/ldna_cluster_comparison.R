#!/usr/bin/env Rscript
# LDna Cluster Comparison - Summary visualizations from existing cluster data
# Uses pre-computed cluster position files (no LD matrix parsing needed)
# Author: Luciano Cosme
# Date: 2025-12-15

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Simple argument parsing
args <- commandArgs(trailingOnly = TRUE)

# Default values
opt <- list(
  `ldna-dir` = "output/ldna",
  outdir = "output/ldna/network_plots"
)

# Parse arguments
for (arg in args) {
  if (grepl("^--ldna-dir=", arg)) opt$`ldna-dir` <- sub("^--ldna-dir=", "", arg)
  if (grepl("^--outdir=", arg)) opt$outdir <- sub("^--outdir=", "", arg)
}

# Create output directory
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

cat("=== LDna Cluster Comparison ===\n")
cat("LDna directory:", opt$`ldna-dir`, "\n")
cat("Output directory:", opt$outdir, "\n\n")

# Populations and chromosomes
populations <- c("AUTO", "NON-AUTO-FIELD")
chromosomes <- c("chr1", "chr2", "chr3")

# Read all cluster position files
cat("Loading cluster position files...\n")
all_clusters <- data.frame()

for (pop in populations) {
  for (chrom in chromosomes) {
    cluster_file <- file.path(opt$`ldna-dir`, "clusters",
                              paste0(pop, "_", chrom, "_cluster_positions.txt"))
    if (file.exists(cluster_file)) {
      df <- read.table(cluster_file, header = TRUE, sep = "\t")
      # Handle column name variations
      if (!"Population" %in% names(df)) df$Population <- pop
      if (!"Chromosome" %in% names(df)) df$Chromosome <- chrom
      all_clusters <- rbind(all_clusters, df)
      cat("  Loaded:", basename(cluster_file), "-", nrow(df), "clusters\n")
    } else {
      cat("  Not found:", cluster_file, "\n")
    }
  }
}

if (nrow(all_clusters) == 0) {
  cat("\nNo cluster data found. Exiting.\n")
  quit(save = "no", status = 1)
}

cat("\nTotal clusters loaded:", nrow(all_clusters), "\n")

# Standardize column names
# nSNPs -> n_snps
if ("nSNPs" %in% names(all_clusters) && !"n_snps" %in% names(all_clusters)) {
  names(all_clusters)[names(all_clusters) == "nSNPs"] <- "n_snps"
}
# Size (in bp) -> Size_Mb
if ("Size" %in% names(all_clusters) && !"Size_Mb" %in% names(all_clusters)) {
  all_clusters$Size_Mb <- all_clusters$Size / 1e6
}
# Ensure Chromosome is character (chr1, chr2, chr3 format)
if (is.numeric(all_clusters$Chromosome)) {
  all_clusters$Chromosome <- paste0("chr", all_clusters$Chromosome)
}

# Summary statistics
cat("\n=== Summary Statistics ===\n")
summary_stats <- all_clusters |>
  group_by(Population, Chromosome) |>
  summarise(
    n_clusters = n(),
    total_snps = sum(n_snps, na.rm = TRUE),
    mean_size_mb = mean(Size_Mb, na.rm = TRUE),
    max_size_mb = max(Size_Mb, na.rm = TRUE),
    median_snps = median(n_snps, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

# Save summary stats
write.csv(summary_stats,
          file.path(opt$outdir, "cluster_summary_stats.csv"),
          row.names = FALSE)
cat("\nSaved: cluster_summary_stats.csv\n")

# Population totals
pop_totals <- summary_stats |>
  group_by(Population) |>
  summarise(
    total_clusters = sum(n_clusters),
    total_snps = sum(total_snps),
    .groups = "drop"
  )
cat("\nPopulation Totals:\n")
print(pop_totals)

# === VISUALIZATION 1: Cluster count comparison ===
cat("\n=== Creating Visualizations ===\n")

p1 <- ggplot(summary_stats, aes(x = Chromosome, y = n_clusters, fill = Population)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "LDna Cluster Count by Population and Chromosome",
       x = "Chromosome",
       y = "Number of Clusters") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("AUTO" = "#E41A1C", "NON-AUTO-FIELD" = "#377EB8"))

ggsave(file.path(opt$outdir, "cluster_count_comparison.png"),
       p1, width = 8, height = 6, dpi = 150)
cat("Saved: cluster_count_comparison.png\n")

# === VISUALIZATION 2: Total SNPs in clusters ===
p2 <- ggplot(summary_stats, aes(x = Chromosome, y = total_snps, fill = Population)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Total SNPs in LDna Clusters",
       x = "Chromosome",
       y = "Total SNPs in Clusters") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("AUTO" = "#E41A1C", "NON-AUTO-FIELD" = "#377EB8"))

ggsave(file.path(opt$outdir, "total_snps_comparison.png"),
       p2, width = 8, height = 6, dpi = 150)
cat("Saved: total_snps_comparison.png\n")

# === VISUALIZATION 3: Cluster size distribution ===
p3 <- ggplot(all_clusters, aes(x = n_snps, fill = Population)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~Chromosome, scales = "free_y") +
  labs(title = "Distribution of Cluster Sizes (SNPs per cluster)",
       x = "Number of SNPs in Cluster",
       y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("AUTO" = "#E41A1C", "NON-AUTO-FIELD" = "#377EB8"))

ggsave(file.path(opt$outdir, "cluster_size_distribution.png"),
       p3, width = 12, height = 8, dpi = 150)
cat("Saved: cluster_size_distribution.png\n")

# === VISUALIZATION 4: Cluster genomic span ===
p4 <- ggplot(all_clusters, aes(x = Size_Mb, fill = Population)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~Chromosome, scales = "free_y") +
  labs(title = "Distribution of Cluster Genomic Span",
       x = "Cluster Size (Mb)",
       y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("AUTO" = "#E41A1C", "NON-AUTO-FIELD" = "#377EB8"))

ggsave(file.path(opt$outdir, "cluster_span_distribution.png"),
       p4, width = 12, height = 8, dpi = 150)
cat("Saved: cluster_span_distribution.png\n")

# === VISUALIZATION 5: Population summary bar chart ===
p5 <- ggplot(pop_totals, aes(x = Population, y = total_clusters, fill = Population)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = total_clusters), vjust = -0.3, size = 5) +
  labs(title = "Total LDna Clusters by Population",
       subtitle = "AUTO shows 12x more LD clustering than NON-AUTO-FIELD",
       x = "Population",
       y = "Total Number of Clusters") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("AUTO" = "#E41A1C", "NON-AUTO-FIELD" = "#377EB8"))

ggsave(file.path(opt$outdir, "population_total_clusters.png"),
       p5, width = 6, height = 6, dpi = 150)
cat("Saved: population_total_clusters.png\n")

# === VISUALIZATION 6: Scatter plot of cluster size vs span ===
p6 <- ggplot(all_clusters, aes(x = Size_Mb, y = n_snps, color = Population)) +
  geom_point(alpha = 0.6, size = 2) +
  facet_wrap(~Chromosome) +
  labs(title = "Cluster Size vs Genomic Span",
       x = "Genomic Span (Mb)",
       y = "Number of SNPs") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("AUTO" = "#E41A1C", "NON-AUTO-FIELD" = "#377EB8"))

ggsave(file.path(opt$outdir, "cluster_size_vs_span.png"),
       p6, width = 12, height = 8, dpi = 150)
cat("Saved: cluster_size_vs_span.png\n")

# === VISUALIZATION 7: Chromosome ideogram-style plot ===
# Show cluster positions along chromosomes
p7 <- ggplot(all_clusters, aes(x = Start/1e6, xend = End/1e6,
                                y = Population, yend = Population,
                                color = Population)) +
  geom_segment(linewidth = 3, alpha = 0.7) +
  facet_wrap(~Chromosome, ncol = 1, scales = "free_x") +
  labs(title = "LDna Cluster Positions Along Chromosomes",
       x = "Position (Mb)",
       y = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold")) +
  scale_color_manual(values = c("AUTO" = "#E41A1C", "NON-AUTO-FIELD" = "#377EB8"))

ggsave(file.path(opt$outdir, "cluster_positions_ideogram.png"),
       p7, width = 14, height = 10, dpi = 150)
cat("Saved: cluster_positions_ideogram.png\n")

cat("\n=== Done! ===\n")
cat("Output directory:", opt$outdir, "\n")
cat("Generated 7 visualization files + summary CSV\n")

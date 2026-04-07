#!/usr/bin/env Rscript
# ============================================================
# Process LDna Results - Extract Cluster SNPs and Create Plots
# ============================================================
#
# This script post-processes LDna results to:
#   1. Extract cluster SNP membership from saved rds files
#   2. Merge with bim file to get genomic positions
#   3. Calculate cluster Start/End/Size
#   4. Create visualizations matching original Rmd
#
# Based on: markdown/08.Linkage_network_analysis.Rmd
# Author: Luciano Cosme
# ============================================================

suppressPackageStartupMessages({
  .libPaths(c("/workspace/lib/ldna_v215", .libPaths()))
  library(tidyverse)
  library(here)
  library(data.table)
  library(LDna)
  library(reshape2)
  library(scales)
  library(viridis)
  library(ggrepel)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
project_dir <- ifelse(length(args) >= 1, args[1], "/project")
min_edges_target <- ifelse(length(args) >= 2, as.integer(args[2]), 100)
min_size_mb <- ifelse(length(args) >= 3, as.numeric(args[3]), 0.1)  # Minimum cluster size in Mb

cat("============================================================\n")
cat("LDna Results Processing and Visualization\n")
cat("============================================================\n")
cat("Project dir:", project_dir, "\n")
cat("Target min.edges:", min_edges_target, "\n")
cat("Minimum cluster size:", min_size_mb, "Mb\n")
cat("============================================================\n")

setwd(project_dir)

# Create output directory
dir.create("output/ldna/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("output/ldna/clusters", recursive = TRUE, showWarnings = FALSE)

# === Helper Functions ===

# Import bim file
import_bim <- function(file_path) {
  fread(file_path, header = FALSE,
        col.names = c("CHR", "SNP", "cM", "Position", "A1", "A2"))
}

# Extract cluster SNPs from extractBranches output
extract_cluster_snps <- function(clusters) {
  # extractBranches returns a list with cluster names
  # Each element contains SNP names
  if (is.null(clusters) || length(clusters) == 0) {
    return(data.frame())
  }

  # Get cluster names (removing the list structure)
  cluster_data <- list()

  for (cluster_name in names(clusters)) {
    snps <- clusters[[cluster_name]]
    if (length(snps) > 0) {
      cluster_data[[cluster_name]] <- data.frame(
        Cluster = cluster_name,
        SNP = snps,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(cluster_data) > 0) {
    return(bind_rows(cluster_data))
  }

  data.frame()
}

# Calculate cluster positions from SNP data
calculate_cluster_positions <- function(cluster_snps_df, bim_df) {
  # Merge cluster SNPs with bim to get positions
  merged <- cluster_snps_df |>
    left_join(bim_df |> select(SNP, CHR, Position), by = "SNP") |>
    na.omit()

  if (nrow(merged) == 0) {
    return(data.frame())
  }

  # Calculate cluster Start, End, Size
  cluster_positions <- merged |>
    group_by(Cluster) |>
    summarise(
      Start = min(Position),
      End = max(Position),
      Size = max(Position) - min(Position),
      nSNPs = n(),
      Chromosome = first(CHR),
      .groups = "drop"
    ) |>
    arrange(Start)

  cluster_positions
}

# === Main Processing ===

cat("\n=== Loading data ===\n")

# Read bim file for all SNPs
bim_all <- import_bim("output/ldna/pop/AUT.bim")
cat("Total SNPs in bim:", nrow(bim_all), "\n")

# Populations and chromosomes (AUT=AUTO, NEW=NON-AUTO-FIELD)
populations <- c("AUT", "NEW")
chromosomes <- 1:3

# Store all cluster data
all_clusters <- list()
all_positions <- list()

# === Process each population and chromosome ===
for (pop in populations) {
  for (chr in chromosomes) {
    cat("\n--- Processing", pop, "chr", chr, "---\n")

    chr_dir <- paste0("output/ldna/pop/chr", chr)
    ldna_file <- paste0(chr_dir, "/", pop, ".rds")
    ld2_file <- paste0(chr_dir, "/", pop, "_ld2.rds")

    if (!file.exists(ldna_file) || !file.exists(ld2_file)) {
      cat("WARNING: LDna files not found for", pop, "chr", chr, "\n")
      next
    }

    # Load LDna objects
    ldna <- readRDS(ldna_file)
    ld2 <- readRDS(ld2_file)

    cat("LD matrix:", nrow(ld2), "x", ncol(ld2), "\n")

    # Extract clusters at target min.edges
    tryCatch({
      clusters <- extractBranches(
        ldna,
        min.edges = min_edges_target,
        merge.min = 0.8,
        plot.tree = FALSE,
        cores = 1
      )

      # Extract SNP membership
      cluster_snps_df <- extract_cluster_snps(clusters)

      if (nrow(cluster_snps_df) > 0) {
        cluster_snps_df$Population <- pop
        cluster_snps_df$Chromosome <- chr

        cat("Extracted", length(unique(cluster_snps_df$Cluster)), "clusters with",
            nrow(cluster_snps_df), "SNPs\n")

        # Get bim for this chromosome
        chr_bim <- bim_all |> filter(CHR == chr)

        # Calculate positions
        cluster_positions <- calculate_cluster_positions(cluster_snps_df, chr_bim)
        cluster_positions$Population <- pop

        if (nrow(cluster_positions) > 0) {
          # Save cluster SNP data
          write.table(cluster_snps_df,
                      paste0("output/ldna/clusters/", pop, "_chr", chr, "_cluster_snps.txt"),
                      row.names = FALSE, sep = "\t", quote = FALSE)

          # Save cluster positions (sushi-style file)
          write.table(cluster_positions,
                      paste0("output/ldna/clusters/", pop, "_chr", chr, "_cluster_positions.txt"),
                      row.names = FALSE, sep = "\t", quote = FALSE)

          all_clusters[[paste(pop, chr, sep = "_chr")]] <- cluster_snps_df
          all_positions[[paste(pop, chr, sep = "_chr")]] <- cluster_positions
        }
      }
    }, error = function(e) {
      cat("ERROR:", e$message, "\n")
    })
  }
}

# Combine all cluster positions
if (length(all_positions) > 0) {
  combined_positions <- bind_rows(all_positions)
  # Remap to manuscript labels before saving and plotting
  combined_positions$Population <- pop_labels[combined_positions$Population]
  write.csv(combined_positions, "output/ldna/clusters/all_cluster_positions.csv", row.names = FALSE)
  cat("\n=== Combined", nrow(combined_positions), "clusters across all populations ===\n")
} else {
  cat("ERROR: No cluster data extracted!\n")
  quit(status = 1)
}

# === VISUALIZATIONS ===

cat("\n=== Creating visualizations ===\n")

# Format function for Mb labels
label_mb <- function(x) {
  sprintf("%.0fMb", x / 1e6)
}

# Manuscript label mapping and color palette
pop_labels  <- c("AUT" = "AUTO", "NEW" = "NON-AUTO-FIELD")
pop_colors  <- c("AUTO" = "#E41A1C", "NON-AUTO-FIELD" = "#377EB8")

# --- Plot 1: Individual chromosome plots per population ---
for (pop in populations) {
  for (chr in chromosomes) {
    pop_chr_data <- combined_positions |>
      filter(Population == pop, Chromosome == chr)

    if (nrow(pop_chr_data) == 0) next

    # Filter by minimum size
    pop_chr_data <- pop_chr_data |>
      filter(Size >= min_size_mb * 1e6)

    if (nrow(pop_chr_data) == 0) next

    # Convert to Mb
    pop_chr_data$Start_Mb <- pop_chr_data$Start / 1e6
    pop_chr_data$End_Mb <- pop_chr_data$End / 1e6
    pop_chr_data$Size_Mb <- pop_chr_data$Size / 1e6

    # Create plot (matching original style)
    p <- ggplot(pop_chr_data, aes(xmin = Start_Mb, xmax = End_Mb, ymin = 0, ymax = Size_Mb)) +
      geom_rect(aes(fill = as.factor(Cluster)), color = "black", linewidth = 0.2) +
      scale_x_continuous(labels = label_number(suffix = "Mb"), breaks = pretty_breaks(n = 10)) +
      scale_y_continuous(labels = label_number(suffix = "Mb"), breaks = pretty_breaks(n = 10)) +
      labs(
        x = paste0("Chromosome ", chr, " (Mb)"),
        y = "Cluster Size (Mb)",
        fill = "Cluster ID",
        title = paste(pop, "- Chromosome", chr)
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 5)),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(color = "gray", linetype = "dotted"),
        panel.grid.minor.x = element_blank(),
        legend.position = "right"
      ) +
      guides(fill = guide_legend(ncol = 2, title.position = "top"))

    ggsave(
      filename = paste0("output/ldna/plots/", pop, "_chr", chr, "_clusters.pdf"),
      plot = p,
      width = 10, height = 6
    )
    ggsave(
      filename = paste0("output/ldna/plots/", pop, "_chr", chr, "_clusters.png"),
      plot = p,
      width = 10, height = 6, dpi = 150
    )

    cat("  Saved:", pop, "chr", chr, "\n")
  }
}

# --- Plot 2: Faceted plot (Population x Chromosome) - Fixed scale ---
large_clusters <- combined_positions |>
  filter(Size >= 1e6)  # >= 1Mb

if (nrow(large_clusters) > 0) {
  p_facet_fixed <- ggplot(large_clusters, aes(xmin = Start, xmax = End, ymin = 0, ymax = Size)) +
    geom_rect(aes(fill = Population), color = "black", linewidth = 0.2, alpha = 0.7) +
    scale_x_continuous(labels = label_mb, breaks = pretty_breaks(n = 3)) +
    scale_y_continuous(labels = label_mb, breaks = pretty_breaks(n = 5)) +
    scale_fill_manual(values = pop_colors) +
    labs(x = "Position", y = "Cluster Size", title = "LD Clusters >= 1Mb (Fixed Scale)") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 5)),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray", linetype = "dotted"),
      strip.background = element_rect(fill = "#e8e8e8", colour = NA),
      panel.grid.minor.x = element_blank(),
      legend.position = "top",
      panel.spacing.x = unit(1, "lines")
    ) +
    facet_grid(Population ~ Chromosome, scales = "fixed", space = "fixed")

  ggsave("output/ldna/plots/all_populations_fixed_scale.pdf", p_facet_fixed, width = 12, height = 8)
  ggsave("output/ldna/plots/all_populations_fixed_scale.png", p_facet_fixed, width = 12, height = 8, dpi = 150)
  cat("  Saved: faceted plot (fixed scale)\n")
}

# --- Plot 3: Faceted plot - Free scale ---
if (nrow(large_clusters) > 0) {
  p_facet_free <- ggplot(large_clusters, aes(xmin = Start, xmax = End, ymin = 0, ymax = Size)) +
    geom_rect(aes(fill = Population), color = "black", linewidth = 0.2, alpha = 0.7) +
    scale_x_continuous(labels = label_mb, breaks = pretty_breaks(n = 3)) +
    scale_y_continuous(labels = label_mb, breaks = pretty_breaks(n = 5)) +
    scale_fill_manual(values = pop_colors) +
    labs(x = "Position", y = "Cluster Size", title = "LD Clusters >= 1Mb (Free Y Scale)") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 5)),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray", linetype = "dotted"),
      panel.grid.minor.x = element_blank(),
      legend.position = "top",
      panel.spacing.x = unit(1, "lines")
    ) +
    facet_wrap(~ Population + Chromosome, scales = "free_y", ncol = 3)

  ggsave("output/ldna/plots/all_populations_free_scale.pdf", p_facet_free, width = 12, height = 8)
  ggsave("output/ldna/plots/all_populations_free_scale.png", p_facet_free, width = 12, height = 8, dpi = 150)
  cat("  Saved: faceted plot (free scale)\n")
}

# --- Plot 4: Clusters 1-10 Mb range ---
mid_clusters <- combined_positions |>
  filter(Size >= 1e6, Size <= 10e6)

if (nrow(mid_clusters) > 0) {
  p_mid <- ggplot(mid_clusters, aes(xmin = Start, xmax = End, ymin = 0, ymax = Size)) +
    geom_rect(aes(fill = Population), color = "black", linewidth = 0.2, alpha = 0.7) +
    scale_x_continuous(labels = label_mb, breaks = pretty_breaks(n = 3)) +
    scale_y_continuous(labels = label_mb, breaks = pretty_breaks(n = 5)) +
    scale_fill_manual(values = pop_colors) +
    labs(x = "Position", y = "Cluster Size", title = "LD Clusters 1-10 Mb") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 5)),
      panel.grid.major.x = element_line(color = "gray", linetype = "dotted"),
      strip.background = element_rect(fill = "#e8e8e8", colour = NA),
      legend.position = "top",
      panel.spacing.x = unit(1, "lines")
    ) +
    facet_grid(Population ~ Chromosome, scales = "fixed", space = "fixed")

  ggsave("output/ldna/plots/clusters_1_to_10Mb.pdf", p_mid, width = 12, height = 8)
  ggsave("output/ldna/plots/clusters_1_to_10Mb.png", p_mid, width = 12, height = 8, dpi = 150)
  cat("  Saved: 1-10Mb clusters plot\n")
}

# --- Plot 5: Comparison bar chart ---
comparison_data <- combined_positions |>
  filter(Size >= 1e6) |>
  group_by(Population, Chromosome) |>
  summarise(
    n_clusters = n(),
    total_size_mb = sum(Size) / 1e6,
    .groups = "drop"
  )

if (nrow(comparison_data) > 0) {
  p_bar <- ggplot(comparison_data, aes(x = factor(Chromosome), y = n_clusters, fill = Population)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = pop_colors) +
    labs(
      x = "Chromosome",
      y = "Number of LD Clusters (>= 1Mb)",
      title = "LD Cluster Counts by Population and Chromosome"
    ) +
    theme_minimal() +
    theme(legend.position = "top")

  ggsave("output/ldna/plots/cluster_count_comparison.pdf", p_bar, width = 8, height = 6)
  ggsave("output/ldna/plots/cluster_count_comparison.png", p_bar, width = 8, height = 6, dpi = 150)
  cat("  Saved: cluster count comparison\n")
}

# === Load pcadapt outliers and overlay ===
cat("\n=== Loading selection scan outliers for overlay ===\n")

pcadapt_file <- "output/selection_scans/pcadapt/outlier_snps.txt"
outflank_file <- "output/selection_scans/outflank/outlier_snps.txt"

if (file.exists(pcadapt_file)) {
  pcadapt_snps <- fread(pcadapt_file, header = FALSE, col.names = c("SNP"))
  pcadapt_snps <- pcadapt_snps |>
    left_join(bim_all |> select(SNP, CHR, Position), by = "SNP") |>
    na.omit()
  cat("pcadapt outliers with positions:", nrow(pcadapt_snps), "\n")

  # Create overlay plot for each population
  for (pop in populations) {
    for (chr in chromosomes) {
      pop_chr_clusters <- combined_positions |>
        filter(Population == pop, Chromosome == chr, Size >= 0.5e6)

      chr_outliers <- pcadapt_snps |> filter(CHR == chr)

      if (nrow(pop_chr_clusters) == 0 || nrow(chr_outliers) == 0) next

      p_overlay <- ggplot() +
        # Cluster rectangles
        geom_rect(data = pop_chr_clusters,
                  aes(xmin = Start, xmax = End, ymin = 0, ymax = Size, fill = as.factor(Cluster)),
                  color = "black", linewidth = 0.2, alpha = 0.6) +
        # pcadapt SNPs
        geom_vline(data = chr_outliers,
                   aes(xintercept = Position),
                   color = "red", linetype = "dashed", linewidth = 0.5, alpha = 0.7) +
        scale_x_continuous(labels = label_mb, breaks = pretty_breaks(n = 8)) +
        scale_y_continuous(labels = label_mb, breaks = pretty_breaks(n = 5)) +
        labs(
          x = paste0("Chromosome ", chr, " (Mb)"),
          y = "Cluster Size (Mb)",
          title = paste(pop, "- Chr", chr, "with pcadapt outliers (red dashed)")
        ) +
        theme_minimal() +
        theme(legend.position = "none")

      ggsave(
        filename = paste0("output/ldna/plots/", pop, "_chr", chr, "_with_pcadapt.pdf"),
        plot = p_overlay, width = 12, height = 6
      )
      ggsave(
        filename = paste0("output/ldna/plots/", pop, "_chr", chr, "_with_pcadapt.png"),
        plot = p_overlay, width = 12, height = 6, dpi = 150
      )
      cat("  Saved:", pop, "chr", chr, "with pcadapt overlay\n")
    }
  }
}

# === Summary statistics ===
cat("\n=== SUMMARY STATISTICS ===\n")

summary_stats <- combined_positions |>
  group_by(Population) |>
  summarise(
    total_clusters = n(),
    clusters_1Mb = sum(Size >= 1e6),
    clusters_5Mb = sum(Size >= 5e6),
    clusters_10Mb = sum(Size >= 10e6),
    mean_size_kb = mean(Size) / 1000,
    max_size_mb = max(Size) / 1e6,
    total_coverage_mb = sum(Size) / 1e6,
    .groups = "drop"
  )

print(summary_stats)
write.csv(summary_stats, "output/ldna/plots/cluster_summary_stats.csv", row.names = FALSE)

cat("\n============================================================\n")
cat("VISUALIZATION COMPLETE\n")
cat("============================================================\n")
cat("Output directory: output/ldna/plots/\n")
cat("Cluster data: output/ldna/clusters/\n")

# Save session info
sink("output/ldna/plots/session_info.txt")
sessionInfo()
sink()

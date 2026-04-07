#!/usr/bin/env Rscript
# ==============================================================================
# Scaffold Cluster Plots - Extracted from Original 08.Linkage_network_analysis.Rmd
# ==============================================================================
# Source: dropbox_original/autogenous/scripts/markdown_files/08.Linkage_network_analysis.Rmd
# Lines: 4475-4800
#
# Creates:
#   - cluster_14_aut_scaffolds.pdf
#   - cluster_6_aut_scaffolds.pdf
#   - cluster_14_aut_scaffolds_with_small_windows.pdf
#   - cluster_6_aut_scaffolds_with_small_windows.pdf
#
# Author: Luciano Cosme
# ==============================================================================

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(scales)
  library(ggrepel)
})

# Set working directory
if (interactive()) {
  setwd(here::here())
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
project_dir <- ifelse(length(args) >= 1, args[1], here::here())
data_dir <- ifelse(length(args) >= 2, args[2], "dropbox_original/autogenous")

cat("==============================================================================\n")
cat("Scaffold Cluster Plots\n")
cat("==============================================================================\n")
cat("Project dir:", project_dir, "\n")
cat("Data dir:", data_dir, "\n")

# Create output directories
dir.create(file.path(project_dir, "output", "ldna", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_dir, "output", "figures_comparison", "replicated"), recursive = TRUE, showWarnings = FALSE)

# Output to comparison directory for side-by-side comparison
comparison_dir <- file.path(project_dir, "output", "figures_comparison", "replicated")

# ==============================================================================
# Import bim file function (from original Rmd)
# ==============================================================================
# Import bim file - MATCHES ORIGINAL import_bim.R exactly
# Column 1 is already the Scaffold ID (e.g., 1.1, 2.5, 3.2)
# DO NOT split or parse - keep as-is
import_bim <- function(file_path) {
  df <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE,
                   colClasses = c("character", "character", "integer", "numeric", "character", "character"),
                   col.names = c("Scaffold", "SNP", "Cm", "Position", "Allele1", "Allele2"))
  df
}

# ==============================================================================
# Load Data
# ==============================================================================
cat("\n--- Loading data ---\n")

# SNP data from bim file
bim_file <- file.path(data_dir, "output", "ldna", "files", "file1.bim")
if (!file.exists(bim_file)) {
  # Try alternative location
  bim_file <- file.path(data_dir, "output", "quality_control", "file4.bim")
}
if (!file.exists(bim_file)) {
  stop("Cannot find bim file. Tried: ", bim_file)
}

cat("Reading bim file:", bim_file, "\n")
snps <- import_bim(bim_file)
cat("Total SNPs:", nrow(snps), "\n")

# Load cluster data from LDna extractBranches output
aut_ch2_file <- file.path(data_dir, "output", "ldna", "pop", "chr2", "AUT_clusters_snps.rds")
aut_ch1_file <- file.path(data_dir, "output", "ldna", "pop", "chr1", "AUT_clusters_snps.rds")

if (!file.exists(aut_ch2_file) || !file.exists(aut_ch1_file)) {
  stop("Cannot find AUT cluster files. Expected: ", aut_ch2_file, " and ", aut_ch1_file)
}

cat("Loading cluster data...\n")
aut_ch2 <- readRDS(aut_ch2_file)
aut_ch1 <- readRDS(aut_ch1_file)

# Check available clusters
cat("\nAvailable clusters in aut_ch2:", names(aut_ch2), "\n")
cat("Available clusters in aut_ch1:", names(aut_ch1), "\n")

# ==============================================================================
# Extract Cluster 14 and Cluster 6 (from original Rmd lines 4475-4530)
# ==============================================================================
cat("\n--- Extracting clusters ---\n")

# Cluster 14: from chromosome 2, cluster 3139_0.78
cluster_14 <- snps[snps$SNP %in% aut_ch2$`3139_0.78`, ]
cat("Cluster 14: ", nrow(cluster_14), " SNPs\n")

# Cluster 6: from chromosome 1, cluster 1942_0.64
cluster_6 <- snps[snps$SNP %in% aut_ch1$`1942_0.64`, ]
cat("Cluster 6: ", nrow(cluster_6), " SNPs\n")

# Get all SNPs on the same scaffolds
cluster_14_all <- snps[snps$Scaffold %in% cluster_14$Scaffold, ]
cluster_6_all <- snps[snps$Scaffold %in% cluster_6$Scaffold, ]

cat("Cluster 14 scaffolds - total SNPs:", nrow(cluster_14_all), "\n")
cat("Cluster 6 scaffolds - total SNPs:", nrow(cluster_6_all), "\n")

# Add linked column
cluster_14_all <- cluster_14_all |>
  mutate(linked = SNP %in% aut_ch2$`3139_0.78`)

cluster_6_all <- cluster_6_all |>
  mutate(linked = SNP %in% aut_ch1$`1942_0.64`)

# Count linked SNPs
cat("\nCluster 14 linked counts:\n")
print(cluster_14_all |> group_by(linked) |> summarise(count = n()))

cat("\nCluster 6 linked counts:\n")
print(cluster_6_all |> group_by(linked) |> summarise(count = n()))

# ==============================================================================
# Load 157 outlier SNPs (from original Rmd line 4589)
# ==============================================================================
outlier_file <- file.path(data_dir, "output", "pcadapt", "outlier_157_SNPs.txt")
if (!file.exists(outlier_file)) {
  # Try data/validation location
  outlier_file <- file.path(project_dir, "data", "validation", "pcadapt", "outlier_157_SNPs_dropbox.txt")
}
if (!file.exists(outlier_file)) {
  cat("WARNING: 157 outlier SNPs file not found. Plots will not include outlier markers.\n")
  snps_157 <- data.frame(V1 = character(0))
} else {
  snps_157 <- read.table(outlier_file, stringsAsFactors = FALSE)
  cat("Loaded", nrow(snps_157), "outlier SNPs\n")
}

# Get 157 SNPs that are in each cluster
snps_157b <- cluster_14_all |> filter(SNP %in% snps_157$V1)
snps_157c <- cluster_6_all |> filter(SNP %in% snps_157$V1)

cat("Outlier SNPs in cluster 14:", nrow(snps_157b), "\n")
cat("Outlier SNPs in cluster 6:", nrow(snps_157c), "\n")

# ==============================================================================
# Helper Functions
# ==============================================================================
label_mb <- function(x) {
  sprintf("%.0fMb", x / 1e6)
}

# ==============================================================================
# Plot 1: Cluster 14 Scaffolds (from original Rmd lines 4620-4658)
# ==============================================================================
cat("\n--- Creating Cluster 14 scaffold plot ---\n")

rect_data <- cluster_14_all |>
  arrange(Scaffold, Position) |>
  mutate(change = linked != lag(linked, default = first(linked))) |>
  group_by(Scaffold) |>
  mutate(group_id = cumsum(change)) |>
  group_by(Scaffold, group_id, linked) |>
  summarize(start = min(Position), end = max(Position), .groups = "drop") |>
  ungroup()

p1 <- ggplot(
  rect_data,
  aes(xmin = start, xmax = end,
      ymin = as.numeric(Scaffold) - 0.1, ymax = as.numeric(Scaffold) + 0.1)
) +
  geom_rect(data = rect_data |> group_by(Scaffold) |>
            summarize(start = min(start), end = max(end), .groups = "drop"),
            aes(xmin = start, xmax = end, ymin = as.numeric(Scaffold) - 0.1, ymax = as.numeric(Scaffold) + 0.1),
            fill = "gray", inherit.aes = FALSE) +
  geom_rect(aes(fill = linked), data = rect_data |> filter(linked == TRUE)) +
  scale_fill_manual(values = c("TRUE" = "green")) +
  facet_wrap(~ Scaffold, scales = "fixed", ncol = 1) +
  geom_vline(data = snps_157b, aes(xintercept = Position), linetype = "solid", color = "red") +
  geom_text(data = snps_157b, aes(x = Position, y = as.numeric(Scaffold) + 0.2, label = SNP),
            inherit.aes = FALSE, angle = 90, vjust = 0, size = 2, check_overlap = TRUE) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(t = 1, b = 5))
  ) +
  scale_x_continuous(labels = label_mb, breaks = scales::pretty_breaks(n = 6)) +
  labs(x = "Position", title = "Continuous Stretches of Linked SNPs by Scaffold - Cluster 14 AUT", fill = "Linked") +
  guides(fill = "none")

# Save to both locations
output_path <- file.path(project_dir, "output", "ldna", "figures", "cluster_14_aut_scaffolds.pdf")
ggsave(output_path, p1, height = 5, width = 8, units = "in")
cat("Saved:", output_path, "\n")

# Also save to comparison directory
comparison_path <- file.path(comparison_dir, "cluster_14_aut_scaffolds.pdf")
ggsave(comparison_path, p1, height = 5, width = 8, units = "in")
cat("Saved:", comparison_path, "\n")

# ==============================================================================
# Plot 2: Cluster 14 with small windows (from original Rmd lines 4665-4704)
# ==============================================================================
cat("\n--- Creating Cluster 14 scaffold plot with small windows ---\n")

min_width <- 1e4

rect_data <- cluster_14_all |>
  arrange(Scaffold, Position) |>
  mutate(change = linked != lag(linked, default = first(linked))) |>
  group_by(Scaffold) |>
  mutate(group_id = cumsum(change)) |>
  group_by(Scaffold, group_id, linked) |>
  summarize(start = min(Position), end = max(Position), .groups = "drop") |>
  ungroup() |>
  mutate(width = end - start,
         adjusted_end = ifelse(width < min_width, start + min_width, end))

p2 <- ggplot(
  rect_data,
  aes(xmin = start, xmax = adjusted_end,
      ymin = as.numeric(Scaffold) - 0.1, ymax = as.numeric(Scaffold) + 0.1)
) +
  geom_rect(data = rect_data |> group_by(Scaffold) |>
            summarize(start = min(start), end = max(adjusted_end), .groups = "drop"),
            aes(xmin = start, xmax = end, ymin = as.numeric(Scaffold) - 0.1, ymax = as.numeric(Scaffold) + 0.1),
            fill = "gray", inherit.aes = FALSE) +
  geom_rect(aes(fill = linked), data = rect_data |> filter(linked == TRUE)) +
  scale_fill_manual(values = c("TRUE" = "green")) +
  facet_wrap(~ Scaffold, scales = "fixed", ncol = 1) +
  theme_minimal() +
  geom_vline(data = snps_157b, aes(xintercept = Position), linetype = "solid", color = "red") +
  geom_text_repel(data = snps_157b, aes(x = Position, y = as.numeric(Scaffold) + 0.2, label = SNP),
                  inherit.aes = FALSE, angle = 45, size = 2,
                  nudge_y = 0.1, direction = "y") +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12, face = "bold", hjust = 0.5)) +
  scale_x_continuous(labels = label_mb, breaks = scales::pretty_breaks(n = 6)) +
  labs(x = "Position", title = "Continuous Stretches of Linked SNPs by Scaffold - Cluster 14 AUT", fill = "Linked") +
  guides(fill = "none")

output_path <- file.path(project_dir, "output", "ldna", "figures", "cluster_14_aut_scaffolds_with_small_windows.pdf")
ggsave(output_path, p2, height = 5, width = 8, units = "in")
ggsave(
  file.path(comparison_dir, "cluster_14_aut_scaffolds_with_small_windows.pdf"),
  p2, height = 5, width = 8, units = "in"
)
cat("Saved:", output_path, "\n")

# ==============================================================================
# Plot 3: Cluster 6 Scaffolds (from original Rmd lines 4708-4754)
# ==============================================================================
cat("\n--- Creating Cluster 6 scaffold plot ---\n")

rect_data <- cluster_6_all |>
  arrange(Scaffold, Position) |>
  mutate(change = linked != lag(linked, default = first(linked))) |>
  group_by(Scaffold) |>
  mutate(group_id = cumsum(change)) |>
  group_by(Scaffold, group_id, linked) |>
  summarize(start = min(Position), end = max(Position), .groups = "drop") |>
  ungroup()

p3 <- ggplot(
  rect_data,
  aes(xmin = start, xmax = end,
      ymin = as.numeric(Scaffold) - 0.1, ymax = as.numeric(Scaffold) + 0.1)
) +
  geom_rect(data = rect_data |> group_by(Scaffold) |>
            summarize(start = min(start), end = max(end), .groups = "drop"),
            aes(xmin = start, xmax = end, ymin = as.numeric(Scaffold) - 0.1, ymax = as.numeric(Scaffold) + 0.1),
            fill = "gray", inherit.aes = FALSE) +
  geom_rect(aes(fill = linked), data = rect_data |> filter(linked == TRUE)) +
  scale_fill_manual(values = c("TRUE" = "green")) +
  facet_wrap(~ Scaffold, scales = "fixed", ncol = 1) +
  geom_vline(data = snps_157c, aes(xintercept = Position), linetype = "solid", color = "red") +
  geom_text(data = snps_157c, aes(x = Position, y = as.numeric(Scaffold) + 0.2, label = SNP),
            inherit.aes = FALSE, angle = 90, vjust = 0, size = 2, check_overlap = TRUE) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(t = 1, b = 5))
  ) +
  scale_x_continuous(labels = label_mb, breaks = scales::pretty_breaks(n = 6)) +
  labs(x = "Position", title = "Continuous Stretches of Linked SNPs by Scaffold - Cluster 6 AUT", fill = "Linked") +
  guides(fill = "none")

output_path <- file.path(project_dir, "output", "ldna", "figures", "cluster_6_aut_scaffolds.pdf")
ggsave(output_path, p3, height = 5, width = 8, units = "in")
ggsave(file.path(comparison_dir, "cluster_6_aut_scaffolds.pdf"), p3, height = 5, width = 8, units = "in")
cat("Saved:", output_path, "\n")

# ==============================================================================
# Plot 4: Cluster 6 with small windows (from original Rmd lines 4760-4800)
# ==============================================================================
cat("\n--- Creating Cluster 6 scaffold plot with small windows ---\n")

rect_data <- cluster_6_all |>
  arrange(Scaffold, Position) |>
  mutate(change = linked != lag(linked, default = first(linked))) |>
  group_by(Scaffold) |>
  mutate(group_id = cumsum(change)) |>
  group_by(Scaffold, group_id, linked) |>
  summarize(start = min(Position), end = max(Position), .groups = "drop") |>
  ungroup() |>
  mutate(width = end - start,
         adjusted_end = ifelse(width < min_width, start + min_width, end))

p4 <- ggplot(
  rect_data,
  aes(xmin = start, xmax = adjusted_end,
      ymin = as.numeric(Scaffold) - 0.1, ymax = as.numeric(Scaffold) + 0.1)
) +
  geom_rect(data = rect_data |> group_by(Scaffold) |>
            summarize(start = min(start), end = max(adjusted_end), .groups = "drop"),
            aes(xmin = start, xmax = end, ymin = as.numeric(Scaffold) - 0.1, ymax = as.numeric(Scaffold) + 0.1),
            fill = "gray", inherit.aes = FALSE) +
  geom_rect(aes(fill = linked), data = rect_data |> filter(linked == TRUE)) +
  scale_fill_manual(values = c("TRUE" = "green")) +
  facet_wrap(~ Scaffold, scales = "fixed", ncol = 1) +
  theme_minimal() +
  geom_vline(data = snps_157c, aes(xintercept = Position), linetype = "solid", color = "red") +
  geom_text_repel(data = snps_157c, aes(x = Position, y = as.numeric(Scaffold) + 0.2, label = SNP),
                  inherit.aes = FALSE, angle = 45, size = 2,
                  nudge_y = 0.1, direction = "y") +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12, face = "bold", hjust = 0.5)) +
  scale_x_continuous(labels = label_mb, breaks = scales::pretty_breaks(n = 6)) +
  labs(x = "Position", title = "Continuous Stretches of Linked SNPs by Scaffold - Cluster 6 AUT", fill = "Linked") +
  guides(fill = "none")

output_path <- file.path(project_dir, "output", "ldna", "figures", "cluster_6_aut_scaffolds_with_small_windows.pdf")
ggsave(output_path, p4, height = 5, width = 8, units = "in")
ggsave(
  file.path(comparison_dir, "cluster_6_aut_scaffolds_with_small_windows.pdf"),
  p4, height = 5, width = 8, units = "in"
)
cat("Saved:", output_path, "\n")

# ==============================================================================
# Save RDS files for reproducibility
# ==============================================================================
cat("\n--- Saving intermediate data ---\n")

saveRDS(cluster_14, file.path(project_dir, "output", "ldna", "cluster_14.rds"))
saveRDS(cluster_6, file.path(project_dir, "output", "ldna", "cluster_6.rds"))
saveRDS(cluster_14_all, file.path(project_dir, "output", "ldna", "cluster_14_all.rds"))
saveRDS(cluster_6_all, file.path(project_dir, "output", "ldna", "cluster_6_all.rds"))

cat("Saved RDS files\n")

# ==============================================================================
# Summary
# ==============================================================================
cat("\n==============================================================================\n")
cat("COMPLETE\n")
cat("==============================================================================\n")
cat("Output files:\n")
cat("  - output/ldna/figures/cluster_14_aut_scaffolds.pdf\n")
cat("  - output/ldna/figures/cluster_14_aut_scaffolds_with_small_windows.pdf\n")
cat("  - output/ldna/figures/cluster_6_aut_scaffolds.pdf\n")
cat("  - output/ldna/figures/cluster_6_aut_scaffolds_with_small_windows.pdf\n")
cat("  - output/ldna/cluster_14.rds\n")
cat("  - output/ldna/cluster_6.rds\n")
cat("  - output/ldna/cluster_14_all.rds\n")
cat("  - output/ldna/cluster_6_all.rds\n")

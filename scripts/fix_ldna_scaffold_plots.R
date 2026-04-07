#!/usr/bin/env Rscript
# Fix LDna scaffold-level cluster plots
# Recreates the original "Continuous Stretches of Linked SNPs by Scaffold" figures

library(tidyverse)
library(here)
library(ggrepel)
library(scales)

# Function to format numbers as Mb
label_mb <- function(x) {
  sprintf("%.0fMb", x / 1e6)
}

# Import BIM file function
import_bim <- function(file_path) {
  read.table(file_path, header = FALSE, stringsAsFactors = FALSE,
             col.names = c("Scaffold", "SNP", "Cm", "Position", "Allele1", "Allele2"))
}

# Load all SNPs with scaffold positions
cat("Loading SNP data from file1.bim...\n")
snps <- import_bim(here("data", "old_qc_reference", "file1.bim"))
cat("Total SNPs:", nrow(snps), "\n")

# Load cluster data from validation
cat("\nLoading cluster 14 data...\n")
cluster_14 <- readRDS(here("data", "validation", "ldna", "cluster_14.rds"))
cat("Cluster 14 SNPs:", nrow(cluster_14), "\n")
cat("Scaffolds in cluster 14:", paste(unique(cluster_14$Scaffold), collapse = ", "), "\n")

cat("\nLoading cluster 6 data...\n")
cluster_6 <- readRDS(here("data", "validation", "ldna", "cluster_6.rds"))
cat("Cluster 6 SNPs:", nrow(cluster_6), "\n")
cat("Scaffolds in cluster 6:", paste(unique(cluster_6$Scaffold), collapse = ", "), "\n")

# Load outlier SNPs
cat("\nLoading 157 outlier SNPs...\n")
snps_157 <- read.table(here("data", "validation", "pcadapt", "outlier_157_SNPs.txt"),
                       stringsAsFactors = FALSE)$V1
cat("Outlier SNPs loaded:", length(snps_157), "\n")

# Get all SNPs on cluster 14 scaffolds
cluster_14_all <- snps[snps$Scaffold %in% cluster_14$Scaffold, ]
cat("\nAll SNPs on Cluster 14 scaffolds:", nrow(cluster_14_all), "\n")

# Get all SNPs on cluster 6 scaffolds
cluster_6_all <- snps[snps$Scaffold %in% cluster_6$Scaffold, ]
cat("All SNPs on Cluster 6 scaffolds:", nrow(cluster_6_all), "\n")

# Mark linked SNPs for cluster 14
cluster_14_all <- cluster_14_all |>
  mutate(linked = SNP %in% cluster_14$SNP)
cat("Linked SNPs in Cluster 14:", sum(cluster_14_all$linked), "\n")

# Mark linked SNPs for cluster 6
cluster_6_all <- cluster_6_all |>
  mutate(linked = SNP %in% cluster_6$SNP)
cat("Linked SNPs in Cluster 6:", sum(cluster_6_all$linked), "\n")

# Get outlier SNPs in cluster 14
snps_157_c14 <- cluster_14_all |> filter(SNP %in% snps_157)
cat("\nOutlier SNPs in Cluster 14:", nrow(snps_157_c14), "\n")

# Get outlier SNPs in cluster 6
snps_157_c6 <- cluster_6_all |> filter(SNP %in% snps_157)
cat("Outlier SNPs in Cluster 6:", nrow(snps_157_c6), "\n")

# Create output directory
dir.create(here("output", "ldna", "figures"), showWarnings = FALSE, recursive = TRUE)

# ============================================
# Plot Cluster 14
# ============================================
cat("\n=== Generating Cluster 14 plot ===\n")

# Calculate the start and end positions for each stretch of TRUE or FALSE
rect_data_14 <- cluster_14_all |>
  arrange(Scaffold, Position) |>
  mutate(change = linked != lag(linked, default = first(linked))) |>
  group_by(Scaffold) |>
  mutate(group_id = cumsum(change)) |>
  group_by(Scaffold, group_id, linked) |>
  summarize(start = min(Position), end = max(Position), .groups = "drop") |>
  ungroup()

# Set a minimum width for visibility
min_width <- 1e4

rect_data_14 <- rect_data_14 |>
  mutate(width = end - start,
         adjusted_end = ifelse(width < min_width, start + min_width, end))

# Create the plot
p14 <- ggplot(rect_data_14, aes(xmin = start, xmax = adjusted_end,
                                 ymin = as.numeric(factor(Scaffold)) - 0.1,
                                 ymax = as.numeric(factor(Scaffold)) + 0.1)) +
  # Gray background for each scaffold
  geom_rect(data = rect_data_14 |> group_by(Scaffold) |>
              summarize(start = min(start), end = max(adjusted_end), .groups = "drop"),
            aes(xmin = start, xmax = end,
                ymin = as.numeric(factor(Scaffold)) - 0.1,
                ymax = as.numeric(factor(Scaffold)) + 0.1),
            fill = "gray", inherit.aes = FALSE) +
  # Green for linked SNPs
  geom_rect(aes(fill = linked), data = rect_data_14 |> filter(linked == TRUE)) +
  scale_fill_manual(values = c("TRUE" = "green")) +
  facet_wrap(~ Scaffold, scales = "fixed", ncol = 1) +
  # Red lines for outlier SNPs
  geom_vline(data = snps_157_c14, aes(xintercept = Position),
             linetype = "solid", color = "red") +
  # Labels for outlier SNPs
  geom_text_repel(data = snps_157_c14,
                  aes(x = Position, y = as.numeric(factor(Scaffold)) + 0.2, label = SNP),
                  inherit.aes = FALSE, angle = 45, size = 2,
                  nudge_y = 0.1, direction = "y") +
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
  labs(x = "Position",
       title = "Continuous Stretches of Linked SNPs by Scaffold - Cluster 14 AUT",
       fill = "Linked") +
  guides(fill = "none")

# Save cluster 14 plot
output_path_14 <- here("output", "ldna", "figures", "cluster_14_aut_scaffolds.png")
ggsave(output_path_14, p14, height = 7, width = 10, units = "in", dpi = 150)
cat("Saved:", output_path_14, "\n")

# Also save PDF
output_path_14_pdf <- here("output", "ldna", "figures", "cluster_14_aut_scaffolds.pdf")
ggsave(output_path_14_pdf, p14, height = 7, width = 10, units = "in")
cat("Saved:", output_path_14_pdf, "\n")

# ============================================
# Plot Cluster 6
# ============================================
cat("\n=== Generating Cluster 6 plot ===\n")

rect_data_6 <- cluster_6_all |>
  arrange(Scaffold, Position) |>
  mutate(change = linked != lag(linked, default = first(linked))) |>
  group_by(Scaffold) |>
  mutate(group_id = cumsum(change)) |>
  group_by(Scaffold, group_id, linked) |>
  summarize(start = min(Position), end = max(Position), .groups = "drop") |>
  ungroup()

rect_data_6 <- rect_data_6 |>
  mutate(width = end - start,
         adjusted_end = ifelse(width < min_width, start + min_width, end))

p6 <- ggplot(rect_data_6, aes(xmin = start, xmax = adjusted_end,
                               ymin = as.numeric(factor(Scaffold)) - 0.1,
                               ymax = as.numeric(factor(Scaffold)) + 0.1)) +
  geom_rect(data = rect_data_6 |> group_by(Scaffold) |>
              summarize(start = min(start), end = max(adjusted_end), .groups = "drop"),
            aes(xmin = start, xmax = end,
                ymin = as.numeric(factor(Scaffold)) - 0.1,
                ymax = as.numeric(factor(Scaffold)) + 0.1),
            fill = "gray", inherit.aes = FALSE) +
  geom_rect(aes(fill = linked), data = rect_data_6 |> filter(linked == TRUE)) +
  scale_fill_manual(values = c("TRUE" = "green")) +
  facet_wrap(~ Scaffold, scales = "fixed", ncol = 1) +
  geom_vline(data = snps_157_c6, aes(xintercept = Position),
             linetype = "solid", color = "red") +
  geom_text_repel(data = snps_157_c6,
                  aes(x = Position, y = as.numeric(factor(Scaffold)) + 0.2, label = SNP),
                  inherit.aes = FALSE, angle = 45, size = 2,
                  nudge_y = 0.1, direction = "y") +
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
  labs(x = "Position",
       title = "Continuous Stretches of Linked SNPs by Scaffold - Cluster 6 AUT",
       fill = "Linked") +
  guides(fill = "none")

output_path_6 <- here("output", "ldna", "figures", "cluster_6_aut_scaffolds.png")
ggsave(output_path_6, p6, height = 5, width = 10, units = "in", dpi = 150)
cat("Saved:", output_path_6, "\n")

output_path_6_pdf <- here("output", "ldna", "figures", "cluster_6_aut_scaffolds.pdf")
ggsave(output_path_6_pdf, p6, height = 5, width = 8, units = "in")
cat("Saved:", output_path_6_pdf, "\n")

cat("\n=== Done! ===\n")

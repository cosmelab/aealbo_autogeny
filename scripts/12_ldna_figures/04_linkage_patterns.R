#!/usr/bin/env Rscript
# ==============================================================================
# 04_linkage_patterns.R
# Figure 3: Linkage patterns (AUT vs NEW) - cluster size bars
#
# EXTRACTED FROM: dropbox_original/autogenous/scripts/markdown_files/08.Linkage_network_analysis.Rmd
# Section: Lines 3626-3725
# ==============================================================================

# Load libraries
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ==============================================================================
# Configuration
# ==============================================================================

# Input: Original data from dropbox
input_base <- here("dropbox_original", "autogenous", "output", "ldna", "pop")

# Output: Comparison directory
output_dir <- here("output", "figures_comparison", "replicated")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# Load Data (EXACTLY as in original Rmd lines 3628-3655)
# ==============================================================================

cat("Loading cluster data for NEW and AUT populations...\n")

# NEW
NEW1 <- readRDS(file.path(input_base, "chr1", "NEW_plot2.rds")) |> mutate(Population = "NEW")
NEW2 <- readRDS(file.path(input_base, "chr2", "NEW_plot2.rds")) |> mutate(Population = "NEW")
NEW3 <- readRDS(file.path(input_base, "chr3", "NEW_plot2.rds")) |> mutate(Population = "NEW")

# AUT
AUT1 <- readRDS(file.path(input_base, "chr1", "AUT_plot2.rds")) |> mutate(Population = "AUT")
AUT2 <- readRDS(file.path(input_base, "chr2", "AUT_plot2.rds")) |> mutate(Population = "AUT")
AUT3 <- readRDS(file.path(input_base, "chr3", "AUT_plot2.rds")) |> mutate(Population = "AUT")

# Bind the objects
albo <- rbind(NEW1, NEW2, NEW3, AUT1, AUT2, AUT3)

cat("Total rows:", nrow(albo), "\n")

# Create table with Size column
table1 <- albo |>
  mutate(Size = End - Start) |>
  dplyr::select(
    Population, Chromosome, Cluster, r2, nSegments, nSNPs, Start, End, Size
  ) |>
  arrange(Population, Chromosome, Start)

cat("Total clusters:", nrow(table1), "\n")

# ==============================================================================
# Filter for clusters >= 1Mb (lines 3660-3664)
# ==============================================================================

table1 <- table1 |>
  dplyr::filter(Size >= 1000000) |>
  dplyr::mutate(`Size (Mb)` = round(Size / 1000000, 2))

cat("Clusters >= 1Mb:", nrow(table1), "\n")

# ==============================================================================
# Create Plot (lines 3689-3725)
# ==============================================================================

cat("Creating linkage patterns plot...\n")

# To plot only clusters equal or bigger than 1Mb
albo2 <- table1 |> dplyr::filter(Size >= 1000000)

# Function to format numbers as Mb
label_mb <- function(x) {
  sprintf("%.0fMb", x / 1e6)
}

# Plot it (EXACTLY as in original)
p <- ggplot(albo2, aes(xmin = Start, xmax = End, ymin = 0, ymax = Size)) +
  geom_rect(aes(fill = as.factor(Cluster)), color = "black", linewidth = 0.2) +
  scale_x_continuous(labels = label_mb, breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(labels = label_mb, breaks = pretty_breaks(n = 5)) +
  labs(x = "Position", y = "Cluster Size") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 5)),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray", linetype = "dotted"),
    strip.background = element_rect(fill = "#e8e8e8", colour = NA),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",
    panel.spacing.x = unit(1, "lines")
  ) +
  facet_grid(Population ~ Chromosome, scales = "fixed", space = "fixed") +
  guides(fill = "none")

print(p)

# ==============================================================================
# Save Output
# ==============================================================================

output_path <- file.path(output_dir, "AUT_NEW_fixed.pdf")
ggsave(
  filename = output_path,
  plot = p,
  device = "pdf",
  width = 10,
  height = 5,
  units = "in"
)

cat("\nSaved:", output_path, "\n")
cat("Done!\n")

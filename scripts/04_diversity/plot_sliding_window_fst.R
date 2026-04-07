#!/usr/bin/env Rscript
# ==============================================================================
# Plot Sliding Window Fst Results
# Replicates: 10.Fst.Rmd (Section 3 - Sliding windows Fst estimates)
# ==============================================================================
# Input: output/fst/*.windowed.weir.fst
# Output: output/fst/figures/*.png
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(scales)
})

# Configuration
OUTPUT_DIR <- "output/fst"
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")

cat("============================================================\n")
cat("Plotting Sliding Window Fst Results\n")
cat("============================================================\n")

# Create figure directory
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Load data
# ==============================================================================
cat("\n=== Loading Fst data ===\n")

# Read the three comparison files
nonAuto_nonAutoField <- read.delim(
  file.path(OUTPUT_DIR, "nonAuto_nonAutoField.windowed.weir.fst"),
  header = TRUE
)
nonAuto_auto <- read.delim(
  file.path(OUTPUT_DIR, "nonAuto_auto.windowed.weir.fst"),
  header = TRUE
)
nonAutoField_auto <- read.delim(
  file.path(OUTPUT_DIR, "nonAutoField_auto.windowed.weir.fst"),
  header = TRUE
)

cat(sprintf("nonAuto_nonAutoField: %d windows\n", nrow(nonAuto_nonAutoField)))
cat(sprintf("nonAuto_auto: %d windows\n", nrow(nonAuto_auto)))
cat(sprintf("nonAutoField_auto: %d windows\n", nrow(nonAutoField_auto)))

# ==============================================================================
# Prepare combined data
# ==============================================================================
cat("\n=== Preparing combined data ===\n")

# Add comparison labels (using original naming for consistency with manuscript)
nonAuto_nonAutoField$comparison <- "MAN_NEW"
nonAuto_auto$comparison <- "MAN_AUT"
nonAutoField_auto$comparison <- "NEW_AUT"

# Combine
combined_data <- rbind(nonAuto_nonAutoField, nonAuto_auto, nonAutoField_auto)

# Calculate middle position
combined_data$mid_pos <- (combined_data$BIN_START + combined_data$BIN_END) / 2

# Reorder comparison factor
combined_data$comparison <- fct_relevel(
  combined_data$comparison,
  "MAN_NEW", "MAN_AUT", "NEW_AUT"
)

cat(sprintf("Total windows: %d\n", nrow(combined_data)))

# ==============================================================================
# Create main plot (Figure from 10.Fst.Rmd)
# ==============================================================================
cat("\n=== Creating sliding window Fst plot ===\n")

# Pastel colors for chromosomes
pastel_colors <- c("#ebd99f", "#B3CDE3", "#CCEBC5")

# Create the facet plot
window_fst <- ggplot(combined_data,
       aes(x = mid_pos, y = WEIGHTED_FST, color = as.factor(CHROM))) +
  geom_point(size = 0.1) +
  facet_grid(comparison ~ CHROM, scales = "free_x") +
  scale_x_continuous(labels = label_number(scale = 1e-6, suffix = "Mb")) +
  scale_color_manual(values = pastel_colors) +
  # Smooth line for MAN_NEW in black
  geom_smooth(
    data = subset(combined_data, comparison == "MAN_NEW"),
    aes(group = CHROM),
    se = FALSE,
    span = 0.3,
    color = "black"
  ) +
  # Smooth line for MAN_AUT in red
  geom_smooth(
    data = subset(combined_data, comparison == "MAN_AUT"),
    aes(group = CHROM),
    se = FALSE,
    span = 0.3,
    color = "red"
  ) +
  # Smooth line for NEW_AUT in red
  geom_smooth(
    data = subset(combined_data, comparison == "NEW_AUT"),
    aes(group = CHROM),
    se = FALSE,
    span = 0.3,
    color = "red"
  ) +
  theme_bw() +
  labs(x = "Position", y = "Weighted Fst", title = "") +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0.5),
    strip.text.y = element_text(angle = 90, hjust = 0.5),
    panel.spacing.x = unit(1.01, "lines"),
    legend.position = "none"
  )

# Save as PNG for HTML embedding
png_path <- file.path(FIGURE_DIR, "pairwise_estimates_windows.png")
ggsave(png_path, window_fst, height = 7, width = 9, dpi = 150)
cat(sprintf("Saved: %s\n", png_path))

# Save as PDF for publication
pdf_path <- file.path(FIGURE_DIR, "pairwise_estimates_windows.pdf")
ggsave(pdf_path, window_fst, height = 7, width = 9, dpi = 300)
cat(sprintf("Saved: %s\n", pdf_path))

# ==============================================================================
# Summary statistics
# ==============================================================================
cat("\n=== Summary Statistics ===\n")

summary_stats <- combined_data |>
  group_by(comparison) |>
  summarise(
    n_windows = n(),
    mean_fst = mean(WEIGHTED_FST, na.rm = TRUE),
    median_fst = median(WEIGHTED_FST, na.rm = TRUE),
    max_fst = max(WEIGHTED_FST, na.rm = TRUE),
    min_fst = min(WEIGHTED_FST, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

# Find high Fst windows (> 0.5)
high_fst <- combined_data |>
  filter(WEIGHTED_FST >= 0.5) |>
  arrange(desc(WEIGHTED_FST))

cat(sprintf("\nWindows with Fst >= 0.5: %d\n", nrow(high_fst)))
if (nrow(high_fst) > 0) {
  print(head(high_fst, 10))
}

# Save combined data for future use
saveRDS(combined_data, file.path(OUTPUT_DIR, "pop_fst_new.rds"))
cat(sprintf("\nSaved RDS: %s\n", file.path(OUTPUT_DIR, "pop_fst_new.rds")))

cat("\n============================================================\n")
cat("Done!\n")
cat("============================================================\n")

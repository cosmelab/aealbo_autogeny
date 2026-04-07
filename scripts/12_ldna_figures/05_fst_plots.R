#!/usr/bin/env Rscript
# ==============================================================================
# 05_fst_plots.R
# Figure 5: Fst sliding window (wild vs autogenous)
# Figure 6: Fst sliding window (males vs females)
#
# EXTRACTED FROM: dropbox_original/autogenous/scripts/markdown_files/10.Fst.Rmd
# Sections: Lines 592-720 (Figure 5) and Lines 892-1033 (Figure 6)
# ==============================================================================

# Load libraries
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(forcats)
})

# ==============================================================================
# Configuration
# ==============================================================================

# Input: Original data from dropbox
input_base <- here("dropbox_original", "autogenous", "output", "fst")

# Output: Comparison directory
output_dir <- here("output", "figures_comparison", "replicated")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# FIGURE 5: Pairwise Fst between populations (lines 592-720)
# ==============================================================================

cat("=== Creating Figure 5: Population Fst ===\n")

# Load data (using saved RDS for convenience)
combined_data <- readRDS(file = file.path(input_base, "pop_fst.rds"))

# Calculate the middle position
combined_data$mid_pos <- with(combined_data, (BIN_START + BIN_END) / 2)

# Reorder the levels of the comparison factor
combined_data$comparison <- fct_relevel(combined_data$comparison, "MAN_NEW", "MAN_AUT", "NEW_AUT")

# Define a set of pastel colors for three chromosomes
pastel_colors <- c("#ebd99f", "#B3CDE3", "#CCEBC5")

# Create the facet plot (EXACTLY as in original)
window_fst <- ggplot(combined_data,
       aes(
         x = mid_pos,
         y = WEIGHTED_FST,
         color = as.factor(CHROM)
       )) +
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
    strip.text.x = element_text(angle = 0, hjust = .5),
    strip.text.y = element_text(angle = 90, hjust = .5),
    panel.spacing.x = unit(1.01, "lines"),
    legend.position = "none"
  )

print(window_fst)

output_path <- file.path(output_dir, "pairwise_estimates_windows.pdf")
ggsave(output_path, window_fst, height = 7, width = 9, dpi = 300)
cat("Saved:", output_path, "\n")

# ==============================================================================
# FIGURE 6: Sex Fst (lines 892-1033)
# ==============================================================================

cat("\n=== Creating Figure 6: Sex Fst ===\n")

# Read individual files
new_sex <- read_delim(file.path(input_base, "new_sex.windowed.weir.fst"),
                      delim = "\t", col_names = TRUE, show_col_types = FALSE)
aut_sex <- read_delim(file.path(input_base, "aut_sex.windowed.weir.fst"),
                      delim = "\t", col_names = TRUE, show_col_types = FALSE)
new_aut_males <- read_delim(file.path(input_base, "new_aut_males.windowed.weir.fst"),
                            delim = "\t", col_names = TRUE, show_col_types = FALSE)
new_aut_females <- read_delim(file.path(input_base, "new_aut_females.windowed.weir.fst"),
                              delim = "\t", col_names = TRUE, show_col_types = FALSE)

# Add comparison labels
new_sex$comparison <- "NEW_sex"
aut_sex$comparison <- "AUT_sex"
new_aut_males$comparison <- "Males"
new_aut_females$comparison <- "Females"

# Combine the tibbles using rbind
combined_data <- rbind(new_sex, aut_sex, new_aut_males, new_aut_females)

# Calculate the middle position
combined_data$mid_pos <- with(combined_data, (BIN_START + BIN_END) / 2)

# Reorder the levels of the comparison factor
combined_data$comparison <- fct_relevel(combined_data$comparison, "NEW_sex", "AUT_sex", "Males", "Females")

# Create the facet plot (EXACTLY as in original)
window_fst <- ggplot(combined_data,
       aes(
         x = mid_pos,
         y = WEIGHTED_FST,
         color = as.factor(CHROM)
       )) +
  geom_point(size = 0.1) +
  facet_grid(comparison ~ CHROM, scales = "free_x") +
  scale_x_continuous(labels = label_number(scale = 1e-6, suffix = "Mb")) +
  scale_color_manual(values = pastel_colors) +
  geom_smooth(
    data = subset(combined_data, comparison == "NEW_sex"),
    aes(group = CHROM),
    se = FALSE,
    span = 0.3,
    color = "black"
  ) +
  geom_smooth(
    data = subset(combined_data, comparison == "AUT_sex"),
    aes(group = CHROM),
    se = FALSE,
    span = 0.3,
    color = "red"
  ) +
  geom_smooth(
    data = subset(combined_data, comparison == "Males"),
    aes(group = CHROM),
    se = FALSE,
    span = 0.3,
    color = "orange"
  ) +
  geom_smooth(
    data = subset(combined_data, comparison == "Females"),
    aes(group = CHROM),
    se = FALSE,
    span = 0.3,
    color = "orange"
  ) +
  theme_bw() +
  labs(x = "Position", y = "Weighted Fst", title = "") +
  theme(
    strip.text.x = element_text(angle = 0, hjust = .5),
    strip.text.y = element_text(angle = 90, hjust = .5),
    panel.spacing.x = unit(1.01, "lines"),
    legend.position = "none"
  )

print(window_fst)

output_path <- file.path(output_dir, "sex_pairwise_estimates_windows.pdf")
ggsave(output_path, window_fst, height = 7, width = 9, dpi = 300)
cat("Saved:", output_path, "\n")

cat("\nDone!\n")

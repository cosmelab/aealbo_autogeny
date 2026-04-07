#!/usr/bin/env Rscript
# ==============================================================================
# 03_genotype_frequencies.R
# Figure 2: Genotype frequencies for 17 outlier SNPs
#
# EXTRACTED FROM: dropbox_original/autogenous/scripts/markdown_files/07.Frequencies.Rmd
# Section: Lines 340-492
# ==============================================================================

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(ggtext)
  library(RColorBrewer)
  library(viridis)
})

# ==============================================================================
# Configuration
# ==============================================================================

# Input: Original data from dropbox
input_base <- here("dropbox_original", "autogenous", "output")

# Output: Comparison directory
output_dir <- here("output", "figures_comparison", "replicated")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# Define Theme (from my_theme2.R)
# ==============================================================================

my_theme <- function() {
  theme_minimal(base_size = 12, base_family = "") +
    theme(
      panel.grid.major = element_line(
        linetype = "dashed",
        linewidth = 0.2,
        color = "pink"
      ),
      panel.grid.minor = element_line(
        linetype = "dashed",
        linewidth = 0.2,
        color = "pink"
      ),
      axis.title.x = element_text(
        angle = 0,
        hjust = 1,
        face = "bold"
      ),
      axis.title.y = element_text(
        angle = 90,
        hjust = 1,
        face = "bold"
      )
    )
}

# ==============================================================================
# Load Data (EXACTLY as in original Rmd lines 112-142)
# ==============================================================================

cat("Loading allele frequency data...\n")

files <- list.files(
  path = file.path(input_base, "frequencies", "snp_files"),
  pattern = "*.txt",
  full.names = TRUE
)

all_data <- lapply(files, function(file) {
  dat <- read.table(
    file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    colClasses = c("character")
  )

  # Storing allele names from columns 3 and 4
  allele_names <- colnames(dat)[3:4]

  # Renaming columns 3 and 4
  colnames(dat)[3:4] <- c("Allele1_Value", "Allele2_Value")

  # Add the allele names as new columns
  dat$Allele1_Name <- allele_names[1]
  dat$Allele2_Name <- allele_names[2]

  # Add SNP column
  dat$SNP <- gsub(".txt", "", basename(file))

  dat
})

# Binding all data frames
all_data_df <- do.call(rbind, all_data)

cat("Total SNPs in data:", length(unique(all_data_df$SNP)), "\n")

# ==============================================================================
# Load Genotype Data (lines 239-338)
# ==============================================================================

cat("Loading genotype frequency data...\n")

geno_files <- list.files(
  path = file.path(input_base, "frequencies", "genotype_files"),
  pattern = "*.txt",
  full.names = TRUE
)

all_geno_data <- lapply(geno_files, function(file) {
  dat <- read.table(
    file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    colClasses = c("character")
  )

  # Get genotype names from columns 3, 4, and 5
  genotype_names <- colnames(dat)[3:5]

  # Rename the columns
  colnames(dat)[3:5] <- c("Genotype1_Value", "Genotype2_Value", "Genotype3_Value")

  # Add genotype names as new columns
  dat$Genotype1_Name <- genotype_names[1]
  dat$Genotype2_Name <- genotype_names[2]
  dat$Genotype3_Name <- genotype_names[3]

  # Add SNP column
  dat$SNP <- gsub(".txt", "", basename(file))

  dat
})

# Binding all data frames
all_data_df <- do.call(rbind, all_geno_data)

# ==============================================================================
# Select the 17 SNPs (lines 344-355)
# ==============================================================================

cat("Loading 17 outlier SNPs...\n")

snps_17 <- read.table(
  file.path(input_base, "pcadapt", "4_way_venn_common_SNPs_pcadapt_outflank.txt"),
  stringsAsFactors = FALSE
)

# Get the 17 SNPs
snps_aut17 <- all_data_df |>
  filter(SNP %in% snps_17$V1)

cat("Selected SNPs:", length(unique(snps_aut17$SNP)), "\n")

# ==============================================================================
# Set Colors (lines 391-405)
# ==============================================================================

# Generate a palette with brewer
palette1 <- brewer.pal(11, "Spectral")
palette2 <- brewer.pal(5, "Set3")
genotype_colors <- unique(c(palette1, palette2))
genotype_colors <- genotype_colors[1:16]

# ==============================================================================
# Create Plot (lines 409-492)
# ==============================================================================

cat("Creating genotype frequency plot...\n")

# Reshape the data to long format for plotting genotypes
geno_long <- snps_aut17 |>
  mutate(
    Genotype1_Name = as.character(Genotype1_Name),
    Genotype2_Name = as.character(Genotype2_Name),
    Genotype3_Name = as.character(Genotype3_Name)
  ) |>
  pivot_longer(
    cols = c(Genotype1_Value, Genotype2_Value, Genotype3_Value),
    names_to = "Genotype_Num",
    values_to = "Value"
  ) |>
  mutate(
    Genotype = case_when(
      Genotype_Num == "Genotype1_Value" ~ Genotype1_Name,
      Genotype_Num == "Genotype2_Value" ~ Genotype2_Name,
      Genotype_Num == "Genotype3_Value" ~ Genotype3_Name
    ),
    Value = as.numeric(Value)
  )

# Convert population labels to manuscript format
pop_labels <- c("AUT" = "AUTO", "MAN" = "NON-AUTO", "NEW" = "NON-AUTO-FIELD")
geno_long$Stratum <- pop_labels[geno_long$Stratum]

# Define the colors for the borders corresponding to each Stratum (manuscript labels)
stratum_border_colors <- c("AUTO" = "red", "NON-AUTO" = "black", "NON-AUTO-FIELD" = "black")

# Create a named vector of colors for the axis text
axis_text_colors <- c("AUTO" = "red", "NON-AUTO" = "black", "NON-AUTO-FIELD" = "black")

# Update the plot with genotypes
plot <- ggplot(geno_long, aes(x = Stratum, y = Value, fill = Genotype, group = Genotype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), aes(color = Stratum), linewidth = 0.5) +
  scale_fill_manual(name = "Genotype", values = genotype_colors) +
  scale_color_manual(name = "", values = stratum_border_colors) +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, by = 0.25)) +
  facet_wrap(~ SNP, scales = "free_x", ncol = 4) +
  geom_text(
    aes(label = Genotype, y = Value + 0.02, group = Genotype),
    position = position_dodge(width = 0.9),
    vjust = -0.25,
    size = 2.5,
    check_overlap = TRUE
  ) +
  labs(x = "Population", y = "Frequency") +
  my_theme() +
  theme(
    strip.text.x = element_text(size = 10, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
    legend.position = "top",
    legend.title.align = 0.5,
    strip.background = element_rect(fill = "#e8e8e8", colour = NA),
    panel.spacing = unit(1, "lines"),
    axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(color = "black")
  ) +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE, title = "Genotypes"),
    color = "none"
  )

# Print the plot
print(plot)

# ==============================================================================
# Save Output
# ==============================================================================

output_path <- file.path(output_dir, "significant_17_snps_genotypes.pdf")
ggsave(output_path, plot, height = 8, width = 7, dpi = 300)

cat("\nSaved:", output_path, "\n")
cat("Done!\n")

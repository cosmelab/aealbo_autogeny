#!/usr/bin/env Rscript
# ============================================================
# Genetic Diversity Analysis - He/Ho per Population
# ============================================================
#
# Purpose: Calculate observed and expected heterozygosity for
#          AUTO, NON-AUTO, and NON-AUTO-FIELD populations
#
# Requested by: Peter (Dec 7, 2025 email)
# Question: Is AUTO differentiation on PCA due to drift?
#
# NOTE: SNP chip data - values are conditional on ascertained
#       SNP set. Comparable BETWEEN populations but NOT to
#       sequence-based estimates from other studies.
#
# Author: Luciano Cosme
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_prefix <- ifelse(length(args) >= 1, args[1], "output/quality_control/file7")
output_dir <- ifelse(length(args) >= 2, args[2], "output/diversity")

# IMPORTANT: Use LD-pruned SNPs for heterozygosity calculation
pruned_snps_file <- "output/quality_control/indepSNP.prune.in"
use_pruned <- file.exists(pruned_snps_file)

cat("============================================================\n")
cat("Genetic Diversity Analysis (He/Ho)\n")
cat("============================================================\n")
cat("Input:", input_prefix, "\n")
cat("Output:", output_dir, "\n")
cat("============================================================\n")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# === STEP 1: Read PLINK files manually ===
cat("\n=== STEP 1: Reading PLINK files ===\n")

# Read .fam file for population info
fam_file <- paste0(input_prefix, ".fam")
fam <- fread(fam_file, header = FALSE,
             col.names = c("FID", "IID", "PID", "MID", "SEX", "PHENO"))
populations <- unique(fam$FID)
n_samples <- nrow(fam)

cat("Samples:", n_samples, "\n")
cat("Populations:", paste(populations, collapse = ", "), "\n")
cat("Population counts:\n")
print(table(fam$FID))

# Read .bim file for SNP info
bim_file <- paste0(input_prefix, ".bim")
bim <- fread(bim_file, header = FALSE,
             col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
n_snps_total <- nrow(bim)
cat("Total SNPs in file:", n_snps_total, "\n")

# Filter to LD-pruned SNPs if available
if (use_pruned) {
  pruned_snps <- fread(pruned_snps_file, header = FALSE)$V1
  cat("LD-pruned SNPs available:", length(pruned_snps), "\n")
  pruned_idx <- which(bim$SNP %in% pruned_snps)
  cat("Using", length(pruned_idx), "LD-pruned SNPs for He calculation\n")
} else {
  pruned_idx <- seq_len(n_snps_total)
  cat("WARNING: No LD-pruned SNP list found, using all SNPs\n")
}
n_snps <- length(pruned_idx)

# Read .bed file (binary genotypes)
bed_file <- paste0(input_prefix, ".bed")
cat("Reading BED file...\n")

# Read raw binary
bed_raw <- readBin(bed_file, what = "raw", n = file.info(bed_file)$size)

# Check magic number
if (bed_raw[1] != as.raw(0x6c) || bed_raw[2] != as.raw(0x1b) || bed_raw[3] != as.raw(0x01)) {
  stop("Invalid BED file format")
}

# Parse genotypes (SNP-major mode: bed_raw[3] == 0x01)
# Each byte encodes 4 samples, 2 bits each
# Coding: 00=hom1, 01=het, 10=missing, 11=hom2
cat("Parsing genotypes (only LD-pruned SNPs)...\n")

bytes_per_snp <- ceiling(n_samples / 4)
geno_matrix <- matrix(NA, nrow = n_snps, ncol = n_samples)

# Only read the pruned SNPs
for (i in seq_along(pruned_idx)) {
  snp_idx <- pruned_idx[i]
  start_byte <- 4 + (snp_idx - 1) * bytes_per_snp
  snp_bytes <- bed_raw[(start_byte):(start_byte + bytes_per_snp - 1)]

  sample_idx <- 1
  for (b in snp_bytes) {
    for (shift in c(0, 2, 4, 6)) {
      if (sample_idx > n_samples) break
      code <- as.integer(b) %/% (2^shift) %% 4
      # 0 = hom1 (0 copies of A1), 3 = hom2 (2 copies of A1), 2 = het (1 copy), 1 = missing
      geno_matrix[i, sample_idx] <- switch(code + 1, 0L, NA_integer_, 1L, 2L)
      sample_idx <- sample_idx + 1
    }
  }

  if (i %% 10000 == 0) cat("  Processed", i, "of", n_snps, "SNPs...\n")
}
cat("  Done! Processed", n_snps, "LD-pruned SNPs\n")

# === STEP 2: Calculate He/Ho per population ===
cat("\n=== STEP 2: Calculating heterozygosity ===\n")

# Function to calculate He and Ho
calc_diversity <- function(geno_vec) {
  # Remove missing
  geno_vec <- geno_vec[!is.na(geno_vec)]
  n <- length(geno_vec)
  if (n == 0) return(c(Ho = NA, He = NA, n = 0))

  # Count genotypes
  n_hom1 <- sum(geno_vec == 0)
  n_het <- sum(geno_vec == 1)
  n_hom2 <- sum(geno_vec == 2)

  # Observed heterozygosity
  Ho <- n_het / n

  # Allele frequency
  p <- (2 * n_hom2 + n_het) / (2 * n)
  q <- 1 - p

  # Expected heterozygosity (Hardy-Weinberg)
  He <- 2 * p * q

  c(Ho = Ho, He = He, n = n)
}

# Calculate per population
diversity_results <- list()

for (pop in populations) {
  cat("Processing population:", pop, "\n")
  pop_idx <- which(fam$FID == pop)
  pop_geno <- geno_matrix[, pop_idx, drop = FALSE]

  # Calculate per-locus
  locus_stats <- t(apply(pop_geno, 1, calc_diversity))

  # Mean across loci
  mean_Ho <- mean(locus_stats[, "Ho"], na.rm = TRUE)
  mean_He <- mean(locus_stats[, "He"], na.rm = TRUE)

  diversity_results[[pop]] <- list(
    Ho = mean_Ho,
    He = mean_He,
    n_samples = length(pop_idx),
    n_loci = sum(!is.na(locus_stats[, "Ho"])),
    locus_stats = locus_stats
  )

  cat("  Ho:", round(mean_Ho, 4), " He:", round(mean_He, 4), "\n")
}

# === STEP 3: Create summary table ===
cat("\n=== STEP 3: Summary ===\n")

# AUT=AUTO, MAN=NON-AUTO, NEW=NON-AUTO-FIELD (manuscript labels)
pop_labels <- c("AUT" = "AUTO", "MAN" = "NON-AUTO", "NEW" = "NON-AUTO-FIELD")

diversity_df <- data.frame(
  Population = pop_labels[names(diversity_results)],
  N_samples = sapply(diversity_results, function(x) x$n_samples),
  Ho = round(sapply(diversity_results, function(x) x$Ho), 4),
  He = round(sapply(diversity_results, function(x) x$He), 4),
  stringsAsFactors = FALSE
)
diversity_df$Fis <- round((diversity_df$He - diversity_df$Ho) / diversity_df$He, 4)

cat("\nPer-population diversity:\n")
print(diversity_df, row.names = FALSE)

# === STEP 4: Overall statistics ===
cat("\n=== STEP 4: Overall statistics ===\n")

# Overall Ho and He across all samples
overall_stats_per_locus <- t(apply(geno_matrix, 1, calc_diversity))
overall_Ho <- mean(overall_stats_per_locus[, "Ho"], na.rm = TRUE)
overall_He <- mean(overall_stats_per_locus[, "He"], na.rm = TRUE)

cat("Overall Ho:", round(overall_Ho, 4), "\n")
cat("Overall He:", round(overall_He, 4), "\n")

# === STEP 5: Save results ===
cat("\n=== STEP 5: Saving results ===\n")

# Population diversity
write_csv(diversity_df, file.path(output_dir, "population_diversity.csv"))
cat("Saved: population_diversity.csv\n")

# Overall stats
overall_df <- data.frame(
  Metric = c("Overall_Ho", "Overall_He", "Overall_Fis"),
  Value = round(c(overall_Ho, overall_He, (overall_He - overall_Ho) / overall_He), 4)
)
write_csv(overall_df, file.path(output_dir, "overall_statistics.csv"))
cat("Saved: overall_statistics.csv\n")

# === STEP 6: Create visualization ===
cat("\n=== STEP 6: Creating visualization ===\n")

plot_file <- file.path(output_dir, "diversity_barplot.png")
png(plot_file, width = 800, height = 600, res = 100)

diversity_long <- diversity_df |>
  pivot_longer(cols = c(Ho, He), names_to = "Metric", values_to = "Value")

diversity_long$Population <- factor(diversity_long$Population,
  levels = c("AUTO", "NON-AUTO", "NON-AUTO-FIELD"))

ggplot(diversity_long, aes(x = Population, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Ho" = "#3182ce", "He" = "#e53e3e"),
                    labels = c("Ho (Observed)", "He (Expected)")) +
  labs(title = "Expected Heterozygosity per Population",
       subtitle = sprintf("LD-pruned SNPs (n=%d)", n_snps),
       y = "Heterozygosity",
       x = "Population",
       caption = "AUT = AUTO, MAN = NON-AUTO, NEW = NON-AUTO-FIELD") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 11),
    legend.position = "top"
  ) +
  geom_text(aes(label = round(Value, 3)),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3)

dev.off()
cat("Saved:", plot_file, "\n")

# === STEP 7: Interpretation for Peter ===
cat("\n============================================================\n")
cat("INTERPRETATION FOR PETER\n")
cat("============================================================\n")

cat("\nQuestion: Is AUTO differentiation due to drift during selection?\n\n")

he_auto <- diversity_df$He[diversity_df$Population == "AUTO"]
he_nonauto <- diversity_df$He[diversity_df$Population == "NON-AUTO"]
he_field <- diversity_df$He[diversity_df$Population == "NON-AUTO-FIELD"]

cat("Expected Heterozygosity (He):\n")
cat("  AUTO:            ", he_auto, "\n")
cat("  NON-AUTO:        ", he_nonauto, "\n")
cat("  NON-AUTO-FIELD:  ", he_field, "\n\n")

if (he_auto < he_nonauto && he_auto < he_field) {
  cat("FINDING: AUTO has LOWER He than both control populations.\n")
  cat("This is CONSISTENT with genetic drift during selection:\n")
  cat("  - Smaller effective population size during selection\n")
  cat("  - Loss of genetic diversity through bottleneck\n")
  cat("  - Or stronger selection reducing diversity at linked loci\n")
} else if (he_auto > he_nonauto || he_auto > he_field) {
  cat("FINDING: AUTO does NOT have lower He than controls.\n")
  cat("This suggests differentiation is NOT primarily due to drift.\n")
} else {
  cat("FINDING: He is similar across populations.\n")
}

cat("\nCAVEAT: SNP chip data - relative comparisons only.\n")

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")

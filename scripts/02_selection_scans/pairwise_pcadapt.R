#!/usr/bin/env Rscript
# ============================================================
# Pairwise pcadapt Analysis - Selection Scans
# ============================================================
#
# Purpose: Run pcadapt for pairwise population comparisons
# Input: QC'd genotype data
# Output: Outlier SNPs for each pairwise comparison
#
# Comparisons:
#   1. NON-AUTO vs AUTO (original: 111 outliers)
#   2. NON-AUTO vs NON-AUTO-FIELD (original: 79 outliers)
#   3. NON-AUTO-FIELD vs AUTO (original: 67 outliers)
#
# Author: Luciano Cosme
# ============================================================

suppressPackageStartupMessages({
  library(pcadapt)
  library(data.table)
  library(tidyverse)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_dir <- ifelse(length(args) >= 1, args[1], "output/quality_control")
output_dir <- ifelse(length(args) >= 2, args[2], "output/pcadapt")

cat("============================================================\n")
cat("Pairwise pcadapt Analysis\n")
cat("============================================================\n")
cat("Input:", input_dir, "\n")
cat("Output:", output_dir, "\n")
cat("============================================================\n")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# === STEP 1: Load data ===
cat("\n=== STEP 1: Loading genotype data ===\n")

bed_file <- file.path(input_dir, "file7")

# Read .fam file for population info
fam <- fread(paste0(bed_file, ".fam"), header = FALSE)
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")
cat("Total samples:", nrow(fam), "\n")
cat("Populations:\n")
print(table(fam$FID))

# Get population indices (AUT=AUTO, MAN=NON-AUTO, NEW=NON-AUTO-FIELD)
pop_indices <- list(
  "MAN" = which(fam$FID == "MAN"),
  "AUT" = which(fam$FID == "AUT"),
  "NEW" = which(fam$FID == "NEW")
)

cat("AUT (AUTO) samples:", length(pop_indices$AUT), "\n")
cat("MAN (NON-AUTO) samples:", length(pop_indices$MAN), "\n")
cat("NEW (NON-AUTO-FIELD) samples:", length(pop_indices$NEW), "\n")

# Read .bim for SNP info
bim <- fread(paste0(bed_file, ".bim"), header = FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "POS", "A1", "A2")
cat("Total SNPs:", nrow(bim), "\n")

# IMPORTANT: Original pcadapt did NOT use LD-pruned SNPs!
# The original code used ALL SNPs from output/outflank/*.bed files
# pcadapt handles LD internally, so we should NOT pre-filter
cat("\n*** IMPORTANT: Using ALL SNPs (no LD pruning) as per original analysis ***\n")
cat("Original analysis used files like 'output/outflank/man_aut.bed' with ALL SNPs\n")
ld_pruned_snps <- NULL  # Set to NULL to use all SNPs

# === STEP 2: Define comparisons ===
cat("\n=== STEP 2: Defining pairwise comparisons ===\n")

# AUT=AUTO, MAN=NON-AUTO, NEW=NON-AUTO-FIELD
comparisons <- list(
  list(
    name = "MAN_vs_AUT",
    pops = c("MAN", "AUT"),
    original_file = "man_aut_SNPs_pcadapt.txt",
    expected_outliers = 111,
    K = 2
  ),
  list(
    name = "MAN_vs_NEW",
    pops = c("MAN", "NEW"),
    original_file = "man_new_SNPs_pcadapt.txt",
    expected_outliers = 79,
    K = 3
  ),
  list(
    name = "NEW_vs_AUT",
    pops = c("NEW", "AUT"),
    original_file = "new_aut_SNPs_pcadapt.txt",
    expected_outliers = 67,
    K = 2
  )
)

# === STEP 3: Run pairwise pcadapt ===
cat("\n=== STEP 3: Running pairwise pcadapt ===\n")

# Function to run pcadapt for a subset of samples
run_pairwise_pcadapt <- function(comparison, bed_file, fam, bim, ld_pruned_snps, output_dir) {
  cat("\n--- Running:", comparison$name, "---\n")

  # Get sample indices for this comparison
  pop1_idx <- which(fam$FID == comparison$pops[1])
  pop2_idx <- which(fam$FID == comparison$pops[2])
  sample_idx <- c(pop1_idx, pop2_idx)

  cat("Population 1 (", comparison$pops[1], "):", length(pop1_idx), "samples\n")
  cat("Population 2 (", comparison$pops[2], "):", length(pop2_idx), "samples\n")
  cat("Total samples:", length(sample_idx), "\n")

  # Create subset of data
  temp_dir <- file.path(output_dir, "temp")
  dir.create(temp_dir, showWarnings = FALSE)

  # Write sample list to file
  sample_file <- file.path(temp_dir, paste0(comparison$name, "_samples.txt"))
  subset_samples <- fam[sample_idx, .(FID, IID)]
  write.table(subset_samples, sample_file, row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Create subset BED file using PLINK
  subset_bed <- file.path(temp_dir, comparison$name)

  if (!is.null(ld_pruned_snps)) {
    # Write LD-pruned SNPs to file
    snp_file <- file.path(temp_dir, "ld_pruned_snps.txt")
    writeLines(ld_pruned_snps, snp_file)

    system(paste(
      "plink2",
      "--bfile", bed_file,
      "--keep", sample_file,
      "--extract", snp_file,
      "--make-bed",
      "--out", subset_bed,
      "--silent"
    ))
  } else {
    system(paste(
      "plink2",
      "--bfile", bed_file,
      "--keep", sample_file,
      "--make-bed",
      "--out", subset_bed,
      "--silent"
    ))
  }

  # Check subset was created
  if (!file.exists(paste0(subset_bed, ".bed"))) {
    cat("ERROR: Failed to create subset\n")
    return(NULL)
  }

  # Read subset bim
  subset_bim <- fread(paste0(subset_bed, ".bim"), header = FALSE)
  cat("SNPs in subset:", nrow(subset_bim), "\n")

  # Read with pcadapt
  geno <- read.pcadapt(paste0(subset_bed, ".bed"), type = "bed")
  cat("Loaded", attr(geno, "K"), "samples and", attr(geno, "L"), "SNPs\n")

  # Run pcadapt with appropriate K (varies by comparison!)
  K_value <- comparison$K
  cat("Running pcadapt with K=", K_value, ", min.maf=0.1...\n")
  pc_result <- pcadapt(geno, K = K_value, min.maf = 0.1)

  # Apply Benjamini-Hochberg correction
  pvalues <- pc_result$pvalues
  pvalues[is.na(pvalues)] <- 1
  adj_pvalues <- p.adjust(pvalues, method = "BH")

  # Get outliers at alpha = 0.05
  alpha <- 0.05
  outlier_idx <- which(adj_pvalues < alpha)
  outlier_snps <- subset_bim$V2[outlier_idx]

  cat("pcadapt outliers (BH adj p <", alpha, "):", length(outlier_snps), "\n")

  # Save results
  result_dir <- file.path(output_dir, comparison$name)
  dir.create(result_dir, showWarnings = FALSE)

  # Save outlier SNPs
  outlier_file <- file.path(result_dir, "outlier_snps.txt")
  writeLines(outlier_snps, outlier_file)
  cat("Saved outliers to:", outlier_file, "\n")

  # Save full results
  results_df <- data.frame(
    SNP = subset_bim$V2,
    CHR = subset_bim$V1,
    POS = subset_bim$V4,
    stat = pc_result$stat,
    pvalue = pc_result$pvalues,
    adj_pvalue = adj_pvalues,
    is_outlier = adj_pvalues < alpha,
    stringsAsFactors = FALSE
  )
  results_file <- file.path(result_dir, "pcadapt_results.csv")
  write_csv(results_df, results_file)

  # Create Manhattan plot
  plot_file <- file.path(result_dir, "manhattan_plot.png")
  png(plot_file, width = 1200, height = 600, res = 100)
  plot(pc_result, option = "manhattan")
  dev.off()
  cat("Saved Manhattan plot to:", plot_file, "\n")

  # Return results
  list(
    name = comparison$name,
    n_outliers = length(outlier_snps),
    outliers = outlier_snps
  )
}

# Run all comparisons
results_list <- list()
for (comp in comparisons) {
  result <- run_pairwise_pcadapt(comp, bed_file, fam, bim, ld_pruned_snps, output_dir)
  if (!is.null(result)) {
    results_list[[comp$name]] <- result
  }
}

# === STEP 4: Summary ===
cat("\n=== STEP 4: Summary ===\n")

summary_df <- data.frame(
  comparison = names(results_list),
  n_outliers = sapply(results_list, function(x) x$n_outliers)
)
print(summary_df)

summary_file <- file.path(output_dir, "pairwise_summary.csv")
write_csv(summary_df, summary_file)
cat("Saved summary to:", summary_file, "\n")

# Get union of all pairwise outliers
all_pairwise_outliers <- unique(unlist(lapply(results_list, function(x) x$outliers)))
cat("\nTotal unique outliers across all pairwise comparisons:", length(all_pairwise_outliers), "\n")

# Save union
union_file <- file.path(output_dir, "all_pairwise_outliers.txt")
writeLines(all_pairwise_outliers, union_file)
cat("Saved union to:", union_file, "\n")

# Clean up temp files
unlink(file.path(output_dir, "temp"), recursive = TRUE)

cat("\n============================================================\n")
cat("PAIRWISE PCADAPT ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("Output directory:", output_dir, "\n")
cat("\nResults:\n")
for (name in names(results_list)) {
  cat("  ", name, ":", results_list[[name]]$n_outliers, "outliers\n")
}
cat("============================================================\n")

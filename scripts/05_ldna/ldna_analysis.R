#!/usr/bin/env Rscript
# ============================================================
# LDna Network Analysis - Following Original Methodology
# ============================================================
#
# Based on: markdown/08.Linkage_network_analysis.Rmd
#
# Methodology (exact replication):
#   1. Start with QC'd data (NOT LD-pruned)
#   2. Create per-population PLINK files (AUTO, NON-AUTO-FIELD) with MAF 0.05
#   3. Get common SNPs shared between populations
#   4. Split by chromosome (1, 2, 3)
#   5. Calculate LD matrix per population per chromosome (plink --r2 triangle gz)
#   6. FORMAT LD matrices using bash (add header + row names)
#   7. Read formatted files into R
#   8. Run LDnaRaw() per population per chromosome
#   9. Run extractBranches() with min.edges=100
#
# NOTE: LDna v2.15 extractBranches() with min.edges=100 produces
#       identical results to original extractClusters() with rm.COCs=TRUE
#       (verified 2025-12-15: all cluster counts match exactly)
#
# Author: Luciano Cosme
# ============================================================

# Load LDna from project library (v2.15)
.libPaths(c("lib/ldna_v215", .libPaths()))

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(data.table)
  library(LDna)
  library(reshape2)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
project_dir <- ifelse(length(args) >= 1, args[1], "/project")
n_cores <- ifelse(length(args) >= 2, as.integer(args[2]), 6)

cat("============================================================\n")
cat("LDna Network Analysis (Original Methodology)\n")
cat("============================================================\n")
cat("Project dir:", project_dir, "\n")
cat("Cores:", n_cores, "\n")
cat("============================================================\n")

setwd(project_dir)

# === STEP 1: Prepare data ===
cat("\n=== STEP 1: Preparing data ===\n")

# Create output directories
dir.create("output/ldna/files", recursive = TRUE, showWarnings = FALSE)
dir.create("output/ldna/pop", recursive = TRUE, showWarnings = FALSE)

# Check if we need to start from scratch
if (!file.exists("output/ldna/files/file1.bed")) {
  cat("Creating initial file...\n")

  # Create file to remove sample 399 (failed heterozygosity)
  writeLines("AUT 399", "output/ldna/files/remove_aut_399.txt")

  # Use QC'd data (NOT LD-pruned - 110353 SNPs)
  system(paste(
    "plink2",
    "--allow-extra-chr",
    "--bfile output/quality_control/file7",
    "--remove output/ldna/files/remove_aut_399.txt",
    "--make-bed",
    "--out output/ldna/files/file1",
    "--silent"
  ))
}

# Check file1
system("grep 'samples\\|variants' output/ldna/files/file1.log 2>/dev/null || wc -l output/ldna/files/file1.bim")

# === STEP 2: Create per-population PLINK files ===
cat("\n=== STEP 2: Creating per-population files ===\n")

# Import bim
import_bim <- function(file_path) {
  fread(file_path, header = FALSE,
        col.names = c("CHR", "SNP", "cM", "Position", "A1", "A2"))
}

# Read fam file
fam <- fread("output/ldna/files/file1.fam", header = FALSE,
             col.names = c("FID", "IID", "PID", "MID", "SEX", "PHENO"))

# Use AUTO and NON-AUTO-FIELD (skip NON-AUTO - only 10 samples)
populations <- c("AUTO", "NON-AUTO-FIELD")
cat("Using populations:", paste(populations, collapse = ", "), "\n")

# Create population-specific files with MAF 0.05
for (pop in populations) {
  pop_fam <- fam |> filter(FID == pop)
  write.table(pop_fam, paste0("output/ldna/files/", pop, ".txt"),
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

  system(paste(
    "plink",
    "--keep-allele-order --allow-no-sex",
    "--bfile output/ldna/files/file1",
    "--make-bed",
    "--keep-fam", paste0("output/ldna/files/", pop, ".txt"),
    "--out", paste0("output/ldna/files/", pop),
    "--geno 0 --maf 0.05",
    "--silent"
  ))

  n_snps <- nrow(import_bim(paste0("output/ldna/files/", pop, ".bim")))
  cat(pop, ":", n_snps, "SNPs\n")
}

# === STEP 3: Get common SNPs ===
cat("\n=== STEP 3: Getting common SNPs ===\n")

pop_snps <- list()
for (pop in populations) {
  pop_snps[[pop]] <- import_bim(paste0("output/ldna/files/", pop, ".bim"))$SNP
}

common_snps <- Reduce(intersect, pop_snps)
cat("Common SNPs:", length(common_snps), "\n")
writeLines(common_snps, "output/ldna/files/common_snps.txt")

# === STEP 4: Filter to common SNPs ===
cat("\n=== STEP 4: Filtering to common SNPs ===\n")

for (pop in populations) {
  system(paste(
    "plink",
    "--keep-allele-order --allow-no-sex",
    "--bfile output/ldna/files/file1",
    "--make-bed",
    "--keep-fam", paste0("output/ldna/files/", pop, ".txt"),
    "--out", paste0("output/ldna/pop/", pop),
    "--extract output/ldna/files/common_snps.txt",
    "--silent"
  ))
  cat(pop, "- filtered\n")
}

# === STEP 5: Get chromosome SNP lists ===
cat("\n=== STEP 5: Getting chromosome-specific SNP lists ===\n")

pop_bim <- import_bim(paste0("output/ldna/pop/", populations[1], ".bim"))

for (chr in 1:3) {
  chr_snps <- pop_bim |> filter(CHR == chr) |> pull(SNP)
  writeLines(chr_snps, paste0("output/ldna/pop/chr", chr, "_snps.txt"))
  cat("Chr", chr, ":", length(chr_snps), "SNPs\n")

  # Create chromosome directory
  dir.create(paste0("output/ldna/pop/chr", chr), showWarnings = FALSE)
  file.copy(paste0("output/ldna/pop/chr", chr, "_snps.txt"),
            paste0("output/ldna/pop/chr", chr, "/"), overwrite = TRUE)
}

# === STEP 6: Calculate LD matrices ===
cat("\n=== STEP 6: Calculating LD matrices ===\n")

for (pop in populations) {
  for (chr in 1:3) {
    cat("Calculating LD:", pop, "chr", chr, "...\n")
    system(paste(
      "plink",
      "--allow-no-sex --keep-allele-order",
      "--bfile", paste0("output/ldna/pop/", pop),
      "--out", paste0("output/ldna/pop/chr", chr, "/", pop, ".chr", chr),
      "--r2 triangle gz",
      "--extract", paste0("output/ldna/pop/chr", chr, "_snps.txt"),
      "--silent"
    ))
  }
}

# === STEP 7: Format LD matrices using bash (like original) ===
cat("\n=== STEP 7: Formatting LD matrices ===\n")

for (chr in 1:3) {
  cat("Formatting chr", chr, "...\n")
  chr_dir <- paste0("output/ldna/pop/chr", chr)

  # Create header and row name files (exactly like original Rmd)
  # snps1.txt = list of SNPs
  # snps2.txt = empty line
  # snps3.txt.gz = empty line + SNPs (for row names column)
  # header.txt.gz = SNPs tab-separated (for column header)
  # IMPORTANT: Original uses tr '\n' ' ' then awk to convert to tabs (handles trailing whitespace)
  system(paste0(
    "cat ", chr_dir, "/../chr", chr, "_snps.txt | awk '{print $1}' > ", chr_dir, "/snps1.txt; ",
    "echo '' > ", chr_dir, "/snps2.txt; ",
    "cat ", chr_dir, "/snps2.txt ", chr_dir, "/snps1.txt | gzip -9 > ", chr_dir, "/snps3.txt.gz; ",
    "cat ", chr_dir, "/snps1.txt | tr '\\n' ' ' | awk -v OFS='\\t' '{$1=$1}1' | gzip -9 > ", chr_dir, "/header.txt.gz"
  ))

  # Add header and row names to each population's LD matrix (exactly like original)
  for (pop in populations) {
    ld_file <- paste0(chr_dir, "/", pop, ".chr", chr, ".ld.gz")
    if (file.exists(ld_file)) {
      # Step 1: Add header to LD matrix
      system(paste0(
        "cat ", chr_dir, "/header.txt.gz ", ld_file, " > ", chr_dir, "/", pop, ".chr", chr, ".ld.txt.gz"
      ))

      # Step 2: Rename
      system(paste0(
        "mv ", chr_dir, "/", pop, ".chr", chr, ".ld.txt.gz ", chr_dir, "/", pop, ".chr", chr, ".txt.gz"
      ))

      # Step 3: Add row names (paste SNP column - empty first line aligns with header)
      cmd3 <- paste0(
        "bash -c 'paste <(gzip -d < ", chr_dir, "/snps3.txt.gz)",
        " <(gzip -d < ", chr_dir, "/", pop, ".chr", chr, ".txt.gz)",
        " > ", chr_dir, "/", pop, ".chr", chr, ".txt'"
      )
      system(cmd3)

      cat("  Formatted:", pop, "chr", chr, "\n")
    }
  }
}

# === STEP 8: Run LDna for each population and chromosome ===
cat("\n=== STEP 8: Running LDna ===\n")

# Parameters from original analysis
# NOTE: v2.15 extractBranches() with min.edges=100 replicates original results
# Verified 2025-12-15: cluster counts match exactly for all pop/chr combinations
edges_primary <- 100  # Primary analysis - matches original extractClusters(rm.COCs=TRUE)
edges_sensitivity <- c(20, 50, 100, 150, 200)  # For sensitivity analysis

run_ldna_analysis <- function(pop, chr, n_cores) {
  cat("\n--- LDna for", pop, "chr", chr, "---\n")

  chr_dir <- paste0("output/ldna/pop/chr", chr)
  txt_file <- paste0(chr_dir, "/", pop, ".chr", chr, ".txt")

  if (!file.exists(txt_file)) {
    cat("WARNING: Formatted LD file not found:", txt_file, "\n")
    return(NULL)
  }

  # Read formatted LD matrix (exactly like original - no fill parameter)
  cat("Reading LD matrix...\n")
  ld1 <- read.delim(txt_file, sep = "\t", header = TRUE, row.names = 1,
                    stringsAsFactors = FALSE, check.names = FALSE,
                    na.strings = c("NA", "nan", "NaN", ""))

  cat("LD matrix:", nrow(ld1), "x", ncol(ld1), "\n")

  # Convert to matrix (like original)
  ld2 <- as.matrix(sapply(ld1, as.numeric))
  names <- rownames(ld1)
  rownames(ld2) <- names

  # Check for problematic values
  n_na <- sum(is.na(ld2))
  n_nan <- sum(is.nan(ld2))
  n_inf <- sum(is.infinite(ld2))
  cat("Before processing - NA:", n_na, "NaN:", n_nan, "Inf:", n_inf, "\n")

  # Replace NaN and Inf with NA
  ld2[is.nan(ld2)] <- NA
  ld2[is.infinite(ld2)] <- NA

  # Remove SNPs with ALL NA values (no valid LD with any SNP)
  # These cause hclust to fail
  row_valid <- apply(ld2, 1, function(x) sum(!is.na(x)) > 0)
  col_valid <- apply(ld2, 2, function(x) sum(!is.na(x)) > 0)
  both_valid <- row_valid & col_valid
  n_removed <- sum(!both_valid)

  if (n_removed > 0) {
    cat("Removing", n_removed, "SNPs with all NA values\n")
    ld2 <- ld2[both_valid, both_valid]
    names <- names[both_valid]
    cat("Remaining SNPs:", nrow(ld2), "\n")
  }

  # Set diagonal to NA (like original)
  diag(ld2) <- NA

  # Keep only lower triangle (like original)
  ld2[!lower.tri(ld2)] <- NA

  # Check for problematic values after processing
  n_na_after <- sum(is.na(ld2))
  n_valid <- sum(!is.na(ld2))
  cat("After processing - NA:", n_na_after, "Valid values:", n_valid, "\n")

  # LDna requires at least some valid LD values
  if (n_valid == 0) {
    cat("ERROR: No valid LD values in matrix!\n")
    return(NULL)
  }

  # Run LDnaRaw
  cat("Running LDnaRaw...\n")
  ldna <- LDnaRaw(ld2, mc.cores = n_cores, method = "single")

  # Save objects
  saveRDS(ldna, paste0(chr_dir, "/", pop, "_ldna.rds"))
  saveRDS(ld2, paste0(chr_dir, "/", pop, "_ld2.rds"))
  cat("Saved LDna objects\n")

  # Create PDF for cluster plots
  pdf(paste0(chr_dir, "/", pop, "_cluster_plots.pdf"), width = 20, height = 12)
  par(mfcol = c(1, 4))

  all_summaries <- list()

  # Run primary analysis (min.edges=100 - matches original)
  cat("  Primary analysis: min.edges =", edges_primary, "\n")

  tryCatch({
    clusters_primary <- extractBranches(
      ldna,
      min.edges = edges_primary,
      plot.tree = TRUE,
      cores = 1
    )

    summary_primary <- summaryLDna(ldna, clusters_primary, ld2)

    if (!is.null(summary_primary) && nrow(summary_primary) > 0) {
      # Save primary results
      write.table(summary_primary,
                  paste0(chr_dir, "/summary_", pop, "_primary.txt"),
                  row.names = FALSE, sep = "\t", quote = FALSE)
      saveRDS(clusters_primary, paste0(chr_dir, "/", pop, "_clusters.rds"))

      summary_primary$min_edges <- edges_primary
      summary_primary$population <- pop
      summary_primary$chromosome <- chr
      all_summaries[["primary"]] <- summary_primary

      cat("    PRIMARY: Found", nrow(summary_primary), "clusters\n")
    }
  }, error = function(e) {
    cat("    Error in primary analysis:", e$message, "\n")
  })

  # Run sensitivity analysis
  cat("  Sensitivity analysis...\n")
  for (edges in edges_sensitivity) {
    tryCatch({
      clusters <- extractBranches(
        ldna,
        min.edges = edges,
        plot.tree = FALSE,
        cores = 1
      )

      summary_df <- summaryLDna(ldna, clusters, ld2)

      if (!is.null(summary_df) && nrow(summary_df) > 0) {
        write.table(summary_df,
                    paste0(chr_dir, "/summary_", pop, "_edges", edges, ".txt"),
                    row.names = FALSE, sep = "\t", quote = FALSE)

        summary_df$min_edges <- edges
        summary_df$population <- pop
        summary_df$chromosome <- chr
        all_summaries[[as.character(edges)]] <- summary_df

        cat("    min.edges=", edges, ": ", nrow(summary_df), " clusters\n", sep = "")
      }
    }, error = function(e) {
      cat("    Error at edges=", edges, ": ", e$message, "\n", sep = "")
    })
  }

  dev.off()

  if (length(all_summaries) > 0) {
    combined <- bind_rows(all_summaries)
    write.csv(combined, paste0(chr_dir, "/", pop, "_all_summaries.csv"), row.names = FALSE)
    return(combined)
  }

  NULL
}

# Run for each population and chromosome
all_results <- list()

for (pop in populations) {
  for (chr in 1:3) {
    cat("\n========================================\n")
    cat("Processing:", pop, "chromosome", chr, "\n")
    cat("========================================\n")

    result <- run_ldna_analysis(pop, chr, n_cores)

    if (!is.null(result)) {
      all_results[[paste(pop, chr, sep = "_chr")]] <- result
    }
  }
}

# === STEP 9: Save combined results ===
cat("\n=== STEP 9: Saving combined results ===\n")

if (length(all_results) > 0) {
  combined_results <- bind_rows(all_results)
  write.csv(combined_results, "output/ldna/all_ldna_results.csv", row.names = FALSE)
  cat("Saved combined results\n")

  cat("\n=== SUMMARY ===\n")
  print(combined_results |>
          group_by(population, chromosome, min_edges) |>
          summarise(n_clusters = n(), total_snps = sum(nSNPs), .groups = "drop"))
}

cat("\n============================================================\n")
cat("LDna ANALYSIS COMPLETE\n")
cat("============================================================\n")

# Save session info
sink("output/ldna/session_info.txt")
sessionInfo()
sink()

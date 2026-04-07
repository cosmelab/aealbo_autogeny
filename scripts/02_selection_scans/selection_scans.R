#!/usr/bin/env Rscript
# ==============================================================================
# Selection Scan Analysis - REPLICATING ORIGINAL ANALYSIS
# ==============================================================================
# This script replicates the exact original analysis workflow:
# 1. Extract intergenic SNPs
# 2. LD prune
# 3. Run OutFLANK with neutral calibration
# 4. Run pcadapt on LD-pruned set
# ==============================================================================

suppressPackageStartupMessages({
  library(OutFLANK)
  library(pcadapt)
  library(qvalue)
  library(vcfR)
  library(tidyverse)
  library(ggplot2)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
bed_file <- ifelse(length(args) >= 1, args[1], "output/quality_control/file7")
output_dir <- ifelse(length(args) >= 2, args[2], "output/selection_scans")
intergenic_file <- ifelse(length(args) >= 3, args[3], "data/files/intergenic_SNPs.txt")

cat("============================================================\n")
cat("Selection Scan - ORIGINAL WORKFLOW REPLICATION\n")
cat("============================================================\n")
cat("Input:", bed_file, "\n")
cat("Output:", output_dir, "\n")
cat("============================================================\n\n")

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "outflank"), showWarnings = FALSE)
dir.create(file.path(output_dir, "pcadapt"), showWarnings = FALSE)

# =============================================================================
# STEP 1: Extract intergenic SNPs and LD prune
# =============================================================================
cat("=== STEP 1: Extract intergenic SNPs and LD prune ===\n")

# Check if intergenic file exists
if (!file.exists(intergenic_file)) {
  stop("Intergenic SNPs file not found: ", intergenic_file)
}

intergenic_dir <- file.path(output_dir, "intermediate")
dir.create(intergenic_dir, showWarnings = FALSE)

# Step 1a: Extract intergenic SNPs
cat("Extracting intergenic SNPs...\n")
system(paste(
  "plink2",
  "--bfile", bed_file,
  "--extract", intergenic_file,
  "--make-bed",
  "--out", file.path(intergenic_dir, "intergenic"),
  "--silent"
))

# Step 1b: LD pruning on intergenic SNPs
cat("LD pruning intergenic SNPs (window=5, step=1, r2=0.1)...\n")
system(paste(
  "plink2",
  "--bfile", file.path(intergenic_dir, "intergenic"),
  "--indep-pairwise 5 1 0.1",
  "--out", file.path(intergenic_dir, "indepSNP"),
  "--silent"
))

# Step 1c: Create LD-pruned intergenic dataset
cat("Creating LD-pruned intergenic dataset...\n")
system(paste(
  "plink2",
  "--bfile", file.path(intergenic_dir, "intergenic"),
  "--extract", file.path(intergenic_dir, "indepSNP.prune.in"),
  "--make-bed",
  "--export vcf",
  "--out", file.path(intergenic_dir, "intergenic_ldpruned"),
  "--silent"
))

# Count SNPs
ld_snps <- readLines(file.path(intergenic_dir, "indepSNP.prune.in"))
cat("Intergenic SNPs after LD pruning:", length(ld_snps), "\n\n")

# =============================================================================
# STEP 2: OutFLANK Analysis
# =============================================================================
cat("=== STEP 2: OutFLANK Analysis ===\n")

# --- STEP 2a: Calibrate neutral distribution using INTERGENIC SNPs ---
cat("Step 2a: Calibrating neutral FST from intergenic SNPs...\n")

vcf_intergenic <- file.path(intergenic_dir, "intergenic_ldpruned.vcf")
if (!file.exists(vcf_intergenic)) {
  stop("VCF file not created: ", vcf_intergenic)
}

vcf <- read.vcfR(vcf_intergenic, verbose = FALSE)
geno <- vcfR::extract.gt(vcf)
sampleNames <- colnames(geno)
popNames <- sapply(strsplit(sampleNames, "_"), "[", 1)

cat("Populations:", paste(unique(popNames), collapse = ", "), "\n")
print(table(popNames))

# Recode genotypes for intergenic
G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
G[geno %in% c("0/0", "0|0")] <- 0
G[geno %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
G <- replace(G, is.na(G), 9)

snpNames_intergenic <- vcfR::getID(vcf)

# Calculate FST for intergenic SNPs
cat("Calculating FST for intergenic SNPs...\n")
fst_intergenic <- MakeDiploidFSTMat(t(G), locusNames = snpNames_intergenic, popNames = popNames)
cat("Intergenic FST summary:\n")
print(summary(fst_intergenic$FST))

# Select neutral SNPs (FST < 0.2)
neutral_idx <- which(fst_intergenic$FST >= -0.1 & fst_intergenic$FST < 0.2)
cat("Neutral SNPs for calibration:", length(neutral_idx), "\n")

# Run OutFLANK to estimate neutral distribution
cat("Running OutFLANK on intergenic SNPs...\n")
out_trim <- OutFLANK(
  fst_intergenic[neutral_idx, ],
  NumberOfSamples = length(unique(popNames)),
  qthreshold = 0.05,
  Hmin = 0.1,
  RightTrimFraction = 0.05,
  LeftTrimFraction = 0.05
)

cat("Inferred df:", out_trim$dfInferred, "\n")
cat("Mean FST (no corr):", out_trim$FSTNoCorrbar, "\n")

# Save neutral SNPs list
neutral_snps <- fst_intergenic$LocusName[neutral_idx]
writeLines(neutral_snps, file.path(output_dir, "outflank", "neutral_snps.txt"))

# --- STEP 2b: Apply neutral parameters to ALL QC'd SNPs ---
cat("\nStep 2b: Applying neutral parameters to ALL SNPs...\n")

# Create VCF with ALL SNPs
all_snps_vcf <- file.path(intergenic_dir, "all_snps.vcf")
system(paste(
  "plink2",
  "--bfile", bed_file,
  "--export vcf",
  "--out", file.path(intergenic_dir, "all_snps"),
  "--silent"
))

# Read ALL SNPs VCF
cat("Reading all SNPs VCF...\n")
vcf_all <- read.vcfR(all_snps_vcf, verbose = FALSE)
geno_all <- vcfR::extract.gt(vcf_all)
sampleNames_all <- colnames(geno_all)
popNames_all <- sapply(strsplit(sampleNames_all, "_"), "[", 1)

# Recode genotypes for ALL SNPs
G_all <- matrix(NA, nrow = nrow(geno_all), ncol = ncol(geno_all))
G_all[geno_all %in% c("0/0", "0|0")] <- 0
G_all[geno_all %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G_all[geno_all %in% c("1/1", "1|1")] <- 2
G_all <- replace(G_all, is.na(G_all), 9)

snpNames_all <- vcfR::getID(vcf_all)
posNames_all <- vcfR::getPOS(vcf_all)
locusNames_all <- paste(snpNames_all, posNames_all, sep = ".")

cat("Total SNPs for outlier detection:", length(snpNames_all), "\n")

# Calculate FST for ALL SNPs
cat("Calculating FST for all SNPs...\n")
fst_all <- MakeDiploidFSTMat(t(G_all), locusNames = locusNames_all, popNames = popNames_all)
cat("All SNPs FST summary:\n")
print(summary(fst_all$FST))

# Apply neutral parameters to ALL SNPs
cat("Finding outliers among all SNPs...\n")
P1 <- pOutlierFinderChiSqNoCorr(
  fst_all,
  Fstbar = out_trim$FSTNoCorrbar,
  dfInferred = out_trim$dfInferred,
  qthreshold = 0.05,
  Hmin = 0.1
)

# Extract SNP ID from LocusName (format: SNP_ID.position)
P1$SNP_ID <- sapply(strsplit(P1$LocusName, "\\."), "[", 1)

# Get outliers
outflank_outliers <- P1 |>
  filter(OutlierFlag == TRUE) |>
  arrange(qvalues)

cat("OutFLANK outliers:", nrow(outflank_outliers), "\n")

# Save results
write_csv(P1, file.path(output_dir, "outflank", "outflank_all_results.csv"))
write_csv(outflank_outliers, file.path(output_dir, "outflank", "outflank_outliers.csv"))
if (nrow(outflank_outliers) > 0) {
  writeLines(outflank_outliers$SNP_ID,
             file.path(output_dir, "outflank", "outlier_snps.txt"))
}

# FST distribution plot
png(file.path(output_dir, "outflank", "fst_distribution.png"), width = 800, height = 600)
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE, NoCorr = TRUE,
                       Hmin = 0.1, binwidth = 0.005)
dev.off()

cat("OutFLANK analysis complete.\n\n")

# =============================================================================
# STEP 3: pcadapt Analysis (uses ALL LD-pruned SNPs, not just intergenic)
# =============================================================================
cat("=== STEP 3: pcadapt Analysis ===\n")

# The original pcadapt used LD-pruned SNPs from QC step
# Use the ORIGINAL prune.in file if it exists, otherwise create new one
# Try both container path and relative path
cat("Creating LD-pruned dataset from all QC'd SNPs for pcadapt...\n")
system(paste(
  "plink2",
  "--bfile", bed_file,
  "--indep-pairwise 5 1 0.1",
  "--out", file.path(intergenic_dir, "all_indepSNP"),
  "--silent"
))

system(paste(
  "plink2",
  "--bfile", bed_file,
  "--extract", file.path(intergenic_dir, "all_indepSNP.prune.in"),
  "--make-bed",
  "--out", file.path(intergenic_dir, "all_ldpruned"),
  "--silent"
))

all_ld_snps <- readLines(file.path(intergenic_dir, "all_indepSNP.prune.in"))
cat("All LD-pruned SNPs for pcadapt:", length(all_ld_snps), "\n")

# Use ALL LD-pruned dataset for pcadapt (same as original)
bed_path <- paste0(file.path(intergenic_dir, "all_ldpruned"), ".bed")
x <- read.pcadapt(bed_path, type = "bed")
n_ind <- attr(x, "n")
n_snps <- attr(x, "p")
cat("Loaded", n_ind, "individuals and", n_snps, "SNPs\n")

# Scree plot
pcadapt_scree <- pcadapt(x, K = 10)
png(file.path(output_dir, "pcadapt", "screeplot.png"), width = 800, height = 600)
plot(pcadapt_scree, option = "screeplot")
dev.off()

# Run pcadapt with K=2, min.maf=0.1 (original parameters)
cat("Running pcadapt with K=2, min.maf=0.1...\n")
pcadapt_res <- pcadapt(
  x,
  K = 2,
  method = "mahalanobis",
  min.maf = 0.1,  # Original used 0.1
  LD.clumping = NULL,
  tol = 1e-04
)

# Read BIM for SNP info (use all_ldpruned for pcadapt)
bim_file <- paste0(file.path(intergenic_dir, "all_ldpruned"), ".bim")
snp_info <- read_delim(
  bim_file,
  col_names = c("chr", "snp_id", "cm", "pos", "a1", "a2"),
  show_col_types = FALSE
)

# Combine results
pcadapt_results <- snp_info |>
  mutate(
    pvalue = pcadapt_res$pvalues,
    stat = pcadapt_res$stat,
    maf = pcadapt_res$maf
  ) |>
  drop_na(pvalue)

# Use Benjamini-Hochberg correction (same as original)
pcadapt_results <- pcadapt_results |>
  mutate(
    padj = p.adjust(pvalue, method = "BH"),
    log10p = -log10(pvalue)
  )

# Get outliers at alpha = 0.05
pcadapt_outliers <- pcadapt_results |>
  filter(padj < 0.05) |>
  arrange(pvalue)

cat("pcadapt outliers (BH adj p < 0.05):", nrow(pcadapt_outliers), "\n")

# Save results
write_csv(pcadapt_results, file.path(output_dir, "pcadapt", "pcadapt_all_results.csv"))
write_csv(pcadapt_outliers, file.path(output_dir, "pcadapt", "pcadapt_outliers.csv"))
writeLines(pcadapt_outliers$snp_id, file.path(output_dir, "pcadapt", "outlier_snps.txt"))

# Q-Q plot
png(file.path(output_dir, "pcadapt", "qq_plot.png"), width = 800, height = 600)
plot(pcadapt_res, option = "qqplot")
dev.off()

# Manhattan plot
png(file.path(output_dir, "pcadapt", "manhattan_plot.png"), width = 1200, height = 600)
pcadapt_results |>
  mutate(chr_factor = as.factor(chr)) |>
  ggplot(aes(x = pos, y = log10p, color = chr_factor)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept = -log10(max(pcadapt_outliers$pvalue)),
             linetype = "dashed", color = "red") +
  scale_color_manual(values = rep(c("#1f77b4", "#ff7f0e"), 100)) +
  labs(
    title = "pcadapt Selection Scan (Original Parameters)",
    subtitle = paste("K=2, MAF>0.1 |", nrow(pcadapt_outliers), "outliers (BH adj p < 0.05)"),
    x = "Position",
    y = "-log10(p-value)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

cat("pcadapt analysis complete.\n\n")

# =============================================================================
# =============================================================================
# Summary
# =============================================================================
cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("Output directory:", output_dir, "\n")
cat("\nResults:\n")
cat("  OutFLANK outliers:", nrow(outflank_outliers), "\n")
cat("  pcadapt outliers:", nrow(pcadapt_outliers), "\n")
cat("============================================================\n")

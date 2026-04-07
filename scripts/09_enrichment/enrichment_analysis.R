#!/usr/bin/env Rscript
# ============================================================================
# Functional Enrichment Analysis for Selection Scan Outliers
# ============================================================================
#
# Analyzes the 157 SNPs from the UNION of autogeny selection scan comparisons:
# - NON-AUTO_AUTO (pairwise): Lab vs Autogenous
# - NON-AUTO-FIELD_AUTO (pairwise): Field vs Autogenous
# - NON-AUTO-FIELD_NON-AUTO_AUTO (3-way): All populations
#
# The 17 SNPs in the intersection represent highest-confidence candidates.
# The 157 union SNPs are all potential autogeny-related outliers.
#
# Population labels (per Peter's request):
# - AUTO: Autogenous population (selected for egg production without blood meal)
# - NON-AUTO: Non-autogenous lab population
# - NON-AUTO-FIELD: Non-autogenous field population
#
# This script replicates the original analysis from:
# - markdown/04.pcadapt.Rmd (Venn diagrams)
# - markdown/06.SNPs_functional_annotation.Rmd (SnpEff + bedtools)
#
# Author: Luciano Cosme
# Date: December 2025
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) > 0) args[1] else "."

cat("=============================================================\n")
cat("Functional Enrichment Analysis - Selection Scan Outliers\n")
cat("=============================================================\n\n")

# Create output directory
output_dir <- file.path(project_dir, "output", "enrichment")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load SNP Sets from Venn Diagram Analysis
# -----------------------------------------------------------------------------
cat("1. Loading SNP sets from selection scan analysis...\n")

# 157 SNPs = UNION of all autogeny comparisons (THE KEY SET)
union_file <- file.path(project_dir, "data", "validation", "pcadapt", "outlier_157_SNPs.txt")
# Also check SNPs_158.txt
if (!file.exists(union_file)) {
  union_file <- file.path(project_dir, "data", "validation", "snpeff", "SNPs_158.txt")
}

# 17 SNPs = INTERSECTION (highest confidence)
intersect_file <- file.path(project_dir, "data", "validation", "pcadapt",
                            "4_way_venn_common_SNPs_pcadapt_outflank.txt")

# Load SNP lists
if (file.exists(union_file)) {
  union_snps <- read.table(union_file, header = FALSE, stringsAsFactors = FALSE)$V1
  cat("   Loaded", length(union_snps), "SNPs from UNION (all autogeny comparisons)\n")
} else {
  stop("ERROR: Union SNP file not found: ", union_file)
}

if (file.exists(intersect_file)) {
  intersect_snps <- read.table(intersect_file, header = FALSE, stringsAsFactors = FALSE)$V1
  cat("   Loaded", length(intersect_snps), "SNPs from INTERSECTION (highest confidence)\n")
}

# Load pairwise comparison files for context
# Note: File names use old labels, display names use Peter's labels
pcadapt_dir <- file.path(project_dir, "data", "validation", "pcadapt")
pairwise_files <- list(
  "NON-AUTO_AUTO" = file.path(pcadapt_dir, "man_aut_common_SNPs_pcadapt_outflank.txt"),
  "NON-AUTO-FIELD_AUTO" = file.path(pcadapt_dir, "new_aut_common_SNPs_pcadapt_outflank.txt"),
  "3-way" = file.path(pcadapt_dir, "common_SNPs_pcadapt_outflank.txt")
)

pairwise_counts <- sapply(pairwise_files, function(f) {
  if (file.exists(f)) length(readLines(f)) else 0
})
cat("\n   Pairwise comparison counts:\n")
for (name in names(pairwise_counts)) {
  cat(sprintf("      %s: %d SNPs\n", name, pairwise_counts[name]))
}

# -----------------------------------------------------------------------------
# 2. Load SnpEff Annotations for 157 SNPs
# -----------------------------------------------------------------------------
cat("\n2. Loading SnpEff annotations...\n")

snpeff_file <- file.path(project_dir, "data", "validation", "snpeff", "snps_autogenous.txt")
if (file.exists(snpeff_file)) {
  # Read SnpEff output (space-delimited, 9 columns)
  snpeff <- read.table(snpeff_file, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
  colnames(snpeff)[1:min(9, ncol(snpeff))] <- c("Scaffold", "Position", "SNP", "Effect", "Impact",
                               "Gene", "Gene_Name", "HGVS_c", "HGVS_p")[1:min(9, ncol(snpeff))]
  cat("   Loaded SnpEff annotations for", nrow(snpeff), "SNPs\n")

  # Impact summary
  cat("\n   Impact breakdown:\n")
  impact_table <- table(snpeff$Impact)
  for (imp in sort(names(impact_table), decreasing = TRUE)) {
    cat(sprintf("      %s: %d SNPs\n", imp, impact_table[imp]))
  }

  # Save clean annotation table (handle missing columns)
  available_cols <- intersect(c("Scaffold", "Position", "SNP", "Effect", "Impact", "Gene", "HGVS_c", "HGVS_p"),
                              names(snpeff))
  snpeff_clean <- snpeff |>
    select(all_of(available_cols)) |>
    filter(SNP %in% union_snps)
  write.csv(snpeff_clean, file.path(output_dir, "snpeff_157_snps.csv"), row.names = FALSE)

} else {
  cat("   WARNING: SnpEff file not found\n")
}

# -----------------------------------------------------------------------------
# 3. Load Genes Within 10kb of SNPs
# -----------------------------------------------------------------------------
cat("\n3. Loading genes within 10kb of SNPs...\n")

genes_10kb_file <- file.path(project_dir, "data", "validation", "snpeff", "genes_10kb_snps.txt")
if (file.exists(genes_10kb_file)) {
  genes_10kb <- read.table(genes_10kb_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(genes_10kb) <- c("SNP", "Gene_ID")

  # Filter to our 157 SNPs
  genes_10kb_filtered <- genes_10kb |> filter(SNP %in% union_snps)

  cat("   SNPs with nearby genes:", n_distinct(genes_10kb_filtered$SNP), "\n")
  cat("   Unique genes within 10kb:", n_distinct(genes_10kb_filtered$Gene_ID), "\n")

  # Save gene list
  unique_genes <- unique(genes_10kb_filtered$Gene_ID)
  writeLines(unique_genes, file.path(output_dir, "genes_10kb_157_snps.txt"))
  write.csv(genes_10kb_filtered, file.path(output_dir, "snp_gene_mapping_10kb.csv"), row.names = FALSE)

} else {
  cat("   WARNING: Genes 10kb file not found\n")
  unique_genes <- character(0)
}

# -----------------------------------------------------------------------------
# 4. Cross-reference with DE Genes (RNAseq)
# -----------------------------------------------------------------------------
cat("\n4. Cross-referencing with differentially expressed genes...\n")

de_file <- file.path(project_dir, "data", "files", "MANvsAUTO_sig_mRNAs.csv")
snps_on_de_file <- file.path(project_dir, "data", "validation", "snpeff", "SNPs_79_DE.txt")

if (file.exists(de_file)) {
  de_genes <- read.csv(de_file)
  cat("   Differentially expressed genes:", nrow(de_genes), "\n")

  # Find overlap between genes near SNPs and DE genes
  if (length(unique_genes) > 0) {
    overlap_genes <- intersect(unique_genes, de_genes$gene)
    cat("   Genes near SNPs AND differentially expressed:", length(overlap_genes), "\n")

    if (length(overlap_genes) > 0) {
      overlap_df <- de_genes |>
        filter(gene %in% overlap_genes) |>
        select(gene, gene.name, log2FoldChange, padj)
      write.csv(overlap_df, file.path(output_dir, "de_genes_near_snps.csv"), row.names = FALSE)

      cat("\n   DE genes near selection scan outliers:\n")
      for (i in seq_len(min(10, nrow(overlap_df)))) {
        row <- overlap_df[i, ]
        cat(sprintf("      %s (%s): log2FC = %.2f\n",
                    row$gene, ifelse(is.na(row$gene.name), "NA", row$gene.name), row$log2FoldChange))
      }
    }
  }
}

if (file.exists(snps_on_de_file)) {
  snps_on_de <- readLines(snps_on_de_file)
  cat("\n   SNPs directly on DE genes (from original analysis):", length(snps_on_de), "\n")
}

# -----------------------------------------------------------------------------
# 5. Key Variants Summary
# -----------------------------------------------------------------------------
cat("\n5. Key variants summary...\n")

if (exists("snpeff")) {
  # MODERATE impact (missense)
  moderate <- snpeff |> filter(Impact == "MODERATE")
  if (nrow(moderate) > 0) {
    cat("\n   MODERATE impact variants (missense):\n")
    for (i in seq_len(nrow(moderate))) {
      row <- moderate[i, ]
      cat(sprintf("      %s (%s:%s) - %s %s\n",
                  row$SNP, row$Scaffold, row$Position, row$Effect, row$HGVS_p))
    }
  }

  # LOW impact (synonymous)
  low <- snpeff |> filter(Impact == "LOW")
  cat(sprintf("\n   LOW impact variants (synonymous): %d\n", nrow(low)))

  # MODIFIER (non-coding)
  modifier <- snpeff |> filter(Impact == "MODIFIER")
  cat(sprintf("   MODIFIER impact (non-coding): %d\n", nrow(modifier)))
}

# -----------------------------------------------------------------------------
# 6. Intersection vs Union Analysis
# -----------------------------------------------------------------------------
cat("\n6. Intersection vs Union analysis...\n")

cat("   17 SNPs in intersection (all autogeny comparisons agree):\n")
cat("      - These are HIGHEST CONFIDENCE selection candidates\n")
cat("      - Present in MAN_AUT AND NEW_AUT AND 3-way\n")

cat("\n   157 SNPs in union (any autogeny comparison):\n")
cat("      - ALL potential selection candidates\n")
cat("      - Different comparisons may capture different signals\n")

# Check which comparisons each intersection SNP appears in
if (length(intersect_snps) > 0 && exists("snpeff")) {
  intersect_annot <- snpeff |> filter(SNP %in% intersect_snps)
  write.csv(intersect_annot, file.path(output_dir, "intersection_17_snps.csv"), row.names = FALSE)
  cat(sprintf("\n   Intersection SNPs with annotations: %d\n", nrow(intersect_annot)))
}

# -----------------------------------------------------------------------------
# 7. Generate Summary Report
# -----------------------------------------------------------------------------
cat("\n7. Generating summary report...\n")

report <- c(
  "=============================================================",
  "Functional Enrichment Analysis Report",
  "Selection Scan Outliers - Autogeny GWAS",
  "=============================================================",
  "",
  "POPULATIONS (per Peter's labels):",
  "  - AUTO: Autogenous (selected for egg production without blood meal)",
  "  - NON-AUTO: Non-autogenous lab population",
  "  - NON-AUTO-FIELD: Non-autogenous field population",
  "",
  "VENN DIAGRAM ANALYSIS",
  "---------------------",
  "Comparisons involving AUTO population:",
  sprintf("  - NON-AUTO vs AUTO (pairwise): %d SNPs", pairwise_counts["NON-AUTO_AUTO"]),
  sprintf("  - NON-AUTO-FIELD vs AUTO (pairwise): %d SNPs", pairwise_counts["NON-AUTO-FIELD_AUTO"]),
  sprintf("  - 3-way comparison: %d SNPs", pairwise_counts["3-way"]),
  "",
  "KEY SNP SETS",
  "------------",
  sprintf("UNION (all autogeny outliers): %d SNPs", length(union_snps)),
  sprintf("INTERSECTION (all comparisons agree): %d SNPs", length(intersect_snps)),
  "",
  "The UNION represents ALL potential selection candidates.",
  "The INTERSECTION represents HIGHEST CONFIDENCE candidates.",
  ""
)

if (exists("snpeff")) {
  report <- c(report,
    "SNPEFF ANNOTATION SUMMARY",
    "-------------------------")
  for (imp in sort(names(impact_table), decreasing = TRUE)) {
    report <- c(report, sprintf("  %s: %d SNPs", imp, impact_table[imp]))
  }
}

if (length(unique_genes) > 0) {
  report <- c(report, "",
    "GENES NEAR SELECTION OUTLIERS",
    "-----------------------------",
    sprintf("SNPs with genes within 10kb: %d", n_distinct(genes_10kb_filtered$SNP)),
    sprintf("Unique genes identified: %d", length(unique_genes)))
}

report <- c(report, "",
  "OUTPUT FILES",
  "------------",
  "- snpeff_157_snps.csv: SnpEff annotations for union SNPs",
  "- genes_10kb_157_snps.txt: Gene IDs within 10kb of SNPs",
  "- snp_gene_mapping_10kb.csv: SNP-to-gene mapping",
  "- de_genes_near_snps.csv: DE genes near selection outliers",
  "- intersection_17_snps.csv: Highest-confidence candidates",
  "- enrichment_report.txt: This report",
  "",
  "EXTERNAL ENRICHMENT TOOLS",
  "-------------------------",
  "Use the gene list (genes_10kb_157_snps.txt) with:",
  "- VectorBase: https://vectorbase.org/vectorbase/app/",
  "- g:Profiler: https://biit.cs.ut.ee/gprofiler/",
  "- DAVID: https://david.ncifcrf.gov/"
)

writeLines(report, file.path(output_dir, "enrichment_report.txt"))

cat(paste(report, collapse = "\n"))
cat("\n\nAnalysis complete!\n")
cat("Output files saved to:", output_dir, "\n")

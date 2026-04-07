#!/usr/bin/env Rscript
# ==============================================================================
# 02_venn_diagrams.R
# Figure 1: 4-way Venn diagram (outlier SNPs)
#
# EXTRACTED FROM: dropbox_original/autogenous/scripts/markdown_files/04.pcadapt.Rmd
# Section: 6. Venn diagram (lines 1740-1802)
# ==============================================================================

# Load libraries
suppressPackageStartupMessages({
  library(here)
  library(ggvenn)
  library(tidyr)
  library(ggplot2)
})

# ==============================================================================
# Configuration
# ==============================================================================

# Input: Original data from dropbox
input_base <- here("dropbox_original", "autogenous", "output", "pcadapt")

# Output: Comparison directory
output_dir <- here("output", "figures_comparison", "replicated")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# Load Data (EXACTLY as in original Rmd)
# ==============================================================================

cat("Loading SNP data...\n")

NEW_MAN_AUT <- read.table(
  file.path(input_base, "common_SNPs_pcadapt_outflank.txt"),
  stringsAsFactors = FALSE
) |> drop_na()

MAN_AUT <- read.table(
  file.path(input_base, "man_aut_common_SNPs_pcadapt_outflank.txt"),
  stringsAsFactors = FALSE
)

NEW_AUT <- read.table(
  file.path(input_base, "new_aut_common_SNPs_pcadapt_outflank.txt"),
  stringsAsFactors = FALSE
)

NEW_MAN <- read.table(
  file.path(input_base, "man_new_SNPs_pcadapt.txt"),
  stringsAsFactors = FALSE
)

cat("NEW_MAN_AUT:", nrow(NEW_MAN_AUT), "SNPs\n")
cat("MAN_AUT:", nrow(MAN_AUT), "SNPs\n")
cat("NEW_AUT:", nrow(NEW_AUT), "SNPs\n")
cat("NEW_MAN:", nrow(NEW_MAN), "SNPs\n")

# ==============================================================================
# Create Venn Diagram (EXACTLY as in original Rmd)
# ==============================================================================

# Create a list with all dataframes
list_of_clusters <- list(
  "NEW_MAN_AUT" = NEW_MAN_AUT$V1,
  "MAN_AUT" = MAN_AUT$V1,
  "NEW_AUT" = NEW_AUT$V1,
  "NEW_MAN" = NEW_MAN$V1
)

# Create Venn diagram
venn_diagram <- ggvenn(list_of_clusters)

# Increase plot margins
venn_diagram <- venn_diagram + theme(plot.margin = margin(60, 60, 60, 60, "points"))

# Resize text
venn_diagram <- venn_diagram + theme(text = element_text(size = 5))

# Print the adjusted plot
print(venn_diagram)

# ==============================================================================
# Save 4-way Venn diagram
# ==============================================================================

output_path <- file.path(output_dir, "4_way_venn_significant_snps.pdf")
ggsave(output_path, venn_diagram, height = 8, width = 8, dpi = 300)
cat("\nSaved:", output_path, "\n")

# ==============================================================================
# 2-way Venn diagrams (PCAdapt vs Outflank for each comparison)
# EXTRACTED FROM: 04.pcadapt.Rmd - individual comparison sections
# ==============================================================================

cat("\nCreating 2-way Venn diagrams...\n")

# Load PCAdapt and Outflank results for each comparison
pcadapt_base <- here("dropbox_original", "autogenous", "output", "pcadapt")
outflank_base <- here("dropbox_original", "autogenous", "output", "outflank")

# --- NEW vs AUT comparison ---
new_aut_pcadapt <- read.table(
  file.path(pcadapt_base, "new_aut_SNPs_pcadapt.txt"),
  stringsAsFactors = FALSE
)
new_aut_outflank <- read.table(
  file.path(outflank_base, "new_aut_SNPs_outflank.txt"),
  stringsAsFactors = FALSE
)

venn_new_aut <- ggvenn(list(
  "PCAdapt" = new_aut_pcadapt$V1,
  "Outflank" = new_aut_outflank$V1
)) + theme(plot.margin = margin(30, 30, 30, 30, "points"))

ggsave(file.path(output_dir, "venn_NEW_AUT_pcadapt_outflank.pdf"),
       venn_new_aut, height = 6, width = 6, dpi = 300)
ggsave(file.path(output_dir, "venn_NEW_AUT_pcadapt_outflank.png"),
       venn_new_aut, height = 6, width = 6, dpi = 300)
cat("Saved: venn_NEW_AUT_pcadapt_outflank.pdf/png\n")

# --- MAN vs AUT comparison ---
man_aut_pcadapt <- read.table(
  file.path(pcadapt_base, "man_aut_SNPs_pcadapt.txt"),
  stringsAsFactors = FALSE
)
man_aut_outflank <- read.table(
  file.path(outflank_base, "man_aut_SNPs_outflank.txt"),
  stringsAsFactors = FALSE
)

venn_man_aut <- ggvenn(list(
  "PCAdapt" = man_aut_pcadapt$V1,
  "Outflank" = man_aut_outflank$V1
)) + theme(plot.margin = margin(30, 30, 30, 30, "points"))

ggsave(file.path(output_dir, "venn_MAN_AUT_pcadapt_outflank.pdf"),
       venn_man_aut, height = 6, width = 6, dpi = 300)
ggsave(file.path(output_dir, "venn_MAN_AUT_pcadapt_outflank.png"),
       venn_man_aut, height = 6, width = 6, dpi = 300)
cat("Saved: venn_MAN_AUT_pcadapt_outflank.pdf/png\n")

# --- 3-pop comparison (NEW_MAN_AUT) ---
three_pop_pcadapt <- read.table(
  file.path(pcadapt_base, "SNPs_pcadapt.txt"),
  stringsAsFactors = FALSE
)
three_pop_outflank <- read.table(
  file.path(outflank_base, "SNPs_outflank.txt"),
  stringsAsFactors = FALSE
)

venn_3_pop <- ggvenn(list(
  "PCAdapt" = three_pop_pcadapt$V1,
  "Outflank" = three_pop_outflank$V1
)) + theme(plot.margin = margin(30, 30, 30, 30, "points"))

ggsave(file.path(output_dir, "venn_3_pop_pcadapt_outflank.pdf"),
       venn_3_pop, height = 6, width = 6, dpi = 300)
ggsave(file.path(output_dir, "venn_3_pop_pcadapt_outflank.png"),
       venn_3_pop, height = 6, width = 6, dpi = 300)
cat("Saved: venn_3_pop_pcadapt_outflank.pdf/png\n")

# Find the intersect of NEW_MAN_AUT, MAN_AUT, and NEW_AUT
common_elements <- Reduce(intersect, list(NEW_MAN_AUT$V1, MAN_AUT$V1, NEW_AUT$V1))
cat("\nCommon elements across 3 groups:", length(common_elements), "\n")

cat("\nDone! Created 4-way and 3x 2-way Venn diagrams.\n")

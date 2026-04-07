#!/usr/bin/env Rscript
# ============================================================
# Allele and Genotype Frequency Plots
# Replicates: File_S7.Frequencies.Rmd Sections 2-4
# ============================================================
# Prerequisites:
#   output/quality_control/file7.bed/bim/fam  (from S1)
#   output/snpeff/SNPs_158.txt                 (from S6)
#   output/pcadapt/4_way_venn_common_SNPs_pcadapt_outflank.txt (17 SNPs, from S3)
#   output/ldna/pop/chr2/clusters/             (cluster 14 SNPs, from S4e)
#
# Output:
#   output/frequencies/figures/significant_17_snps_alleles.pdf
#   output/frequencies/figures/significant_17_snps_genotypes.pdf
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(viridis)
  library(here)
})

args <- commandArgs(trailingOnly = TRUE)
output_dir <- ifelse(length(args) >= 1, args[1], "output/frequencies")

cat("============================================================\n")
cat("Allele and Genotype Frequency Plots\n")
cat("============================================================\n")
cat("Output:", output_dir, "\n")
cat("============================================================\n\n")

dir.create(file.path(output_dir, "snp_files"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)

# ============================================================
# STEP 1: Compute allele frequencies per population
# ============================================================
cat("=== Step 1: Computing allele frequencies ===\n")

pops <- unique(read.table("output/quality_control/file7.fam")$V1)
cat("Populations:", paste(pops, collapse = ", "), "\n")

for (pop in pops) {
  system(paste(
    "echo", pop, "|",
    "plink2",
    "--extract output/snpeff/SNPs_158.txt",
    "--bfile output/quality_control/file7",
    "--keep-fam /dev/stdin",
    "--out", file.path(output_dir, pop),
    "--freq",
    "--silent"
  ))
}

# ============================================================
# STEP 2: Build per-SNP frequency files
# ============================================================
cat("=== Step 2: Building per-SNP files ===\n")

afreq_files <- list.files(output_dir, pattern = "\\.afreq$", full.names = TRUE)
unique_snps <- unique(unlist(lapply(afreq_files, function(f) {
  d <- read.table(f, header = TRUE)
  d$ID
})))

for (snp in unique_snps) {
  snp_rows <- list()
  allele_header <- NULL
  for (f in afreq_files) {
    pop_name <- sub("\\.afreq$", "", basename(f))
    d <- read.table(f, header = TRUE)
    row <- d[d$ID == snp, ]
    if (nrow(row) == 1) {
      if (is.null(allele_header)) allele_header <- c(row$ALT, row$REF)
      snp_rows[[pop_name]] <- data.frame(
        SNP = snp, Stratum = pop_name,
        Allele1_Value = 1 - row$ALT_FREQS,
        Allele2_Value = row$ALT_FREQS
      )
    }
  }
  if (length(snp_rows) > 0 && !is.null(allele_header)) {
    df <- do.call(rbind, snp_rows)
    colnames(df)[3:4] <- allele_header
    write.table(df, file.path(output_dir, "snp_files", paste0(snp, ".txt")),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

# ============================================================
# STEP 3: Load all SNP data
# ============================================================
cat("=== Step 3: Loading SNP frequency data ===\n")

files <- list.files(file.path(output_dir, "snp_files"), pattern = "\\.txt$", full.names = TRUE)

all_data <- lapply(files, function(file) {
  dat <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                    colClasses = "character")
  allele_names <- colnames(dat)[3:4]
  colnames(dat)[3:4] <- c("Allele1_Value", "Allele2_Value")
  dat$Allele1_Name <- allele_names[1]
  dat$Allele2_Name <- allele_names[2]
  dat$SNP <- gsub(".txt", "", basename(file))
  dat
})

all_data_df <- do.call(rbind, all_data)

# Remap to manuscript labels: AUT=AUTO, MAN=NON-AUTO, NEW=NON-AUTO-FIELD
pop_labels <- c("AUT" = "AUTO", "MAN" = "NON-AUTO", "NEW" = "NON-AUTO-FIELD")
all_data_df$Stratum <- factor(
  pop_labels[all_data_df$Stratum],
  levels = c("AUTO", "NON-AUTO", "NON-AUTO-FIELD")
)

# ============================================================
# STEP 4: Plot allele frequencies for 17 SNPs
# ============================================================
cat("=== Step 4: Plotting allele frequencies (17 SNPs) ===\n")

venn_file <- "output/pcadapt/4_way_venn_common_SNPs_pcadapt_outflank.txt"
if (!file.exists(venn_file)) {
  cat("WARNING: 17-SNP list not found at", venn_file, "— skipping plot\n")
} else {
  snps_17 <- read.table(venn_file, stringsAsFactors = FALSE)
  snps_aut17 <- all_data_df |> filter(SNP %in% snps_17$V1)

  allele_colors <- c("A" = "#ffccd1", "T" = "#e3fcc2", "C" = "#ffec8f", "G" = "#bdbdfa")
  stratum_border_colors <- c("AUTO" = "red", "NON-AUTO" = "black", "NON-AUTO-FIELD" = "black")

  all_data_long <- snps_aut17 |>
    pivot_longer(cols = c(Allele1_Value, Allele2_Value),
                 names_to = "Allele", values_to = "Value",
                 names_pattern = "Allele(\\d)_Value") |>
    mutate(
      Allele = if_else(Allele == "1", Allele1_Name, Allele2_Name),
      Value  = as.numeric(Value)
    )

  plot_alleles <- ggplot(all_data_long, aes(x = Stratum, y = Value, fill = Allele, group = Allele)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9),
             aes(color = Stratum), linewidth = 0.5) +
    scale_fill_manual(name = "", values = allele_colors) +
    scale_color_manual(name = "", values = stratum_border_colors) +
    facet_wrap(~ SNP, scales = "free_x", ncol = 4) +
    geom_text(aes(label = Allele, y = Value + 0.02, group = Allele),
              position = position_dodge(width = 0.6), vjust = -0.25, size = 3,
              check_overlap = TRUE) +
    scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, by = 0.25)) +
    labs(x = "Population", y = "Frequency") +
    theme_bw() +
    theme(
      strip.text.x     = element_text(size = 10, face = "bold"),
      legend.position  = "top",
      strip.background = element_rect(fill = "#e8e8e8", colour = NA),
      panel.spacing    = unit(1, "lines")
    ) +
    guides(fill = guide_legend(order = 1), color = "none")

  out <- file.path(output_dir, "figures", "significant_17_snps_alleles.pdf")
  ggsave(out, plot_alleles, height = 8, width = 6, dpi = 300)
  cat("Saved:", out, "\n")
}

# ============================================================
# STEP 5: Compute and plot genotype frequencies for 17 SNPs
# ============================================================
cat("=== Step 5: Computing genotype frequencies ===\n")

for (pop in pops) {
  system(paste(
    "echo", pop, "|",
    "plink2",
    "--extract output/snpeff/SNPs_158.txt",
    "--bfile output/quality_control/file7",
    "--keep-fam /dev/stdin",
    "--out", file.path(output_dir, pop),
    "--geno-counts",
    "--silent"
  ))
}

gcount_files <- list.files(output_dir, pattern = "\\.gcount$", full.names = TRUE)

geno_data <- lapply(gcount_files, function(f) {
  pop_name <- sub("\\.gcount$", "", basename(f))
  d <- read.table(f, header = TRUE)
  total <- d$HOM_ALT_CT + d$HET_ALT_CT + d$TWO_ALT_GENO_CT
  data.frame(
    SNP     = d$ID,
    Stratum = pop_name,
    HOM_REF = d$HOM_ALT_CT / total,
    HET     = d$HET_ALT_CT / total,
    HOM_ALT = d$TWO_ALT_GENO_CT / total,
    REF     = d$REF,
    ALT     = d$ALT
  )
})
geno_df <- do.call(rbind, geno_data)

# Remap to manuscript labels
geno_df$Stratum <- factor(
  pop_labels[geno_df$Stratum],
  levels = c("AUTO", "NON-AUTO", "NON-AUTO-FIELD")
)

if (file.exists(venn_file)) {
  snps_17 <- read.table(venn_file, stringsAsFactors = FALSE)
  geno_17 <- geno_df |>
    filter(SNP %in% snps_17$V1) |>
    pivot_longer(cols = c(HOM_REF, HET, HOM_ALT), names_to = "Genotype", values_to = "Freq") |>
    mutate(Genotype = recode(Genotype,
      HOM_REF = paste0(REF, REF),
      HET     = paste0(REF, ALT),
      HOM_ALT = paste0(ALT, ALT)
    ))

  geno_colors <- c(viridis(8), RColorBrewer::brewer.pal(8, "Set2"))
  stratum_border_colors <- c("AUTO" = "red", "NON-AUTO" = "black", "NON-AUTO-FIELD" = "black")

  plot_geno <- ggplot(geno_17, aes(x = Stratum, y = Freq, fill = Genotype)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9),
             aes(color = Stratum), linewidth = 0.5) +
    scale_color_manual(values = stratum_border_colors) +
    facet_wrap(~ SNP, scales = "free_x", ncol = 4) +
    scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, by = 0.25)) +
    labs(x = "Population", y = "Frequency") +
    theme_bw() +
    theme(
      strip.text.x     = element_text(size = 10, face = "bold"),
      legend.position  = "top",
      strip.background = element_rect(fill = "#e8e8e8", colour = NA),
      panel.spacing    = unit(1, "lines")
    ) +
    guides(fill = guide_legend(order = 1), color = "none")

  out <- file.path(output_dir, "figures", "significant_17_snps_genotypes.pdf")
  ggsave(out, plot_geno, height = 8, width = 6, dpi = 300)
  cat("Saved:", out, "\n")
}

cat("\n============================================================\n")
cat("Frequency plots complete\n")
cat("============================================================\n")

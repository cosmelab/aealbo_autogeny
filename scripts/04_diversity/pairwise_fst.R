#!/usr/bin/env Rscript
# ============================================================
# Pairwise Fst with StAMPP
# Replicates: File_S8.Fst.Rmd Section 1
# ============================================================
# Input:  output/ldna/files/file1 (PLINK bed/bim/fam)
# Output: output/fst/pop_pairwise2.rds
#         output/fst/pop_pairwise2_df.csv
#         output/fst/figures/pairwise_estimates.pdf
# ============================================================

suppressPackageStartupMessages({
  library(adegenet)
  library(StAMPP)
  library(reshape2)
  library(ggplot2)
  library(here)
})

args <- commandArgs(trailingOnly = TRUE)
input_prefix <- ifelse(length(args) >= 1, args[1], "output/ldna/files/file1")
output_dir   <- ifelse(length(args) >= 2, args[2], "output/fst")

cat("============================================================\n")
cat("Pairwise Fst Analysis (StAMPP)\n")
cat("============================================================\n")
cat("Input:", input_prefix, "\n")
cat("Output:", output_dir, "\n")
cat("============================================================\n\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)

# Step 1: Convert PLINK to raw format
cat("=== Step 1: Converting PLINK to raw format ===\n")
system(paste(
  "plink",
  "--keep-allele-order",
  "--bfile", input_prefix,
  "--recodeA",
  "--out", file.path(output_dir, "pop_pairwise"),
  "--silent"
))

# Step 2: Read into genlight object
cat("=== Step 2: Reading genlight object ===\n")
pop_pairwise <- read.PLINK(
  file.path(output_dir, "pop_pairwise.raw"),
  quiet = FALSE,
  chunkSize = 1000,
  parallel = require("parallel"),
  n.cores = 4
)
summary(pop_pairwise)

# Step 3: StAMPP pairwise Fst (10 bootstraps, 95% CI, 8 threads)
cat("=== Step 3: Calculating pairwise Fst (StAMPP) ===\n")
pop_pairwise2 <- stamppConvert(pop_pairwise, type = "genlight")
pop_pairwise2 <- stamppFst(pop_pairwise2, 10, 95, 8)

saveRDS(pop_pairwise2, file.path(output_dir, "pop_pairwise2.rds"))
cat("Saved:", file.path(output_dir, "pop_pairwise2.rds"), "\n")

# Step 4: Export as CSV
cat("=== Step 4: Exporting Fst table ===\n")
pop_pairwise2_df <- data.frame(pop_pairwise2)
write.csv(pop_pairwise2_df, file.path(output_dir, "pop_pairwise2_df.csv"))

# Extract Fst values and build symmetric matrix
fst_columns <- grep("Fsts\\.", names(pop_pairwise2_df), value = TRUE)
fst_df <- pop_pairwise2_df[, fst_columns]
names(fst_df) <- gsub("Fsts\\.", "", names(fst_df))
aa <- as.matrix(fst_df)
aa[upper.tri(aa)] <- t(aa)[upper.tri(t(aa))]
aa[lower.tri(aa)] <- NA

cat("Fst matrix:\n")
print(round(aa, 3))

# Remap to manuscript labels for plot axes
pop_labels <- c("AUT" = "AUTO", "MAN" = "NON-AUTO", "NEW" = "NON-AUTO-FIELD")
rownames(aa) <- pop_labels[rownames(aa)]
colnames(aa) <- pop_labels[colnames(aa)]

# Step 5: Heatmap plot
cat("=== Step 5: Plotting pairwise Fst heatmap ===\n")
pairfst.long <- melt(aa)

pairfst.f <- ggplot(pairfst.long, aes(Var1, Var2)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(
    low = "white", high = "#71b6ff",
    name = "Fst", na.value = "white", limits = c(0, 0.5)
  ) +
  scale_x_discrete(position = "top") +
  theme_bw() +
  geom_text(aes(label = ifelse(is.na(value), "", formatC(value, digits = 2, format = "f"))), size = 6) +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, size = 20),
    axis.title   = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y  = element_text(hjust = 0, size = 20)
  )

output_path <- file.path(output_dir, "figures", "pairwise_estimates.pdf")
ggsave(output_path, pairfst.f, height = 4, width = 4, dpi = 300)
cat("Saved:", output_path, "\n")

cat("\n============================================================\n")
cat("Pairwise Fst complete\n")
cat("============================================================\n")

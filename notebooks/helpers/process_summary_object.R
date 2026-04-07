process_summary_object <-
function(summary_object_name) {
  # Select only the relevant columns
  dt <- get(summary_object_name)[, .(SNP_id, REF_mismatch, ALT_mismatch, Zigo_mismatch)]

  # Reshape data to long format
  dt_long <- reshape2::melt(dt, id.vars = "SNP_id", variable.name = "type", value.name = "count")

  # Convert to data.table if it's not already
  setDT(dt_long)

  # Convert count to numeric if it's not already
  dt_long[, count := as.numeric(count)]

  # Count occurrences per count value
  dt_long <- dt_long[, .(n = .N), by = .(type, count)]

  # Calculate total count of unique SNPs
  total_SNP <- length(unique(dt$SNP_id))

  # Add a new column for the percentage
  dt_long[, perc := n / total_SNP * 100]

  # Define new labels
  new_labels <- c(
    "Reference Allele" = "REF_mismatch",
    "Alternative Allele" = "ALT_mismatch",
    "Zygosity Mismatch" = "Zigo_mismatch"
  )

  # Apply new labels
  dt_long$type <- forcats::fct_recode(dt_long$type, !!!new_labels)

  # Return the processed data.table
  return(dt_long)
}

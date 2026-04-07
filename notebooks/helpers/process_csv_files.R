process_csv_files <-
function(csv_file_ending) {
  # Read the CSV file using fread() function
  csv_file <- paste0("output/wgs_vs_chip/combined_comparison_", csv_file_ending, ".csv")
  data_dt <- data.table::fread(csv_file)

  # Get all column names that end with '_gcomp'
  gcomp_cols <- grep("_gcomp$", names(data_dt), value = TRUE)

  # Convert data.frame to data.table
  setDT(data_dt)

  # Iterate over the '_gcomp' columns and create new '_REF' and '_ALT' columns
  for (col in gcomp_cols) {
    # Split each '_gcomp' column into '_REF' and '_ALT'
    ref_col <- paste0(col, "_REF")
    alt_col <- paste0(col, "_ALT")
    data_dt[, c(ref_col, alt_col) := tstrsplit(get(col), ", ", fixed = TRUE)]

    # Remove unwanted characters from each new column
    data_dt[, (ref_col) := gsub("\\[|\\]|'", "", get(ref_col))]
    data_dt[, (alt_col) := gsub("\\[|\\]|'", "", get(alt_col))]
  }

  # Rename columns to remove '_gcomp'
  new_names <- names(data_dt)
  new_names <- gsub("_gcomp_ALT$", "_ALT", new_names)
  new_names <- gsub("_gcomp_REF$", "_REF", new_names)
  setnames(data_dt, new_names)
  setnames(data_dt, new_names)

  # Return the processed data.table
  return(data_dt)
}

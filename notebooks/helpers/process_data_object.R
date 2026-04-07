process_data_object <-
function(object_name) {
  # Get the data.table object based on the input name
  data_dt <- get(object_name)

  # Create columns for match and mismatch count for columns ending with _REF
  cols_REF <- grep("_REF$", names(data_dt), value = TRUE)
  data_dt[, c("REF_match_count", "REF_mismatch_count") := .(
    rowSums(.SD == "match", na.rm = TRUE),
    rowSums(.SD == "mismatch", na.rm = TRUE)
  ), .SDcols = cols_REF]

  # Create columns for match and mismatch count for columns ending with _ALT
  cols_ALT <- grep("_ALT$", names(data_dt), value = TRUE)
  data_dt[, c("ALT_match_count", "ALT_mismatch_count") := .(
    rowSums(.SD == "match", na.rm = TRUE),
    rowSums(.SD == "mismatch", na.rm = TRUE)
  ), .SDcols = cols_ALT]

  # Create columns for match and mismatch count for columns ending with _zcomp
  cols_Zigo <- grep("_zcomp$", names(data_dt), value = TRUE)
  data_dt[, c("Zigo_match_count", "Zigo_mismatch_count") := .(
    rowSums(.SD == "match", na.rm = TRUE),
    rowSums(.SD == "mismatch", na.rm = TRUE)
  ), .SDcols = cols_Zigo]

  # Summarize the data for each SNP_id
  summary_dt <- data_dt[, .(
    REF_match = sum(REF_match_count, na.rm = TRUE),
    REF_mismatch = sum(REF_mismatch_count, na.rm = TRUE),
    ALT_match = sum(ALT_match_count, na.rm = TRUE),
    ALT_mismatch = sum(ALT_mismatch_count, na.rm = TRUE),
    Zigo_match = sum(Zigo_match_count, na.rm = TRUE),
    Zigo_mismatch = sum(Zigo_mismatch_count, na.rm = TRUE)
  ), by = SNP_id]

  # Sort the summarized data by SNP_id
  setorder(summary_dt, SNP_id)

  # Return the processed summary data.table
  return(summary_dt)
}

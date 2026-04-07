read_fam_file <-
function(file_path) {
  fam_data <- readr::read_delim(
    file = file_path,
    col_types = "ccccccc",
    col_names = FALSE
  )
  return(fam_data)
}

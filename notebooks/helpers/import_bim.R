import_bim <-
function(file) {
  # import as a tibble and set columns as integers
  bim <-
    read_delim(
      file,
      col_names      = FALSE,
      show_col_types = FALSE,
      col_types      = "ccidcc"
    )
  # rename the columns by index
  bim <- bim |>
    rename(
      Scaffold = 1,
      SNP      = 2,
      Cm       = 3,
      Position = 4,
      Allele1  = 5,
      Allele2  = 6
    )
  return(bim)
}

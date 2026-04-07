read_fam_file2 <-
function(file_path) {
  library(readr)
  library(dplyr)

  fam <- read_delim(
    file_path,
    col_names      = FALSE,
    show_col_types = FALSE,
    col_types      = "cnnnnn"
  ) |>
  rename(
    family       = 1,
    id           = 2,
    father       = 3,
    mother       = 4,
    sex          = 5,
    phenotype    = 6
  )

  return(fam)
}

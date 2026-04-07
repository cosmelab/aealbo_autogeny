import_ld_files <-
function(directory) {
  ld_files <-
    list.files(
      path = directory,
      pattern = "\\.ld$",
      full.names = TRUE
    )

  ld_tibbles <- list()

  for (file in ld_files) {
    ld_name <- gsub(".ld", "", basename(file))
    ld_data <- read_delim(
      file,
      col_names = TRUE,
      delim = " ",
      show_col_types = FALSE
    ) |>
      select(c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2"))
    ld_tibbles[[ld_name]] <- ld_data
  }

  return(ld_tibbles)
}

write_chunk <-
function(rmd_file_path, label, output_directory, script_name) {
  # Read the Rmd file
  rmd_lines <- readLines(rmd_file_path)

  # Find the line numbers of the chunk start and end
  start_line <- grep(paste0("```\\{python ", label), rmd_lines)
  end_line <- grep("```", rmd_lines[(start_line + 1):length(rmd_lines)]) + start_line

  # Extract the content of the chunk
  chunk_content <- rmd_lines[(start_line + 1):(end_line - 1)]

  # Save the content of the chunk to a file
  file_path <- file.path(output_directory, script_name)
  cat(chunk_content, sep = "\n", file = file_path)
  message("Python script saved as: ", file_path)
}

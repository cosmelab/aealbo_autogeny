get_half_pop <-
function(ld_tibble) {
  # Split the tibble into sub-tibbles by scaffold
  ld_split <- split(ld_tibble, ld_tibble[, 1])
  # Remove rows with missing data
  ld_no_nan <- lapply(ld_split, function(x) {
    x[complete.cases(x), ]
  })
  # Calculate the distance between SNPs
  ld_dist <-
    lapply(ld_no_nan, function(x) {
      cbind(x, x[, 5] - x[, 2])
    })
  # Sort the sub-tibbles by distance
  ld_dist_order <- lapply(ld_dist, function(x) {
    x[order(x[, 8]), ]
  })
  # Combine the sub-tibbles into one tibble
  ld_proc <- do.call(rbind, ld_dist_order)
  # Sort the tibble by distance
  ld_sorted <- ld_proc[order(ld_proc[, 8]), ]
  # Bin the data by distance
  bin.1000 <-
    cut(ld_sorted[, 8], (c(0:1000) * 1000), include.lowest = TRUE)
  # Calculate the mean R2 for each bin
  ld_bin_1000_sorted <-
    tapply(ld_sorted[, 7], bin.1000, function(x) {
      mean(x)
    })
  # Fit a LOESS curve to the binned data
  lo.1000.ld <- loess(ld_bin_1000_sorted ~ c(1:1000))
  # Find the value of the binned R2 that corresponds to half the maximum R2 value
  d.pop <- max(ld_bin_1000_sorted, na.rm = TRUE)
  half.pop <- which(ld_bin_1000_sorted <= (d.pop / 2))[1]
  return(half.pop)
}

# This file contains a few helper functions. We are at the stage where we do not need a full package and a simple file being sourced is sufficient for now.

total_momentum_str_to_vec <- function (total_momentum_str) {
  parts <- stringr::str_match(total_momentum_str, '(-?\\d)(-?\\d)(-?\\d)')
  sapply(parts[1, 2:4], as.integer)
}

get_workdir <- function () {
  workdir_candidates <- c(
    '~/Lattice/NG2',
    '~/NG2'
  )
  
  for (path in workdir_candidates) {
    if (dir.exists(path)) {
      return (path)
    }
  }
  
  return (NA)
}
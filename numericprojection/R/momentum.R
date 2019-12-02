momentum_str_to_vec <- function (total_momentum_str) {
  parts <- stringr::str_match(total_momentum_str, '(-?\\d)(-?\\d)(-?\\d)')
  sapply(parts[1, 2:4], as.integer, USE.NAMES = FALSE)
}

momentum_vec_to_str <- function (momentum) {
  paste0(sprintf("%d", momentum), collapse = '')
}

momentum_vec_to_sq <- function (momentum) {
  sum(momentum^2)
}

total_momentum_ref_vec <- function (total_momentum_sq) {
  if (total_momentum_sq == 0)
    c(0, 0, 0)
  else if (total_momentum_sq == 1)
    c(0, 0, 1)
  else if (total_momentum_sq == 2)
    c(1, 1, 0)
  else if (total_momentum_sq == 3)
    c(1, 1, 1)
  else if (total_momentum_sq == 4)
    c(0, 0, 2)
}

sort_momenta_vec_list <- function (momenta_vec) {
  momenta_str <- sapply(momenta_vec, momentum_vec_to_str, USE.NAMES = FALSE)
  lapply(sort(momenta_str), momentum_str_to_vec)
}

# This file contains a few helper functions. We are at the stage where we do not need a full package and a simple file being sourced is sufficient for now.

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

drop_empty <- function (l) {
  lengths <- sapply(l, length)
  l[lengths > 0]
}

lapplyn <- function (container, f, n, ...) {
  if (n == 1) {
    result <- lapply(container, f, ...)
  } else {
    result <- lapply(container, function (x) lapplyn(x, f, n - 1, ...))
  }
  return (result)
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

get_parts_long <- function (prescriptions) {
  data <- do.call(c, prescriptions)
  parts <- list()
  
  for (total_momentum_idx in 1:length(data)) {
    total_momentum_str <- names(data)[[total_momentum_idx]]
    irreps <- data[[total_momentum_idx]]
    
    total_momentum_vec <- momentum_str_to_vec(total_momentum_str)
    total_momentum_sq <- sum(total_momentum_vec^2)
    
    for (irrep_idx in 1:length(irreps)) {
      irrep <- names(irreps)[[irrep_idx]]
      irrep_cols <- irreps[[irrep_idx]]
      if (length(irrep_cols) == 0) {
        next
      }
      
      for (irrep_col_idx in 1:length(irrep_cols)) {
        irrep_col_str <- names(irrep_cols)[[irrep_col_idx]]
        irrep_rows <- irrep_cols[[irrep_col_idx]]
        if (length(irrep_rows) == 0) {
          next
        }
        
        for (irrep_row_idx in 1:length(irrep_rows)) {
          irrep_row_str <- names(irrep_rows)[[irrep_row_idx]]
          corr_rows <- irrep_rows[[irrep_row_idx]]
          if (length(corr_rows) == 0) {
            next
          }
          
          for (corr_row_idx in 1:length(corr_rows)) {
            corr_row_str <- names(corr_rows)[[corr_row_idx]]
            corr_cols <- corr_rows[[corr_row_idx]]
            if (length(corr_cols) == 0) {
              next
            }
            
            for (corr_col_idx in 1:length(corr_cols)) {
              corr_col_str <- names(corr_cols)[[corr_col_idx]]
              corr <- corr_cols[[corr_col_idx]]
              
              part <- tibble(
                total_momentum_sq = total_momentum_sq,
                total_momentum_str = total_momentum_str,
                irrep = irrep,
                irrep_col = as.integer(irrep_col_str),
                irrep_row = as.integer(irrep_row_str),
                rel_momenta_str = corr_col_str)
              
              parts <- c(parts, list(part))
            }
          }
        }
      }
    }
  }
  
  reference_momenta_str <- c('000', '001', '110', '111', '002')
  filtered <- do.call(rbind, parts) %>%
    filter(total_momentum_str %in% reference_momenta_str) %>%
    distinct()
}

sort_momenta_vec_list <- function (momenta_vec) {
  momenta_str <- sapply(momenta_vec, momentum_vec_to_str, USE.NAMES = FALSE)
  lapply(sort(momenta_str), momentum_str_to_vec)
}

recursive_merge <- function (data, config_numbers) {
  first <- data[[1]]
  
  if (is.list(first)) {
    names_first <- names(first)
    res <- list()
    for (name in names_first) {
      elements <- lapply(data, function (elem) elem[[name]])
      merged <- recursive_merge(elements, config_numbers)
      res[[name]] <- merged
    }
    return (res)
  } else if (is.vector(first)) {
    res <- do.call(rbind, data)
    if (!missing(config_numbers)) {
      rownames(res) <- config_numbers
    }
    return (res)
  } else {
    stop()
  }
}

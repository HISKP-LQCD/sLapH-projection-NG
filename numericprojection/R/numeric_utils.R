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

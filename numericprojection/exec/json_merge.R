#!/usr/bin/env Rscript

library(dplyr)
library(numericprojection)
library(tibble)

all_projected_paths <- Sys.glob('projected/resolved_*_*_*.js')

parts <- do.call(rbind, lapply(all_projected_paths, split_filename))
colnames(parts) <- c('total_momentum', 'irrep', 'config_number')
df <- as_tibble(parts)
df$path <- all_projected_paths

make_merge <- function (stuff) {
  data <- lapply(stuff$path, jsonlite::read_json, simplifyVector = TRUE)
  merged <- recursive_merge(data, stuff$config_number)[[1]]
  
  dirpath <- sprintf('%s/correlator_matrices', workdir)
  if (!dir.exists(dirpath)) {
    dir.create(dirpath)
  }
                          
  tibble(merged = merged)
}

result <- df %>%
  arrange(total_momentum, irrep, config_number) %>%
  group_by(total_momentum, irrep) %>%
  do(make_merge(.))

merged2 <- result$merged %>%
  lapply(drop_empty)
l <- sapply(merged2, length)
result2 <- result[l > 0, ]

saveRDS(result2, 'correlator_matrices/corr_matrix.Rdata')

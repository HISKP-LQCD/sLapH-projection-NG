library(tibble)
library(dplyr)

source('numeric_utils.R')

split_filename <- function (path) {
  filename <- basename(path)
  parts <- stringr::str_match(filename, 'resolved_(.*)_(.*)_(.*)\\.js')
  parts[1, 2:4]
}

all_projected_paths <- Sys.glob('projected/resolved_*_*_*.js')

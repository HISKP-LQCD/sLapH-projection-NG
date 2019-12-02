#!/usr/bin/env Rscript

library(dplyr)
library(numericprojection)
library(tibble)

paths <- Sys.glob('prescriptions/prescription_*_*.js')
prescriptions <- lapply(paths, jsonlite::read_json)

res <- list()

for (prescription in prescriptions) {
  for (frame in names(prescription)) {
    if (is.null(res[[frame]])) {
      res[[frame]] <- list()
    }
    for (irrep in names(prescription[[frame]])) {
      if (length(unlist(prescription[[frame]][[irrep]])) > 0) {
      res[[frame]][[irrep]] <- prescription[[frame]][[irrep]]
      }
    }
  }
}

path_out <- '/prescriptions/all_prescriptions.js'
jsonlite::write_json(res, path_out)

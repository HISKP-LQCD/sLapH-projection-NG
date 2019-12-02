#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(magrittr)
library(numericprojection)
library(rhdf5)

args <- commandArgs(trailingOnly = TRUE)
cat('Command line arguments:\n  ')
print(args)

if (!exists('total_momentum')) {
  stopifnot(length(args) == 5)
  total_momentum <- as.integer(args[1:3])
  irrep <- args[4]
  config_number <- as.integer(args[5])
}

file_pattern <- sprintf('correlators/*_cnfg%04d.h5', config_number)
files <- Sys.glob(file_pattern)

diagrams <- sapply(files, function (file) strsplit(basename(file), '_')[[1]][1])
names(diagrams) <- NULL

files_list <- as.list(files)
names(files_list) <- diagrams

total_momentum_sq <- sum(total_momentum^2)
total_momentum_str <- paste0(sprintf('%d', total_momentum), collapse = '')

prescription_filename <- sprintf('prescriptions/prescription_%s_%s.js', total_momentum_str, irrep)
cat('Reading all prescriptions …\n')
all_prescriptions <- jsonlite::read_json(prescription_filename)
cat('  Done\n')

cat('Building list of needed dataset names …\n')
needed_names <- unique(unlist(lapplyn(all_prescriptions, function (rule) rule$datasetname, 7)))
cat('  Done\n')

cat('Opening HDF5 files …\n')
file_handles <- list()
for (diagram in diagrams) {
  path <- files_list[[diagram]]
  stopifnot(file.exists(path))
  print(path)
  file_handles[[diagram]] <- H5Fopen(path, flags = 'H5F_ACC_RDONLY')
}
cat('  Done\n')

load_dataset <- function (datasetname) {
    diagram <- strsplit(datasetname, '_')[[1]][1]
    
    skip_h5_errors <- FALSE
    
    if (skip_h5_errors) {
      dataset <- tryCatch(
        h5read(file_handles[[diagram]], datasetname),
        error = function (e) {
          warning(e)
          NA
        })
    } else {
      dataset <- h5read(file_handles[[diagram]], datasetname)
    }
    
    if (length(dataset) == 1 && is.na(dataset)) {
      return (0)
    }
    
    if (ncol(dataset) == 2) {
      return (dataset$re + 1i * dataset$im)
    } else if (colnames(dataset)[1] == 'rere') {
      return (dataset$rere - dataset$imim + 1i * dataset$reim + 1i * dataset$imre)
    } else {
      stop()
    }
}

cat('Loading correlators from HDF5 files …\n')
needed_raw <- lapply(needed_names, load_dataset)
names(needed_raw) <- needed_names
cat('  Done\n')

cat('Closing HDF5 files …\n')
invisible(lapply(file_handles, H5Fclose))
cat('  Done\n')

drop_zero_length <- function (x) {
  x[sapply(x, length) > 0]
}

filtered_prescriptions <- lapplyn(all_prescriptions, drop_zero_length, 5)
filtered_prescriptions <- lapplyn(filtered_prescriptions, drop_zero_length, 4)

resolve <- function (prescription) {
  if (length(prescription) == 0) {
    return (NA)
  }

  summands <- list(rep(NA, length(prescription)))

  for (i in 1:length(prescription)) {
    rule <- prescription[[i]]
    datasetname <- rule$datasetname
    diagram <- strsplit(datasetname, '_')[[1]][1]
    
    if (diagram == 'C4cV') {
      next
    }
    
    correlator <- needed_raw[[datasetname]]
    if (rule$conj) {
      correlator <- Conj(correlator)
    }
    weight <- rule$re + 1i * rule$im
    summands[[i]] <- Re(weight * correlator)
  }
  
  m <- do.call(rbind, na.omit(summands))
  apply(m, 2, sum)
}

cat('Projecting …\n')
resolved <- lapplyn(filtered_prescriptions, resolve, 6)
cat('  Done\n')

filtered <- resolved %>%
  lapplyn(drop_small, 5) %>%
  lapplyn(drop_empty, 4) %>%
  lapplyn(drop_empty, 3)

path <- 'projected'
if (!dir.exists(path)) {
  dir.create(path)
}
cat('Writing JSON …\n')
resolved_filename <- sprintf('%s/resolved_%s_%s_%04d.js', path, total_momentum_str, irrep, config_number)
jsonlite::write_json(filtered, resolved_filename, pretty = TRUE, digits = NA)
cat('  Done\n')

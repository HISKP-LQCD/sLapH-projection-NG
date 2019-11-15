#!/usr/bin/env Rscript

lapplyn <- function (container, f, n, ...) {
  if (n == 1) {
    result <- lapply(container, f, ...)
  } else {
    result <- lapply(container, function (x) lapplyn(x, f, n - 1, ...))
  }
  return (result)
}

print(getwd())

library(ggplot2)
library(dplyr)
library(magrittr)

theme_set(theme_light())

args <- commandArgs(trailingOnly = TRUE)
print(args)

if (!exists('config_number')) {
  stopifnot(length(args) == 1)
  config_number <- as.integer(args[1])
}

file_pattern <- sprintf('correlators/*_cnfg%04d.h5', config_number)
files <- Sys.glob(file_pattern)

if (!dir.exists('/storage/ueding/correlators')) {
    dir.create('/storage/ueding/correlators', recursive = TRUE)
}
temps <- sapply(files, function (x) tempfile(tmpdir = '/storage/ueding/correlators'))

print(files)
cat('Copying files to scratch …\n')
for (i in 1:length(files)) {
    cat(' ', files[i], '\n')
    file.copy(files[i], temps[i])
}

cat('Done.\n')

diagrams <- sapply(files, function (file) strsplit(basename(file), '_')[[1]][1])
names(diagrams) <- NULL

files_list <- as.list(temps)
names(files_list) <- diagrams

prescription_filename <- sprintf('prescriptions/all_prescriptions.js')
all_prescriptions <- jsonlite::read_json(prescription_filename)

needed_names <- unique(unlist(lapplyn(all_prescriptions, function (rule) rule$datasetname, 7)))

load_dataset <- function (datasetname) {
    diagram <- strsplit(datasetname, '_')[[1]][1]
    
    filename <- files_list[[diagram]]
    stopifnot(file.exists(filename))
    
    #print(sprintf('Fetching %s (a %s) from %s.\n', datasetname, diagram, filename))
    
    skip_h5_errors <- FALSE
    
    if (skip_h5_errors) {
      dataset <- tryCatch(
        rhdf5::h5read(filename, datasetname),
        error = function (e) {
          warning(e)
          NA
        })
    } else {
      dataset <- rhdf5::h5read(filename, datasetname)
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

needed_raw <- lapply(needed_names, load_dataset)

names(needed_raw) <- needed_names

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

resolved <- lapplyn(filtered_prescriptions, resolve, 6)

drop_small <- function (corr_list) {
  firsts <- sapply(corr_list, function (corr) corr[1])
  corr_list[abs(firsts) > 1.0e-8]
}

drop_empty <- function (l) {
  lengths <- sapply(l, length)
  l[lengths > 0]
}

filtered <- resolved %>%
  lapplyn(drop_small, 5) %>%
  lapplyn(drop_empty, 4) %>%
  lapplyn(drop_empty, 3)

path <- 'projected'
if (!dir.exists(path)) {
  dir.create(path)
}
resolved_filename <- sprintf('%s/all_resolved_%04d.js', path, config_number)
jsonlite::write_json(filtered, resolved_filename, pretty = TRUE, digits = NA)

for (i in length(temps)) {
    file.remove(temps[i])
}

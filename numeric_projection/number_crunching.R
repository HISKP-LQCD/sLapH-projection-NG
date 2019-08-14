#!/usr/bin/env Rscript

options(echo = TRUE)

print(getwd())

library(ggplot2)
library(dplyr)
library(magrittr)

theme_set(theme_light())

lapplyn <- function (container, f, n, ...) {
  if (n == 1) {
    result <- lapply(container, f, ...)
  } else {
    result <- lapply(container, function (x) lapplyn(x, f, n - 1, ...))
  }
  return (result)
}

args <- commandArgs(trailingOnly = TRUE)
print(args)

if (!exists('total_momentum')) {
  stopifnot(length(args) == 5)
  total_momentum <- as.integer(args[1:3])
  irrep <- args[4]
  config_number <- args[5]
}

total_momentum_sq <- sum(total_momentum^2)
total_momentum_str <- paste0(sprintf('%d', total_momentum), collapse = '')

total_momentum_ref <- (if (total_momentum_sq == 0) c(0, 0, 0)
else if (total_momentum_sq == 1) c(0, 0, 1)
else if (total_momentum_sq == 2) c(1, 1, 0)
else if (total_momentum_sq == 3) c(1, 1, 1)
else if (total_momentum_sq == 1) c(0, 0, 2))

prescription_filename <- sprintf('prescriptions/gevp-rho-%s-%s.js', total_momentum_str, irrep)
all_prescriptions <- jsonlite::read_json(prescription_filename)

needed_names <- unique(unlist(lapplyn(all_prescriptions, function (rule) rule$datasetname, 7)))

file_pattern <- sprintf('correlators/*_cnfg%d.h5', config_number)
files <- Sys.glob(file_pattern)

diagrams <- sapply(files, function (file) strsplit(basename(file), '_')[[1]][1])
names(diagrams) <- NULL

files_list <- as.list(files)
names(files_list) <- diagrams

needed_raw <- lapply(
  needed_names,
  function (datasetname) {
    #print(datasetname)
    diagram <- strsplit(datasetname, '_')[[1]][1]
    
    if (diagram == 'C4cV') {
      return (NA)
    }
    
    filename <- files_list[[diagram]]
    
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
      return (dataset$rere - 0*dataset$imim + 1i * dataset$reim + 1i * dataset$imre)
    } else {
      stop()
    }
  })

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

dir.create('projected')
resolved_filename <- sprintf('projected/resolved-rho-%s-%s-%04d.js', total_momentum_str, irrep, config_number)
jsonlite::write_json(resolved, resolved_filename, pretty = TRUE)

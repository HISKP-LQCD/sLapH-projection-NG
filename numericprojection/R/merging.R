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
    rownames(res) <- config_numbers
    return (res)
  } else {
    stop()
  }
}

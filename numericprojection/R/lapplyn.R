#' @export
drop_empty <- function (l) {
  lengths <- sapply(l, length)
  l[lengths > 0]
}

#' @export
lapplyn <- function (container, f, n, ...) {
  if (n == 1) {
    result <- lapply(container, f, ...)
  } else {
    result <- lapply(container, function (x) lapplyn(x, f, n - 1, ...))
  }
  return (result)
}

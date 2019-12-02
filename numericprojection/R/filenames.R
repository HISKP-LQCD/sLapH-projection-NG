#' @export
split_filename <- function (path) {
  filename <- basename(path)
  parts <- stringr::str_match(filename, 'resolved_(.*)_(.*)_(.*)\\.js')
  parts[1, 2:4]
}

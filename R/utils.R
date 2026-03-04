with_cache <- function(path, expr, force = FALSE) {
  if (file.exists(path) && !force) {
    return(readRDS(path))
  }
  result <- expr
  saveRDS(result, path)
  result
}

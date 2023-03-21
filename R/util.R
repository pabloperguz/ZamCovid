ZamCovid_file <- function(...) {
  system.file(..., package = "ZamCovid", mustWork = TRUE)
}

read_csv <- function(...) {
  read.csv(..., stringsAsFactors = FALSE)
}

`%||%` <- function(a, b) { # nolint
  if (is.null(a)) b else a
}

assert_proportion <- function(x, name = deparse(substitute(x))) {
  if (any(x < 0 | x > 1)) {
    stop(sprintf("'%s' must lie in [0, 1]", name),
         call. = FALSE)
  }
  invisible(x)
}

assert_non_negative <- function(x, name = deparse(substitute(x))) {
  if (any(x < 0)) {
    stop(sprintf("'%s' must have only non-negative values", name),
         call. = FALSE)
  }
  invisible(x)
}

assert_scalar <- function(x, name = deparse(substitute(x))) {
  if (length(x) != 1L) {
    stop(sprintf("'%s' must be a scalar", name), call. = FALSE)
  }
  invisible(x)
}

rename <- function(x, from, to, name = deparse(substitute(x))) {
  verify_names(x, required = from, allow_extra = TRUE, name = name)
  i <- match(from, names(x))
  names(x)[i] <- to
  x
}

verify_names <- function(x, required = NULL, optional = NULL,
                         allow_extra = FALSE,
                         name = deparse(substitute(x))) {
  nms <- names(x)
  if (anyDuplicated(nms)) {
    dups <- unique(nms[duplicated(nms)])
    stop(sprintf("Duplicate element names in '%s': %s",
                 name, paste(squote(dups), collapse = ", ")))
  }
  if (!allow_extra) {
    extra <- setdiff(nms, c(required, optional))
    if (length(extra) > 0) {
      stop(sprintf("Extra elements in '%s': %s",
                   name, paste(squote(extra), collapse = ", ")))
    }
  }
  msg <- setdiff(required, nms)
  if (length(msg) > 0) {
    stop(sprintf("Elements missing from '%s': %s",
                 name, paste(squote(msg), collapse = ", ")))
  }
  invisible(x)
}

squote <- function(x) {
  sprintf("'%s'", x)
}

vnapply <- function(.x, .fun, ...) {
  vapply(.x, .fun, numeric(1), ...)
}

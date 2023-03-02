zamcovid_file <- function(...) {
  system.file(..., package = "ZamCovid", mustWork = TRUE)
}

read_csv <- function(...) {
  read.csv(..., stringsAsFactors = FALSE)
}

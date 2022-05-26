# load all the files from the path
sensitivity_data_files <- list.files(
  path = system.file("sensitivity_data", package = "DAISIEprepExtra")
)
if (length(sensitivity_data_files) == 0) {
  stop("No results are in the results directory")
} else {
  file_paths <- as.list(paste0(
    system.file("sensitivity_data", package = "DAISIEprepExtra"),
    "/",
    sensitivity_data_files
  ))
  sensitivity_data <- lapply(file_paths, readRDS)
}


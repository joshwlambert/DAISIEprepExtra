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

sensitivity_dna <- lapply(sensitivity_data, "[[", "sensitivity_dna")
sensitivity_complete <- lapply(sensitivity_data, "[[", "sensitivity_complete")

sensitivity_dna <- Reduce(rbind, sensitivity_dna)
sensitivity_complete <- Reduce(rbind, sensitivity_complete)

sensitivity_data <- list(
  sensitivity_dna = sensitivity_dna,
  sensitivity_complete = sensitivity_complete
)

sensitivity <- DAISIEprep:::plot_sensitivity(
  sensitivity_data = sensitivity_data$sensitivity_dna
)

ggplot2::ggsave(
  plot = sensitivity,
  filename = file.path("inst", "plots", "sensitivity.png"),
  device = "png",
  width = 300,
  height = 200,
  units = "mm",
  dpi = 600
)

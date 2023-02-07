args <- commandArgs(TRUE)
args <- as.numeric(args)
parameter_index <- args[1]

set.seed(
  2639688,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

message("Parameter set: ", parameter_index)

performance_min <- DAISIEprep::benchmark(
  phylod = NULL,
  tree_size_range = c(10, 10000),
  num_points = 15,
  prob_on_island = c(0.2, 0.5),
  prob_endemic = c(0.2, 0.8),
  replicates = 10,
  extraction_method = "min",
  asr_method = NA,
  tie_preference = NA,
  log_scale = TRUE,
  parameter_index = parameter_index,
  verbose = TRUE
)

performance_asr <- DAISIEprep::benchmark(
  phylod = NULL,
  tree_size_range = c(10, 10000),
  num_points = 15,
  prob_on_island = c(0.2, 0.5),
  prob_endemic = c(0.2, 0.8),
  replicates = 10,
  extraction_method = "asr",
  asr_method = "parsimony",
  tie_preference = "island",
  log_scale = TRUE,
  parameter_index = parameter_index,
  verbose = TRUE
)

performance <- list(
  performance_min = performance_min,
  performance_asr = performance_asr
)

output_name <- paste0("performance_param_set_", parameter_index, ".rds")

output_folder <- file.path("results")

output_file_path <- file.path(output_folder, output_name)

saveRDS(object = performance, file = output_file_path)

message("Finished")

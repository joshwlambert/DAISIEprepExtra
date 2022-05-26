args <- commandArgs(TRUE)

args <- as.numeric(args)

parameter_index <- args[1]

parameter_space <- expand.grid(
  tree_size = c(10, 50, 100, 500, 1000, 5000, 10000),
  prob_on_island = c(0.2, 0.5),
  prob_endemic = c(0.2, 0.8)
)

set.seed(
  2639688,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

replicates <- 10

times_list <- list()

message("Parameter set: ", parameter_index)

mean_times_min <- c()
mean_times_asr <- c()
for (i in seq_len(replicates)) {

  message("Replicate: ", i, " of ", replicates)

  # simulate phylogeny
  phylo <- ape::rcoal(n = parameter_space$tree_size[parameter_index])

  # generate a set of unique tip labels that conform to standard
  tip_labels <- expand.grid(letters, letters, letters)
  tip_labels <- do.call(paste0, tip_labels)
  tip_labels <- tip_labels[1:parameter_space$tree_size[parameter_index]]
  tip_labels <- paste("bird", tip_labels, sep = "_")
  phylo$tip.label <- tip_labels

  prob_not_present <- 1 - parameter_space$prob_on_island[parameter_index]
  prob_endemic <-
    parameter_space$prob_endemic[parameter_index] * parameter_space$prob_on_island[parameter_index]
  prob_nonendemic <-
    (1 - parameter_space$prob_endemic[parameter_index]) * parameter_space$prob_on_island[parameter_index]


  empty_island <- TRUE
  while (empty_island) {
    # generate tip states under uniform sampling
    endemicity_status <- sample(
      x = c("endemic", "nonendemic", "not_present"),
      size = parameter_space$tree_size[parameter_index],
      replace = TRUE,
      prob = c(prob_endemic, prob_nonendemic, prob_not_present)
    )
    if (any(endemicity_status != "not_present")) {
      empty_island <- FALSE
    }
  }

  # add not present outgroup
  phylo <- DAISIEprep::add_outgroup(phylo)
  endemicity_status <- c("not_present", endemicity_status)

  # format data for DAISIEprep
  phylod <- phylobase::phylo4d(phylo, as.data.frame(endemicity_status))
  phylod <- DAISIEprep::add_asr_node_states(
    phylod = phylod,
    asr_method = "parsimony",
    tie_preference = "mainland",
    earliest_col = FALSE
  )


  # run extraction
  min_time <- system.time(for (n in 1:3) {
    island_tbl_min <- DAISIEprep::extract_island_species(
      phylod = phylod,
      extraction_method = "min",
      island_tbl = NULL,
      include_not_present = FALSE
    )
  })

  asr_time <- system.time(for (n in 1:3) {
    island_tbl_asr <- DAISIEprep::extract_island_species(
      phylod = phylod,
      extraction_method = "asr",
      island_tbl = NULL,
      include_not_present = FALSE
    )
  })

  mean_times_min[i] <- min_time["elapsed"] / 3
  mean_times_asr[i] <- asr_time["elapsed"] / 3
}
times_list <- list(
  min = mean_times_min,
  asr = mean_times_asr,
  parameter_index = parameter_index
)

output_name <- paste0("performance_param_set_", parameter_index, ".rds")

output_folder <- file.path("results")

output_file_path <- file.path(output_folder, output_name)

saveRDS(object = times_list, file = output_file_path)

message("Finished")

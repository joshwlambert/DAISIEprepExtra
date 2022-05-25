# load madagascar mammals species data table
data("madagascar_mammals", package = "DAISIEprepExtra")

# load the DNA only and complete posterior distribution of trees
phylos_dna <- ape::read.nexus(file = system.file(
  "extdata/Upham_dna_posterior_100.nex",
  package = "DAISIEprepExtra"
))
phylos_complete <- ape::read.nexus(file = system.file(
  "extdata/Upham_complete_posterior_100.nex",
  package = "DAISIEprepExtra"
))

# convert trees to phylo4 objects
phylos_dna <- lapply(phylos_dna, phylobase::phylo4)
phylos_complete <- lapply(phylos_complete, phylobase::phylo4)

# create endemicity status data frame
endemicity_status_dna <- lapply(
  phylos_dna,
  DAISIEprep::create_endemicity_status,
  island_species = madagascar_mammals
)

endemicity_status_complete <- lapply(
  phylos_complete,
  DAISIEprep::create_endemicity_status,
  island_species = madagascar_mammals
)

# combine tree and endemicity status
multi_phylods_dna <- list()
multi_phylods_complete <- list()
for (i in seq_along(phylos_dna)) {
  message("Converting phylo ", i, " of ", length(phylos_dna))
  multi_phylods_dna[[i]] <- phylobase::phylo4d(
    phylos_dna[[i]],
    endemicity_status_dna[[i]]
  )
  multi_phylods_complete[[i]] <- phylobase::phylo4d(
    phylos_complete[[i]],
    endemicity_status_complete[[i]]
  )
}

multi_phylods_dna <- lapply(
  multi_phylods_dna,
  DAISIEprep::add_asr_node_states,
  asr_method = "parsimony"
)

multi_phylods_complete <- lapply(
  multi_phylods_complete,
  DAISIEprep::add_asr_node_states,
  asr_method = "parsimony"
)

times_list <- list()
median_times_min_dna <- c()
median_times_min_complete <- c()
median_times_asr_dna <- c()
median_times_asr_complete <- c()
for (i in seq_along(multi_phylods_dna)) {

  message("Extracting ", i, " of ", length(multi_phylods_dna))

  # run extraction
  min_time_dna <- microbenchmark::microbenchmark(
    island_tbl_min <- DAISIEprep::extract_island_species(
      phylod = multi_phylods_dna[[i]],
      extraction_method = "min",
      island_tbl = NULL,
      include_not_present = FALSE
    ),
    times = 10L
  )

  min_time_dna <- microbenchmark::microbenchmark(
    island_tbl_min <- DAISIEprep::extract_island_species(
      phylod = multi_phylods_complete[[i]],
      extraction_method = "min",
      island_tbl = NULL,
      include_not_present = FALSE
    ),
    times = 10L
  )

  asr_time_dna <- microbenchmark::microbenchmark(
    island_tbl_asr <- DAISIEprep::extract_island_species(
      phylod = multi_phylods_dna[[i]],
      extraction_method = "asr",
      island_tbl = NULL,
      include_not_present = FALSE
    ),
    times = 10L
  )

  asr_time_complete <- microbenchmark::microbenchmark(
    island_tbl_asr <- DAISIEprep::extract_island_species(
      phylod = multi_phylods_complete[[i]],
      extraction_method = "asr",
      island_tbl = NULL,
      include_not_present = FALSE
    ),
    times = 10L
  )

  median_times_min_dna[i] <- median(min_time_dna$time)
  median_times_min_complete[i] <- median(min_time_complete$time)
  median_times_asr_dna[i] <- median(asr_time_dna$time)
  median_times_asr_complete[i] <- median(asr_time_complete$time)
}

results <- data.frame(
  median_times_min_dna = median_times_min_dna,
  median_times_min_complete = median_times_min_complete,
  median_times_asr_dna = median_times_asr_dna,
  median_times_asr_complete = median_times_asr_complete
)

output_name <- "performance_empirical.rds"

output_folder <- file.path("results")

output_file_path <- file.path(output_folder, output_name)

saveRDS(object = results, file = output_file_path)

message("Finished")


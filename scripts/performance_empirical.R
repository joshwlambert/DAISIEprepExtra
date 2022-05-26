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

mean_times_min_dna <- c()
mean_times_min_complete <- c()
mean_times_asr_dna <- c()
mean_times_asr_complete <- c()
for (i in seq_along(multi_phylods_dna)) {

  message("Extracting ", i, " of ", length(multi_phylods_dna))

  # run extraction
  min_time_dna <- system.time(for (n in 1:3) {
    island_tbl_min <- DAISIEprep::extract_island_species(
      phylod = multi_phylods_dna[[i]],
      extraction_method = "min",
      island_tbl = NULL,
      include_not_present = FALSE
    )
  })

  min_time_complete <- system.time(for (n in 1:3) {
    island_tbl_min <- DAISIEprep::extract_island_species(
      phylod = multi_phylods_complete[[i]],
      extraction_method = "min",
      island_tbl = NULL,
      include_not_present = FALSE
    )
  })

  asr_time_dna <- system.time(for (n in 1:3) {
    island_tbl_asr <- DAISIEprep::extract_island_species(
      phylod = multi_phylods_dna[[i]],
      extraction_method = "asr",
      island_tbl = NULL,
      include_not_present = FALSE
    )
  })

  asr_time_complete <- system.time(for (n in 1:3) {
    island_tbl_asr <- DAISIEprep::extract_island_species(
      phylod = multi_phylods_complete[[i]],
      extraction_method = "asr",
      island_tbl = NULL,
      include_not_present = FALSE
    )
  })

  mean_times_min_dna[i] <- min_time_dna["elapsed"] / 3
  mean_times_min_complete[i] <- min_time_complete["elapsed"] / 3
  mean_times_asr_dna[i] <- asr_time_dna["elapsed"] / 3
  mean_times_asr_complete[i] <- asr_time_complete["elapsed"] / 3
}

results <- data.frame(
  mean_times_min_dna = mean_times_min_dna,
  mean_times_min_complete = mean_times_min_complete,
  mean_times_asr_dna = mean_times_asr_dna,
  mean_times_asr_complete = mean_times_asr_complete
)

output_name <- "performance_empirical.rds"

output_folder <- file.path("results")

output_file_path <- file.path(output_folder, output_name)

saveRDS(object = results, file = output_file_path)

message("Finished")


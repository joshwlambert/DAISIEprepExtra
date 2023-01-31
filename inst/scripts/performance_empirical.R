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

times_min_dna <- list()
times_min_complete <- list()
times_asr_dna <- list()
times_asr_complete <- list()
for (i in seq_along(multi_phylods_dna)) {

  times_min_dna[[i]] <- DAISIEprep::benchmark(
    phylod = multi_phylods_dna[[i]],
    tree_size_range = NA,
    num_points = NA,
    prob_on_island = NA,
    prob_endemic = NA,
    replicates = 1,
    extraction_method = "min",
    asr_method = NA,
    tie_preference = NA,
    verbose = TRUE
  )

  times_min_complete[[i]] <- DAISIEprep::benchmark(
    phylod = multi_phylods_complete[[i]],
    tree_size_range = NA,
    num_points = NA,
    prob_on_island = NA,
    prob_endemic = NA,
    replicates = 1,
    extraction_method = "min",
    asr_method = NA,
    tie_preference = NA,
    verbose = TRUE
  )

  times_asr_dna[[i]] <- DAISIEprep::benchmark(
    phylod = multi_phylods_dna[[i]],
    tree_size_range = NA,
    num_points = NA,
    prob_on_island = NA,
    prob_endemic = NA,
    replicates = 1,
    extraction_method = "asr",
    asr_method = "parsimony",
    tie_preference = "island",
    verbose = TRUE
  )

  times_asr_complete[[i]] <- DAISIEprep::benchmark(
    phylod = multi_phylods_complete[[i]],
    tree_size_range = NA,
    num_points = NA,
    prob_on_island = NA,
    prob_endemic = NA,
    replicates = 1,
    extraction_method = "asr",
    asr_method = "parsimony",
    tie_preference = "island",
    verbose = TRUE
  )
}

performance <- list(
  times_min_dna = times_min_dna,
  times_min_complete = times_min_complete,
  times_asr_dna = times_asr_dna,
  times_asr_complete = times_asr_complete
)

output_name <- "performance_empirical.rds"

output_folder <- file.path("results")

output_file_path <- file.path(output_folder, output_name)

saveRDS(object = performance, file = output_file_path)

message("Finished")


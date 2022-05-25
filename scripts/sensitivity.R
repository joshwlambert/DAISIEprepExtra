parameter_space <- data.frame(
  extraction_method = "min",
  asr_method = NA_character_,
  tie_preference = NA_character_
)

parameter_space <- rbind(
  parameter_space,
  expand.grid(
    extraction_method = c("asr"),
    asr_method = c("parsimony", "mk"),
    tie_preference = c("island", "mainland")
  )
)

# load madagascar mammals species data table
data("madagascar_mammals", package = "DAISIEprepExtra")

# load the DNA only and complete posterior distribution of trees
phylos_dna <- ape::read.nexus(file = system.file(
  "inst/extdata/Upham_dna_posterior_100.nex",
  package = "DAISIEprepExtra"
))
phylos_complete <- ape::read.nexus(file = system.file(
  "inst/extdata/Upham_complete_posterior_100.nex",
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

ml_list <- list()
for (i in seq_len(nrow(parameter_space))) {

  message("Parameter set: ", i)

  if (parameter_space$extraction_method[i] == "asr") {
    multi_phylods_dna <- lapply(
      multi_phylods_dna,
      DAISIEprep::add_asr_node_states,
      asr_method = parameter_space$asr_method[i]
    )

    multi_phylods_complete <- lapply(
      multi_phylods_complete,
      DAISIEprep::add_asr_node_states,
      asr_method = parameter_space$asr_method[i]
    )
  }

  # extract island community
  multi_island_tbl_dna <- DAISIEprep::multi_extract_island_species(
    multi_phylod = multi_phylods_dna,
    extraction_method = parameter_space$extraction_method[i],
    verbose = TRUE
  )

  multi_island_tbl_complete <- DAISIEprep::multi_extract_island_species(
    multi_phylod = multi_phylods_complete,
    extraction_method = parameter_space$extraction_method[i],
    verbose = TRUE
  )

  # convert to daisie data table
  daisie_datatable_dna <- lapply(
    multi_island_tbl_dna,
    DAISIEprep::as_daisie_datatable,
    island_age = 88
  )

  daisie_datatable_complete <- lapply(
    multi_island_tbl_complete,
    DAISIEprep::as_daisie_datatable,
    island_age = 88
  )

  # convert to daisie data list
  daisie_data_list_dna <- lapply(
    daisie_datatable_dna,
    DAISIEprep::create_daisie_data,
    island_age = 88,
    num_mainland_species = 1000
  )

  daisie_data_list_complete <- lapply(
    daisie_datatable_complete,
    DAISIEprep::create_daisie_data,
    island_age = 88,
    num_mainland_species = 1000
  )

  # fit DAISIE model to data
  ml_dna <- list()
  ml_complete <- list()
  for (j in seq_along(phylos_dna)) {
    message(
      "Fitting DAISIE model to data set ",
      j,
      " of ",
      length(daisie_data_list_dna)
    )

    ml_dna[[j]] <- DAISIE::DAISIE_ML_CS(
      datalist = daisie_data_list_dna[[j]],
      initparsopt = c(1, 1, 100, 0.1, 1),
      idparsopt = 1:5,
      parsfix = NULL,
      idparsfix = NULL,
      ddmodel = 11,
      jitter = 1e-5
    )

    ml_complete[[j]] <- DAISIE::DAISIE_ML_CS(
      datalist = daisie_data_list_complete[[j]],
      initparsopt = c(1, 1, 100, 0.1, 1),
      idparsopt = 1:5,
      parsfix = NULL,
      idparsfix = NULL,
      ddmodel = 11,
      jitter = 1e-5
    )
  }

  ml_list[[i]] <- list(dna = ml_dna, complete = ml_complete)
}

output_name <- "sensitivity.rds"

output_folder <- file.path("results")

output_file_path <- file.path(output_folder, output_name)

saveRDS(object = ml_list, file = output_file_path)

message("Finished")


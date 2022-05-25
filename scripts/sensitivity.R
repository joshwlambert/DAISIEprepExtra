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
load(file = system.file(
  "data/madagascar_mammals.rda",
  package = "DAISIEprepExtra"
))

# load the DNA only and complete posterior distribution of trees
dna_phylos <- ape::read.nexus(file = system.file(
  "inst/extdata/Upham_dna_posterior_100.nex",
  package = "DAISIEprepExtra"
))
complete_phylos <- ape::read.nexus(file = system.file(
  "inst/extdata/Upham_complete_posterior_100.nex",
  package = "DAISIEprepExtra"
))

# convert trees to phylo4 objects
dna_phylos <- lapply(dna_phylos, phylobase::phylo4)
complete_phylos <- lapply(complete_phylos, phylobase::phylo4)

# create endemicity status data frame
endemicity_status_dna <- lapply(
  dna_phylos,
  DAISIEprep::create_endemicity_status,
  island_species = madagascar_mammals
)

endemicity_status_complete <- lapply(
  complete_phylos,
  DAISIEprep::create_endemicity_status,
  island_species = madagascar_mammals
)

# combine tree and endemicity status
dna_multi_phylods <- list()
complete_multi_phylods <- list()
for (i in seq_along(dna_phylos)) {
  message("Converting phylo ", i, " of ", length(dna_phylos))
  dna_multi_phylods[[i]] <- phylobase::phylo4d(
    dna_phylos[[i]],
    endemicity_status_dna[[i]]
  )
  complete_multi_phylods[[i]] <- phylobase::phylo4d(
    complete_phylos[[i]],
    endemicity_status_complete[[i]]
  )
}

ml_list <- list()
for (i in seq_len(nrow(parameter_space))) {

  message("Parameter set: ", i)

  if (parameter_space$extraction_method[i] == "asr") {
    dna_multi_phylods <- lapply(
      dna_multi_phylods,
      DAISIEprep::add_asr_node_states,
      asr_method = parameter_space$asr_method[i]
    )

    complete_multi_phylods <- lapply(
      complete_multi_phylods,
      DAISIEprep::add_asr_node_states,
      asr_method = parameter_space$asr_method[i]
    )
  }

  # extract island community
  dna_multi_island_tbl <- DAISIEprep::multi_extract_island_species(
    multi_phylod = dna_multi_phylods,
    extraction_method = parameter_space$extraction_method[i],
    verbose = TRUE
  )
  complete_multi_island_tbl <- DAISIEprep::multi_extract_island_species(
    multi_phylod = complete_multi_phylods,
    extraction_method = parameter_space$extraction_method[i],
    verbose = TRUE
  )

  # convert to daisie data table
  daisie_datatable_dna <- lapply(
    dna_multi_island_tbl,
    DAISIEprep::as_daisie_datatable,
    island_age = 88
  )

  daisie_datatable_complete <- lapply(
    complete_multi_island_tbl,
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
  for (j in seq_along(dna_phylos)) {
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

output_name <- paste0("sensitivity.rds")

output_folder <- file.path("results")

output_file_path <- file.path(output_folder, output_name)

saveRDS(object = ml_list, file = output_file_path)

message("Finished")


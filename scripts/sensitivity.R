library(DAISIEprep)

parameter_space <- data.frame(
  extraction_method = "min",
  asr_method = NA_character_,
  tie_preference = NA_character_,
  earliest_col = NA_character_
)

parameter_space <- rbind(
  parameter_space,
  expand.grid(
    extraction_method = c("asr"),
    asr_method = c("parsimony", "mk"),
    tie_preference = c("island", "mainland"),
    earliest_col = c(TRUE, FALSE)
  )
)

# load madagascar mammals species data table
load(file = system.file(
  "data/madagascar_mammals.rda",
  package = "DAISIEprepExtra"
))

# load the DNA only and complete trees
dna_phylo <- ape::read.nexus(file = system.file(
  "inst/extdata/Upham_dna_mcc.tre",
  package = "DAISIEprepExtra"
))
complete_phylo <- ape::read.nexus(file = system.file(
  "inst/extdata/Upham_complete_mcc.tre",
  package = "DAISIEprepExtra"
))

ml_list <- list()
for (i in seq_len(nrow(parameter_space))) {

  message("Parameter set: ", i)

  # convert trees to phylo4 objects
  dna_phylo <- phylobase::phylo4(dna_phylo)
  complete_phylo <- phylobase::phylo4(complete_phylo)

  # create endemicity status data frame
  endemicity_status_dna <- DAISIEprep::create_endemicity_status(
    phylo = dna_phylo,
    island_species = madagascar_mammals
  )
  endemicity_status_complete <- DAISIEprep::create_endemicity_status(
    phylo = complete_phylo,
    island_species = madagascar_mammals
  )

  # combine tree and endemicity status
  dna_phylod <- phylobase::phylo4d(dna_phylo, endemicity_status_dna)
  complete_phylod <- phylobase::phylo4d(complete_phylo, endemicity_status_complete)

  if (parameter_space$extraction_method == "asr") {
    dna_phylod <- add_asr_node_states(
      phylod = phylod,
      asr_method = parameter_space$asr_method[i]
    )
    complete_phylod <- add_asr_node_states(
      phylod = phylod,
      asr_method = parameter_space$asr_method[i]
    )
  }

  # extract island community using min algorithm
  dna_island_tbl <- DAISIEprep::extract_island_species(
    phylod = dna_phylod,
    extraction_method = parameter_space$extraction_method[i]
  )
  complete_island_tbl <- DAISIEprep::extract_island_species(
    phylod = complete_phylod,
    extraction_method = parameter_space$extraction_method[i]
  )

  # convert to daisie data table
  daisie_datatable_dna <- DAISIEprep::as_daisie_datatable(
    island_tbl = dna_island_tbl,
    island_age = 88
  )
  daisie_datatable_complete <- DAISIEprep::as_daisie_datatable(
    island_tbl = complete_island_tbl,
    island_age = 88
  )

  # convert to daisie data list
  daisie_data_list_dna <- DAISIEprep::create_daisie_data(
    daisie_datatable = daisie_datatable_dna,
    island_age = 88,
    num_mainland_species = 1000
  )
  daisie_data_list_complete <- DAISIEprep::create_daisie_data(
    daisie_datatable = daisie_datatable_complete,
    island_age = 88,
    num_mainland_species = 1000
  )

  ml_dna <- DAISIE::DAISIE_ML_CS(
    datalist = daisie_data_list_dna,
    initparsopt = c(1, 1, 100, 0.1, 1),
    idparsopt = 1:5,
    parsfix = NULL,
    idparsfix = NULL,
    ddmodel = 11,
    jitter = 1e-5
  )

  ml_complete <- DAISIE::DAISIE_ML_CS(
    datalist = daisie_data_list_complete,
    initparsopt = c(1, 1, 100, 0.1, 1),
    idparsopt = 1:5,
    parsfix = NULL,
    idparsfix = NULL,
    ddmodel = 11,
    jitter = 1e-5
  )
  ml_list[[i]] <- list(ml_dna, ml_complete)
}


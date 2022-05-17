# load madagascar mammals species data table
load(file = system.file(
  "data/madagascar_mammals.rda",
  package = "DAISIEprepExtra"
))

# load the DNA only and complete trees
mammal_posterior_dna <- ape::read.nexus(system.file(
  "inst/extdata/Upham_dna_posterior_100.nex",
  package = "DAISIEprepExtra"
))

mammal_posterior_complete <- ape::read.nexus(system.file(
  "inst/extdata/Upham_complete_posterior_100.nex",
  package = "DAISIEprepExtra"
))

# convert trees to phylo4 objects
dna_phylos <- lapply(mammal_posterior_dna, phylobase::phylo4)
complete_phylos <- lapply(mammal_posterior_complete, phylobase::phylo4)

# remove phylo objects to free up some memory
rm(mammal_posterior_dna)
rm(mammal_posterior_complete)
gc()

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

# ensure dna_phylos and complete_phylos are the same length
length(dna_phylos) == length(complete_phylos)

# combine tree and endemicity status
dna_multi_phylods <- list()
complete_multi_phylods <- list()
for (i in seq_along(dna_phylos)) {
  message("Converting phylo ", i, " of ", length(dna_phylos))
  dna_multi_phylods[[i]] <- phylobase::phylo4d(dna_phylos[[i]], endemicity_status_dna[[i]])
  complete_multi_phylods[[i]] <- phylobase::phylo4d(complete_phylos[[i]], endemicity_status_complete[[i]])
}

# extract island community using min algorithm
multi_island_tbl_dna <- DAISIEprep::multi_extract_island_species(
  multi_phylod = dna_multi_phylods,
  extraction_method = "min"
)

multi_island_tbl_complete <- DAISIEprep::multi_extract_island_species(
  multi_phylod = complete_multi_phylods,
  extraction_method = "min"
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


# ensure dna data list and complete data list are the same length
length(daisie_data_list_dna) == length(daisie_data_list_complete)

# combine tree and endemicity status
ml_dna <- list()
ml_complete <- list()
for (i in seq_along(dna_phylos)) {
  message(
    "Fitting DAISIE model to data set ", i, " of ", length(daisie_data_list_dna)
  )

  ml_dna[[i]] <- DAISIE::DAISIE_ML_CS(
    datalist = daisie_data_list_dna[[i]],
    initparsopt = c(1, 1, 100, 0.1, 1),
    idparsopt = 1:5,
    parsfix = NULL,
    idparsfix = NULL,
    ddmodel = 11,
    jitter = 1e-5,
    optimmethod = "simplex"
  )

  ml_complete[[i]] <- DAISIE::DAISIE_ML_CS(
    datalist = daisie_data_list_complete[[i]],
    initparsopt = c(1, 1, 100, 0.1, 1),
    idparsopt = 1:5,
    parsfix = NULL,
    idparsfix = NULL,
    ddmodel = 11,
    jitter = 1e-5,
    optimmethod = "simplex"
  )
}

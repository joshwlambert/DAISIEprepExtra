library(DAISIEprep)
library(DAISIEprepExtra)

# It is useful to know from the outset the requirements for this scripts
# it uses the madagascar_mammals dataset stored in the data/ folder. This is a
# two column data frame with the species names (tip labels), formatted as
# genus_species in the first column and the species endemicity status, formatted
# as all lowercase underscore separated (use DAISIEprep::translate_status()) to
# convert endemicity status to the correct format. Column names must be
# "tip_labels" and "tip_endemicity_status"
# here is the head of the madagascar mammals data
#                  tip_labels tip_endemicity_status
# 1         Echinops_telfairi               endemic
# 2            Geogale_aurita               endemic
# 3    Hemicentetes_nigriceps               endemic
# 4 Hemicentetes_semispinosus               endemic
# 5    Microgale_brevicaudata               endemic
# 6          Microgale_cowani               endemic
# the phylogeny of mammal species is also required and here we use the
# posterior distribution of 100 mammal phylogenetic trees from Upham et al.,
# (2019) https://doi.org/10.1371/journal.pbio.3000494 stored as
# Upham_dna_posterior_100.tre and Upham_complete_posterior_100.tre files loaded
# from extdata/
# the variables required to be given by the user for this script are:
# the island age (island_age) in million years, and number of species in the
# mainland pool (num_mainland_species)
# We set these at the start for our estimates for the mammals of Madagascar
island_age <- 88
num_mainland_species <- 1000

# load madagascar mammals species data table
data("madagascar_mammals", package = "DAISIEprepExtra")

# load the DNA only and complete trees
mammal_posterior_dna <- ape::read.nexus(system.file(
  "extdata/Upham_dna_posterior_100.nex",
  package = "DAISIEprepExtra"
))

mammal_posterior_complete <- ape::read.nexus(system.file(
  "extdata/Upham_complete_posterior_100.nex",
  package = "DAISIEprepExtra"
))

# for this example we subset the posterior distribution of 100 phylogenetic
# trees to 5 trees to make the script faster to run. If you want to run the
# script for the entire set of 100 trees you can not run the next two lines
mammal_posterior_dna <- mammal_posterior_dna[1:5]
mammal_posterior_complete <- mammal_posterior_complete[1:5]

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
  dna_multi_phylods[[i]] <- phylobase::phylo4d(
    dna_phylos[[i]],
    endemicity_status_dna[[i]]
  )
  complete_multi_phylods[[i]] <- phylobase::phylo4d(
    complete_phylos[[i]],
    endemicity_status_complete[[i]]
  )
}

# extract island community using min algorithm
multi_island_tbl_dna <- DAISIEprep::multi_extract_island_species(
  multi_phylod = dna_multi_phylods,
  extraction_method = "min",
  verbose = TRUE
)

multi_island_tbl_complete <- DAISIEprep::multi_extract_island_species(
  multi_phylod = complete_multi_phylods,
  extraction_method = "min",
  verbose = TRUE
)

# convert to daisie data table
daisie_datatable_dna <- lapply(
  multi_island_tbl_dna,
  DAISIEprep::as_daisie_datatable,
  island_age = island_age
)

daisie_datatable_complete <- lapply(
  multi_island_tbl_complete,
  DAISIEprep::as_daisie_datatable,
  island_age = island_age
)

# convert to daisie data list
daisie_data_list_dna <- lapply(
  daisie_datatable_dna,
  DAISIEprep::create_daisie_data,
  island_age = island_age,
  num_mainland_species = num_mainland_species
)

daisie_data_list_complete <- lapply(
  daisie_datatable_complete,
  DAISIEprep::create_daisie_data,
  island_age = island_age,
  num_mainland_species = num_mainland_species
)

# ensure dna data list and complete data list are the same length
length(daisie_data_list_dna) == length(daisie_data_list_complete)

# to run the maximum likelihood DAISIE inference model we now have the data
# required by the model, but now need the model settings to ensure we are
# fitting the correct model
# the initparsopt is the initial parameter estimates that are used by the
# optimisation algorithm as a starting point from which to optimise. In an
# ideal situation the maximum likelihood parameter estimates and model
# likelihood would be independent of the starting conditions by finding the
# global optimum. However, this is not always the case so we pick sensible
# values that should prevent us getting trapped in local optima. To avoid local
# optima in empirical analyses run the model with different initial parameter
# values to check whether a global optimum has likely been found. The parameters
# in the vector are, in order: cladogenesis, extinction, carrying capacity,
# colonisation and anagenesis
# idparsopt is selecting which parameters in the we want to optimise and here
# we chosen the five parameters listed above and thus idparsopt = 1:5
# parsfix and idparsfix are the parameters that are to be fixed and thus not
# optimised, however, here we optimise all five parameters and so these are both
# NULL
# ddmodel is set to 11 to set both cladogenesis and colonisation to diversity-
# dependent given the carrying capacity parameter.
# the jitter is set to 1e-5 to help initial likelihood calculations but is not
# important for the model fitting and can be left as 1e-5 for all model runs

# run the DAISIE model for all data sets
ml_dna <- list()
ml_complete <- list()
for (i in seq_along(dna_phylos)) {
  message(
    "Fitting DAISIE model to data set ", i, " of ", length(daisie_data_list_dna)
  )

  ml_dna[[i]] <- DAISIE::DAISIE_ML_CS(
    datalist = daisie_data_list_dna[[i]],
    initparsopt = c(1, 1, 200, 0.1, 1),
    idparsopt = 1:5,
    parsfix = NULL,
    idparsfix = NULL,
    ddmodel = 11,
    jitter = 1e-5
  )

  ml_complete[[i]] <- DAISIE::DAISIE_ML_CS(
    datalist = daisie_data_list_complete[[i]],
    initparsopt = c(1, 1, 200, 0.1, 1),
    idparsopt = 1:5,
    parsfix = NULL,
    idparsfix = NULL,
    ddmodel = 11,
    jitter = 1e-5
  )
}

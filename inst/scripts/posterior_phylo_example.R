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

# load the checklist
checklist <- utils::read.csv(
  file = system.file(
    "extdata", "madagascar_mammal_checklist.csv",
    package = "DAISIEprepExtra"
  ),
  header = TRUE,
)

missing_species <- DAISIEprep::count_missing_species(
  checklist = checklist,
  phylo_name_col = "Name_In_Tree",
  genus_name_col = "Genus",
  in_phylo_col = "Sampled",
  endemicity_status_col = "DAISIE_Status_Species",
  rm_species_col = "Remove_Species"
)

# load the complete tree
complete_phylos <- ape::read.nexus(
  system.file(
    "extdata", "Upham_complete_posterior_100.nex",
    package = "DAISIEprepExtra"
  )
)

# for this example we subset the posterior distribution of 100 phylogenetic
# trees to 3 trees to make the script faster to run. If you want to run the
# script for the entire set of 100 trees you can not run the next line
complete_phylos <- complete_phylos[1:3]

# convert trees to phylo4 objects
complete_phylos <- lapply(complete_phylos, phylobase::phylo4)

# create endemicity status data frame
endemicity_status <- lapply(
  complete_phylos,
  DAISIEprep::create_endemicity_status,
  island_species = madagascar_mammals
)

# combine tree and endemicity status
multi_phylod <- mapply(
  phylobase::phylo4d,
  complete_phylos,
  endemicity_status
)

# reconstruct geographic ancestral states for extraction with asr
multi_phylod <- lapply(
  multi_phylod,
  DAISIEprep::add_asr_node_states,
  asr_method = "mk",
  tie_preference = "mainland"
)

# extract island community using asr algorithm
multi_island_tbl <- DAISIEprep::multi_extract_island_species(
  multi_phylod = multi_phylod,
  extraction_method = "asr",
  verbose = TRUE
)

# determine which island clade the missing species should be assigned to
missing_genus <- lapply(
  multi_island_tbl,
  DAISIEprep::unique_island_genera
)

# add missing species that match genera found in the island tbl
multi_island_tbl <- mapply(
  DAISIEprep::add_multi_missing_species,
  missing_genus,
  multi_island_tbl,
  missing_species = list(missing_species)
)

# remove missing species that have already been inserted into the island tbl
no_island_tbl_missing_species <- mapply(
  DAISIEprep::rm_multi_missing_species,
  missing_genus,
  multi_island_tbl,
  missing_species = list(missing_species),
  SIMPLIFY = FALSE
)

# add the Archaeoindris, Babakotia, Hadropithecus, Mesopropithecus, Pachylemur
# as a missing species of the clade with
# Megaladapis_edwardsi in it as it is a extinct lemur species
multi_island_tbl <- lapply(
  multi_island_tbl,
  DAISIEprep::add_missing_species,
  num_missing_species = 8,
  species_name = "Megaladapis_edwardsi"
)

# add the Plesiorycteropus as a missing_species of the clade with
# Tenrec_ecaudatus in it as it is a tenrec species
multi_island_tbl <- lapply(
  multi_island_tbl,
  DAISIEprep::add_missing_species,
  num_missing_species = 2,
  species_name = "Tenrec_ecaudatus"
)

# add Chaerephon species as separate colonisation
Chaerephon_stem_age <- list()
for (i in seq_along(multi_phylod)) {
  Chaerephon_stem_age[[i]] <- DAISIEprep::extract_stem_age(
    genus_name = "Chaerephon",
    phylod = multi_phylod[[i]],
    stem = "genus"
  )
}

# add the Chaerephon as an stem age max age given the stem age in the tree
multi_island_tbl <- mapply(
  DAISIEprep::add_island_colonist,
  multi_island_tbl,
  clade_name = list("Chaerephon_leucogaster"),
  status = list("nonendemic"),
  missing_species = list(0),
  col_time = Chaerephon_stem_age,
  col_max_age = list(TRUE),
  branching_times = list(NA_real_),
  min_age = list(NA_real_),
  species = list("Chaerephon_leucogaster"),
  clade_type = list(1)
)

# add the Macronycteris as a missing species of the clade with
# Hipposideros_commersoni in it as it is a bat species
multi_island_tbl <- lapply(
  multi_island_tbl,
  DAISIEprep::add_missing_species,
  num_missing_species = 2,
  species_name = "Hipposideros_commersoni"
)


# convert to daisie data list
daisie_data_list <- lapply(
  multi_island_tbl,
  DAISIEprep::create_daisie_data,
  island_age = island_age,
  num_mainland_species = num_mainland_species
)

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
mls <- list()
for (i in seq_along(multi_phylod)) {
  message(
    "Fitting DAISIE model to data set ", i, " of ", length(daisie_data_list)
  )

  mls[[i]] <- DAISIE::DAISIE_ML_CS(
    datalist = daisie_data_list[[i]],
    initparsopt = c(1, 1, 200, 0.1, 1),
    idparsopt = 1:5,
    parsfix = NULL,
    idparsfix = NULL,
    ddmodel = 11,
    jitter = 1e-5
  )
}

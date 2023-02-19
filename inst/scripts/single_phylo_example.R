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
# Upham_dna_mcc.tre and Upham_complete_mcc.tre files loaded from extdata/
# See Upham et al., (2019) https://doi.org/10.1371/journal.pbio.3000494
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
complete_phylo <- ape::read.nexus(
  file = system.file(
    "extdata", "Upham_complete_mcc.tre",
    package = "DAISIEprepExtra"
  )
)

# remove tags after species name in tip labels
complete_phylo$tip.label <- gsub(
  pattern = "([A-Za-z]+_[A-Za-z]+).*",
  replacement = "\\1",
  x = complete_phylo$tip.label
)

# convert tree to phylo4 objects
complete_phylo <- phylobase::phylo4(complete_phylo)

# create endemicity status data frame
endemicity_status <- DAISIEprep::create_endemicity_status(
  phylo = complete_phylo,
  island_species = madagascar_mammals
)

# combine tree and endemicity status
phylod <- phylobase::phylo4d(complete_phylo, endemicity_status)

# reconstruct geographic ancestral states for extraction with asr
phylod <- DAISIEprep::add_asr_node_states(
  phylod = phylod,
  asr_method = "mk",
  tie_preference = "mainland"
)

# extract island community using asr algorithm
island_tbl <- DAISIEprep::extract_island_species(
  phylod = phylod,
  extraction_method = "asr"
)

# determine which island clade the missing species should be assigned to
missing_genus <- DAISIEprep::unique_island_genera(island_tbl = island_tbl)

# add missing species that match genera found in the island tbl
island_tbl <- DAISIEprep::add_multi_missing_species(
  missing_species = missing_species,
  missing_genus = missing_genus,
  island_tbl = island_tbl
)

# remove missing species that have already been inserted into the island tbl
no_island_tbl_missing_species <- DAISIEprep::rm_multi_missing_species(
  missing_species = missing_species,
  missing_genus = missing_genus,
  island_tbl = island_tbl
)

# add the Archaeoindris, Babakotia, Hadropithecus, Mesopropithecus, Pachylemur
# as a missing species of the clade with
# Megaladapis_edwardsi in it as it is a extinct lemur species
island_tbl <- DAISIEprep::add_missing_species(
  island_tbl = island_tbl,
  num_missing_species = 8,
  species_to_add_to = "Megaladapis_edwardsi"
)

# add the Plesiorycteropus as a missing_species of the clade with
# Tenrec_ecaudatus in it as it is a tenrec species
island_tbl <- DAISIEprep::add_missing_species(
  island_tbl = island_tbl,
  num_missing_species = 2,
  species_to_add_to = "Tenrec_ecaudatus"
)

# add Chaerephon species as separate colonisation
Chaerephon_stem_age <- DAISIEprep::extract_stem_age( # nolint
  genus_name = "Chaerephon",
  phylod = phylod,
  stem = "genus"
)

# add the Chaerephon as an stem age max age given the stem age in the tree
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Chaerephon_leucogaster",
  status = "nonendemic",
  missing_species = 0,
  col_time = Chaerephon_stem_age,
  col_max_age = TRUE,
  branching_times = NA_real_,
  min_age = NA_real_,
  species = "Chaerephon_leucogaster",
  clade_type = 1
)

# add the Macronycteris as a missing species of the clade with
# Hipposideros_commersoni in it as it is a bat species
island_tbl <- DAISIEprep::add_missing_species(
  island_tbl = island_tbl,
  num_missing_species = 2,
  species_to_add_to = "Hipposideros_commersoni"
)

# convert to daisie data list
daisie_data_list <- DAISIEprep::create_daisie_data(
  data = island_tbl,
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

ml <- DAISIE::DAISIE_ML_CS(
  datalist = daisie_data_list,
  initparsopt = c(1, 1, 200, 0.1, 1),
  idparsopt = 1:5,
  parsfix = NULL,
  idparsfix = NULL,
  ddmodel = 11,
  jitter = 1e-5
)

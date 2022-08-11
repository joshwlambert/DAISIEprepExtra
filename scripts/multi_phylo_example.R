library(DAISIEprep)
library(DAISIEprepExtra)

# It is useful to know from the outset the requirements for this scripts
# it uses the hawaii_asters dataset stored in the data/ folder. This is a
# two column data frame with the species names (tip labels), formatted as
# genus_species in the first column and the species endemicity status, formatted
# as all lowercase underscore separated (use DAISIEprep::translate_status()) to
# convert endemicity status to the correct format. Column names must be
# "tip_labels" and "tip_endemicity_status"
# here is the head of the Hawaiian asteraceae data
#             tip_labels tip_endemicity_status
# 1 Adenostemma_viscosum            nonendemic
# 2 Artemisia_kauaiensis               endemic
# 3  Artemisia_mauiensis               endemic
# 4  Artemisia_australis               endemic
# 5    Bidens_amplectens               endemic
# 6   Bidens_asymmetrica               endemic
# multiple phylogenies of asteraceae species are also required and here we use
# the silverswords phylogenetic tree from Landis et al. (2018)
# https://doi.org/10.1111/evo.13594 stored as Landis_2018_g4_mcc.tre, the Bidens
# phylogenetic tree from Knope et al. (2020) https://doi.org/10.1111/jse.12704
# stored as Pacific_Bidens_Knope_2020.tre, and the Hesperomannia phylogenetic
# tree from Keeley et al. (2021) https://doi:10.1002/ajb2.1614 stored as
# Keeley_2021.tre, all files are stored in extdata/
# the variables required to be given by the user for this script are:
# the island age (island_age) in million years, and number of species in the
# mainland pool (num_mainland_species)
# We set these at the start for our estimates for the Asteraceae of Hawaii
island_age <- 6.15
num_mainland_species <- 7500

# load Hawaiian asteraceae species data table
data("hawaii_asters", package = "DAISIEprepExtra")

# load the checklist
checklist <- utils::read.csv(
  file = system.file(
    "extdata", "hawaii_asteraceae_checklist.csv",
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
  rm_species_col = NULL
)

# load the phylogenies
hesperomannia <- ape::read.nexus(
  file = system.file(
    "extdata", "Keeley_2021.tre",
    package = "DAISIEprepExtra",
    mustWork = TRUE
  )
)

silversword <- ape::read.nexus(
  file = system.file(
    "extdata", "Landis_2018_g4_mcc.tre",
    package = "DAISIEprepExtra",
    mustWork = TRUE
  )
)

# remove subspecies from the silversword to not inflate species richness on the
# archipelago, subspecies can be optionally kept if deemed important for
# the analysis of the diversification rate
split_names <- strsplit(x = silversword$tip.label, split = "_")
genus_names <- sapply(split_names, "[[", 1)
species_names <- sapply(split_names, "[[", 2)
genus_species_names <- paste(genus_names, species_names, sep = "_")
island_multi_tip <- genus_species_names[which(duplicated(genus_species_names))]
island_multi_tip <- unique(island_multi_tip)
tip_position <- list()
for (i in seq_along(island_multi_tip)) {
  tip_position[[i]] <- grep(
    pattern = island_multi_tip[i],
    x = silversword$tip.label
  )
}

drop_random_subspecies <- unlist(lapply(tip_position, function(x) {
  sample(x = x, size = length(x) - 1, replace = FALSE)
}))

silversword <- ape::drop.tip(
  phy = silversword,
  tip = drop_random_subspecies
)

# load biden tree
bidens <- ape::read.nexus(
  file = system.file(
    "extdata", "Pacific_Bidens_Knope_2020.tre",
    package = "DAISIEprepExtra",
    mustWork = TRUE
  )
)

# deal with unnamed species
# there is one unnamed species from the tree on Hawaii (BidspnovLF12810G90)
# there are another six species that are that are unnamed and are not on Hawaii
# these can be given temporary, fictional names for the purposes of naming
# standardisation but will not be included in the island data set once extracted
# name is Bid- rather than Bidens_- because it will be changed below with the
# other tips in the tree
bidens$tip.label[grep(pattern = "Bidsp", x = bidens$tip.label)] <- c(
  "Bidbutaud",
  "Bidstarbuckis_a",
  "Bidstarbuckis_b",
  "Bidtaputuarai",
  "Bidspa",
  "Bidspb",
  "Bidspnovoahu"
)

# rename tip labels to conform to DAISIEprep
# split the tip labels by capital letters as they differentiate the genus name
# from the collector name and collector tag
species <- gsub('([[:upper:]])', ' \\1', bidens$tip.label)

# split names by white space
species_split <- strsplit(x = species, split = " ")

# keep only the genus and species name
species_name <- lapply(species_split, "[[", 2)

# rename groups by identifiers
species_name <- gsub(
  pattern = "Bid",
  replacement = "Bidens_",
  x = species_name
)
species_name <- gsub(
  pattern = "Coreop",
  replacement = "Coreopsis_",
  x = species_name
)
species_name <- gsub(
  pattern = "Cos",
  replacement = "Cosmos_",
  x = species_name
)
species_name <- gsub(
  pattern = "Fit",
  replacement = "Fitchia_",
  x = species_name
)
species_name <- gsub(
  pattern = "Opa",
  replacement = "Oparanthus_",
  x = species_name
)
species_name <- gsub(
  pattern = "Coreoc",
  replacement = "Coreocarpus_",
  x = species_name
)
species_name <- gsub(
  pattern = "Hen",
  replacement = "Hendricksonia_",
  x = species_name
)

# rename those species that do not fit identifiers
species_name[grep(pattern = "Glossocardia", x = species_name)] <-
  "Glossocardia_bidens"
species_name[grep(pattern = "Corcyclocarpa", x = species_name)] <-
  "Coreopsis_cyclocarpa"
species_name[grep(pattern = "Corveticde", x = species_name)] <-
  "Coreopsis_verticillata"

# remove any numbers left in the species names
species_name <- gsub(pattern = "[0-9]", replacement = "", x = species_name)

bidens$tip.label <- species_name

# remove multiple samples of the same species, some of these do not form
# monophyletic species so for this example we choose one of the species samples
# in the tree at random and drop the remaining samples (tips). This also stops
# species richness being inflated, as with the silverswords above, by only having
# one tip in the tree per species
island_multi_tip <- bidens$tip.label[which(duplicated(bidens$tip.label))]
island_multi_tip <- unique(island_multi_tip)
tip_position <- list()
for (i in seq_along(island_multi_tip)) {
  tip_position[[i]] <- grep(
    pattern = island_multi_tip[i],
    x = bidens$tip.label
  )
}

drop_random_subspecies <- unlist(lapply(tip_position, function(x) {
  sample(x = x, size = length(x) - 1, replace = FALSE)
}))

bidens <- ape::drop.tip(
  phy = bidens,
  tip = drop_random_subspecies
)

# convert trees to phylo4 objects
hesperomannia <- phylobase::phylo4(hesperomannia)
silversword <- phylobase::phylo4(silversword)
bidens <- phylobase::phylo4(bidens)

# create endemicity status data frame
endemicity_status_hesperomannia <- DAISIEprep::create_endemicity_status(
  phylo = hesperomannia,
  island_species = hawaii_asters
)

endemicity_status_silversword <- DAISIEprep::create_endemicity_status(
  phylo = silversword,
  island_species = hawaii_asters
)

endemicity_status_bidens <- DAISIEprep::create_endemicity_status(
  phylo = bidens,
  island_species = hawaii_asters
)

# combine tree and endemicity status
hesperomannia_phylod <- phylobase::phylo4d(
  hesperomannia,
  endemicity_status_hesperomannia
)

silversword_phylod <- phylobase::phylo4d(
  silversword,
  endemicity_status_silversword
)

bidens_phylod <- phylobase::phylo4d(
  bidens,
  endemicity_status_bidens
)

# extract island community using min algorithm
island_tbl <- DAISIEprep::extract_island_species(
  phylod = hesperomannia_phylod,
  extraction_method = "min"
)

island_tbl <- DAISIEprep::extract_island_species(
  phylod = silversword_phylod,
  extraction_method = "min",
  island_tbl = island_tbl
)

island_tbl <- DAISIEprep::extract_island_species(
  phylod = bidens_phylod,
  extraction_method = "min",
  island_tbl = island_tbl
)

# add Artemisia as an endemic_MaxAge; we have the stem age (3.93 Ma) and crown
# (1.45 Ma) from in text
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Artemisia",
  status = "endemic",
  missing_species = 1,
  col_time = 3.93,
  col_max_age = TRUE,
  branching_times = 1.45,
  min_age = NA_real_,
  species = c(
    "Artemisia_kauaiensis",
    "Artemisia_mauiensis",
    "Artemisia_australis"
  ),
  clade_type = 1
)

# add Lipochaeta-Melanthera alliance, in-text stem age (1.26 Ma)
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Lipochaeta-Melanthera",
  status = "endemic",
  missing_species = 19,
  col_time = 1.26,
  col_max_age = TRUE,
  branching_times = NA_real_,
  min_age = NA_real_,
  species = c(
    "Lipochaeta_connata_connata",
    "Lipochaeta_degeneri",
    "Lipochaeta_heterophylla",
    "Lipochaeta_lobata_lobata",
    "Lipochaeta_rockii",
    "Lipochaeta_succulenta",
    "Melanthera_bryanii",
    "Melanthera_fauriei",
    "Melanthera_integrifolia",
    "Melanthera_kamolensis",
    "Melanthera_lavarum",
    "Melanthera_micrantha_micrantha",
    "Melanthera_perdita",
    "Melanthera_populifolia",
    "Melanthera_remyi",
    "Melanthera_subcordata",
    "Melanthera_tenuifolia",
    "Melanthera_tenuis",
    "Melanthera_venosa",
    "Melanthera_waimeaensis"
  ),
  clade_type = 1
)

# add Pseudognaphalium
# max age (phylogeny in paper but no estimates of colonisations available)
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Pseudognaphalium",
  status = "endemic",
  missing_species = 0,
  col_time = NA_real_,
  col_max_age = TRUE,
  branching_times = NA_real_,
  min_age = NA_real_,
  species = "Pseudognaphalium_sandwicensium",
  clade_type = 1
)


# add Tetramolopium
# diversity data from taxonomic sources
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Tetramolopium",
  status = "endemic",
  missing_species = 10,
  col_time = NA_real_,
  col_max_age = TRUE,
  branching_times = NA_real_,
  min_age = NA,
  species = c(
    "Tetramolopium_humile",
    "Tetramolopium_lepidotum",
    "Tetramolopium_remyi",
    "Tetramolopium_rockii",
    "Tetramolopium_arenarium",
    "Tetramolopium_capillare",
    "Tetramolopium_consanguineum",
    "Tetramolopium_conyzoides",
    "Tetramolopium_filiforme",
    "Tetramolopium_sylvae",
    "Tetramolopium_tenerrimum"
  ),
  clade_type = 1
)

# add Keysseria
# diversity data from taxonomic sources
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Keysseria",
  status = "endemic",
  missing_species = 2,
  col_time = NA_real_,
  col_max_age = TRUE,
  branching_times = NA_real_,
  min_age = NA_real_,
  species = c(
    "Keysseria_maviensis",
    "Keysseria_erici",
    "Keysseria_helena"
  ),
  clade_type = 1
)

# add Remya
# diversity data from taxonomic sources
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Remya",
  status = "endemic",
  missing_species = 2,
  col_time = NA_real_,
  col_max_age = TRUE,
  branching_times = NA_real_,
  min_age = NA,
  species = c(
    "Remya_mauiensis",
    "Remya_kauaiensis",
    "Remya_montgomeryi"
  ),
  clade_type = 1
)

# add Adenostemma
# diversity data from taxonomic sources
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Adenostemma",
  status = "nonendemic",
  missing_species = 0,
  col_time = NA_real_,
  col_max_age = TRUE,
  branching_times = NA_real_,
  min_age = NA_real_,
  species = "Adenostemma_viscosum",
  clade_type = 1
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
  initparsopt = c(1, 1, 100, 0.1, 1),
  idparsopt = 1:5,
  parsfix = NULL,
  idparsfix = NULL,
  ddmodel = 11,
  jitter = 1e-5
)

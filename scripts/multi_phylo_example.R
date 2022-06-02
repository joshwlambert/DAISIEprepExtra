library(DAISIEprep)
library(DAISIEprepExtra)

# load Hawaiian asteraceae species data table
data("hawaii_asters")

# load the DNA only and complete trees
hesperomannia <- ape::read.nexus(file = system.file(
  "extdata/Keeley_2021.tre",
  package = "DAISIEprepExtra"
))

silversword <- ape::read.nexus(file = system.file(
  "extdata/Landis_2018_g4_mcc.tre",
  package = "DAISIEprepExtra"
))

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
bidens <- ape::read.nexus(file = system.file(
  "extdata/Pacific_Bidens_Knope_2020.tre",
  package = "DAISIEprepExtra"
))

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

# add Artemisia
# In-text example
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Artemisia",
  status = "endemic",
  missing_species = 3,
  branching_times = c(6.15),
  min_age = 1.45
)

# add Lipochaeta-Melanthera alliance
# In-text
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Lipochaeta-Melanthera",
  status = "endemic",
  missing_species = 22,
  branching_times = c(1.26),
  min_age = NA
)

# add Pseudognaphalium
# max age (phylogeny in paper but no estimates of colonisations available)
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Pseudognaphalium",
  status = "endemic",
  missing_species = 1,
  branching_times = c(6.15),
  min_age = NA
)

# add Tetramolopium
# diversity data from taxonomic sources
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Tetramolopium",
  status = "endemic",
  missing_species = 11,
  branching_times = c(6.15),
  min_age = NA
)

# add Keysseria
# diversity data from taxonomic sources
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Keysseria",
  status = "endemic",
  missing_species = 3,
  branching_times = c(6.15),
  min_age = NA
)

# add Remya
# diversity data from taxonomic sources
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Remya",
  status = "endemic",
  missing_species = 3,
  branching_times = c(6.15),
  min_age = NA
)

# add Adenostemma
# diversity data from taxonomic sources
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Adenostemma",
  status = "nonendemic",
  missing_species = 0,
  branching_times = c(6.15),
  min_age = NA
)

# convert to daisie data table
daisie_datatable <- DAISIEprep::as_daisie_datatable(
  island_tbl = island_tbl,
  island_age = 6.15
)

# convert to daisie data list
daisie_data_list <- DAISIEprep::create_daisie_data(
  daisie_datatable = daisie_datatable,
  island_age = 6.15,
  num_mainland_species = 1000
)

ml <- DAISIE::DAISIE_ML_CS(
  datalist = daisie_data_list,
  initparsopt = c(1, 1, 100, 0.1, 1),
  idparsopt = 1:5,
  parsfix = NULL,
  idparsfix = NULL,
  ddmodel = 11,
  jitter = 1e-5
)


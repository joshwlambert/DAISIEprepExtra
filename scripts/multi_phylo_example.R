# load Hawaiian asteraceae species data table
load(file = system.file(
  "data/hawaii_asters.rda",
  package = "DAISIEprepExtra"
))

# load the DNA only and complete trees
hesperomannia <- ape::read.nexus(file = system.file(
  "inst/extdata/Keeley_2021.tre",
  package = "DAISIEprepExtra"
))

silversword <- ape::read.nexus(file = system.file(
  "inst/extdata/Landis_2018_g4_mcc.tre",
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

# load biden tree when available

# convert trees to phylo4 objects
hesperomannia <- phylobase::phylo4(hesperomannia)
silversword <- phylobase::phylo4(silversword)

# create endemicity status data frame
endemicity_status_hesperomannia <- DAISIEprep::create_endemicity_status(
  phylo = hesperomannia,
  island_species = hawaii_asters
)

endemicity_status_silversword <- DAISIEprep::create_endemicity_status(
  phylo = silversword,
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
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Artemisia",
  status = "endemic",
  missing_species = 3,
  branching_times = c(6.15),
  min_age = 1.45
)

# add Lipochaeta-Melanthera alliance
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Lipochaeta-Melanthera",
  status = "endemic",
  missing_species = 22,
  branching_times = c(1.26),
  min_age = NA
)

# add Pseudognaphalium
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Pseudognaphalium",
  status = "endemic",
  missing_species = 1,
  branching_times = c(6.15),
  min_age = NA
)

# add Tetramolopium
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Tetramolopium",
  status = "endemic",
  missing_species = 11,
  branching_times = c(6.15),
  min_age = NA
)

# add Keysseria
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Keysseria",
  status = "endemic",
  missing_species = 3,
  branching_times = c(6.15),
  min_age = NA
)

# add Remya
island_tbl <- DAISIEprep::add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Remya",
  status = "endemic",
  missing_species = 3,
  branching_times = c(6.15),
  min_age = NA
)

# add Adenostemma
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


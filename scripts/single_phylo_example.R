# load madagascar mammals species data table
data("madagascar_mammals", package = "DAISIEprepExtra")

# load the DNA only and complete trees
dna_phylo <- ape::read.nexus(file = system.file(
  "inst/extdata/Upham_dna_mcc.tre",
  package = "DAISIEprepExtra"
))
complete_phylo <- ape::read.nexus(file = system.file(
  "inst/extdata/Upham_complete_mcc.tre",
  package = "DAISIEprepExtra"
))

# convert trees to phylo4 objects
dna_phylo <- phylobase::phylo4(dna_phylo)
complete_phylo <- phylobase::phylo4(complete_phylo)

# create endemicity status data frame
endemicity_status <- DAISIEprep::create_endemicity_status(
  phylo = dna_phylo,
  island_species = madagascar_mammals
)

# combine tree and endemicity status
phylod <- phylobase::phylo4d(dna_phylo, endemicity_status)

# extract island community using min algorithm
island_tbl <- DAISIEprep::extract_island_species(
  phylod = phylod,
  extraction_method = "min"
)

# convert to daisie data table
daisie_datatable <- DAISIEprep::as_daisie_datatable(
  island_tbl = island_tbl,
  island_age = 88
)

# convert to daisie data list
daisie_data_list <- DAISIEprep::create_daisie_data(
  daisie_datatable = daisie_datatable,
  island_age = 88,
  num_mainland_species = 1000
)

ml <- DAISIE::DAISIE_ML_CS(
  datalist = daisie_data_list,
  initparsopt = c(1, 1, 200, 0.1, 1),
  idparsopt = 1:5,
  parsfix = NULL,
  idparsfix = NULL,
  ddmodel = 11,
  jitter = 1e-5
)


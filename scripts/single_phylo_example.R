# load madagascar mammals species data table
load("data/madagascar_mammals.rda")

# load the DNA only and complete trees
dna_phylo <- ape::read.nexus(file = "data/Upham_dna_mcc.tre")
complete_phylo <- ape::read.nexus(file = "data/Upham_complete_mcc.tre")

# convert trees to phylo4 objects
dna_phylo <- phylobase::phylo4(dna_phylo)
complete_phylo <- phylobase::phylo4(complete_phylo)

# create endemicity status data frame
endemicity_status <- DAISIEprep::create_endemicity_status(
  phylo = dna_phylo,
  island_species = madagascar_mammals
)

phylod <- phylobase::phylo4d(dna_phylo, endemicity_status)

island_tbl <- DAISIEprep::extract_island_species(
  phylod = phylod,
  extraction_method = "min"
)

daisie_datatable <- DAISIEprep::as_daisie_datatable(
  island_tbl = island_tbl,
  island_age = 88
)

daisie_data_list <- DAISIEprep::create_daisie_data(
  daisie_datatable = daisie_datatable,
  island_age = 88,
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

simplex <- DAISIE::DAISIE_ML_CS(
  datalist = daisie_data_list,
  initparsopt = c(1, 1, 100, 0.1, 1),
  idparsopt = 1:5,
  parsfix = NULL,
  idparsfix = NULL,
  ddmodel = 11,
  optimmethod = "simplex"
)


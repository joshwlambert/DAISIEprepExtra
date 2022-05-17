# load Hawaiian asteraceae species data table
hawaii_asters <- read.csv(file = system.file(
  "inst/extdata/hawaii_asteraceae_species.csv",
  package = "DAISIEprepExtra"
))

aster_tip_labels <- paste(hawaii_asters$Genus, hawaii_asters$Species, sep = "_")

aster_endemicity_status <- c()
for (i in seq_len(nrow(hawaii_asters))) {
  aster_endemicity_status[i] <- DAISIEprep::translate_status(
    hawaii_asters$Endemic[i]
  )

}
hawaii_asters <- data.frame(
  tip_labels = aster_tip_labels,
  tip_endemicity_status =aster_endemicity_status
)

# load the DNA only and complete trees
hesperomannia <- ape::read.nexus(file = system.file(
  "inst/extdata/Keeley_2021.tre",
  package = "DAISIEprepExtra"
))

silversword <- ape::read.nexus(file = system.file(
  "inst/extdata/Landis_2018_g4_mcc.tre",
  package = "DAISIEprepExtra"
))

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
hesperomannia_phylod <- phylobase::phylo4d(hesperomannia, endemicity_status_hesperomannia)
silversword_phylod <- phylobase::phylo4d(silversword, endemicity_status_silversword)

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


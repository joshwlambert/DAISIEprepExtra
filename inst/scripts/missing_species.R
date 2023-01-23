## Script to give an example of each function that influences the number of
## missing species in an island-community data set

# set seed to control for stochasticity
set.seed(1)

# simulate phylogeny and check its ultrametric
tree <- ape::rphylo(n = 100, birth = 1, death = 0)
ape::is.ultrametric(tree)

# rename tip label in line with expectations of DAISIEprep
all_species <- gsub(
  pattern = "^t",
  replacement = "species_",
  x = tree$tip.label
)
tree$tip.label <- all_species

# sample species to remove from the tree to be missing species
missing_species <- sample(x = tree$tip.label, size = 15, replace = FALSE)

# prune tree of missing species
tree <- ape::drop.tip(phy = tree, tip = missing_species)

# sample which species are on the island
island_tip_labels <- sample(all_species, size = 25, replace = FALSE)

# work out whether a species is sampled in the phylogeny
phylo_sampled <- island_tip_labels %in% tree$tip.label

island_species <- data.frame(
  genus = "species",
  species = gsub(pattern = "^species_", replacement = "", x = island_tip_labels),
  tip_labels = island_tip_labels,
  tip_endemicity_status = sample(
    c("endemic", "nonendemic"), size = 25, replace = TRUE),
  sampled = phylo_sampled
)

# add an outgroup species to the tree
tree <- DAISIEprep::add_outgroup(tree)

endemicity_status <- DAISIEprep::create_endemicity_status(
  phylo = tree,
  island_species = island_species[, c("tip_labels", "tip_endemicity_status")]
)

# convert tree to phylo4 objects
tree <- phylobase::phylo4(tree)

# combine tree and endemicity status
phylod <- phylobase::phylo4d(tree, endemicity_status)

island_tbl <- DAISIEprep::extract_island_species(
  phylod = phylod,
  extraction_method = "min"
)

# count missing species for each genera
DAISIEprep::count_missing_species(
  checklist = island_species,
  phylo_name_col = "tip_labels",
  genus_name_col = "genus",
  in_phylo_col = "sampled",
  endemicity_status_col = "tip_endemicity_status"
)

# determine which island clade the missing species should be assigned to
missing_genus <- DAISIEprep::unique_island_genera(island_tbl = island_tbl)

# add missing species that match genera found in the island tbl
island_tbl <- DAISIEprep::add_multi_missing_species(
  missing_species = missing_species,
  missing_genus = missing_genus,
  island_tbl = island_tbl
)

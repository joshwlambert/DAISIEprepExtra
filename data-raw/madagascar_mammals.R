## code to prepare `madagascar_mammals` dataset goes here

tbl <- utils::read.csv(
  file = system.file(
    "extdata", "madagascar_mammal_checklist.csv",
    package = "DAISIEprepExtra"
  ),
  header = TRUE
)

# some species may be included in the checklist because they are present
# on Madagascar, but are introduced by humans, or are migratory species
# that do not comply with the assumptions of rare dispersal to and from the
# island, therefore these species can be removed automatically using the
# rm_species column
rm_species <- which(tbl$Remove_Species)
if (length(rm_species) > 0) {
  tbl <- tbl[-rm_species, ]
}

# subset to name in the phylogeny and the endemicity status
# use endemicity status corrected for DAISIE
island_species <- tbl[, c("Name_In_Tree", "DAISIE_Status_Species")]

names(island_species) <- c("tip_labels", "tip_endemicity_status")

# remove species that have NA in the data frame (most likely from not being in
# the phylogeny)
island_species <- stats::na.omit(island_species)

# replace white space with underscores
name_in_tree <- gsub(
  pattern = " ",
  replacement = "_",
  x = island_species$tip_labels
)
island_species$tip_labels <- name_in_tree


status_species <- c()
for (i in seq_along(island_species$tip_endemicity_status)) {
  status_species[i] <- DAISIEprep::translate_status(
    island_species$tip_endemicity_status[i]
  )
}
island_species$tip_endemicity_status <- status_species

madagascar_mammals <- island_species

usethis::use_data(madagascar_mammals, overwrite = TRUE)

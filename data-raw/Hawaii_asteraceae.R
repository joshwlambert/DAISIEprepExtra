## code to prepare `Hawaii_asteraceae` dataset goes here

tbl <- read.csv(
  file = system.file(
    "inst/extdata/hawaii_asteraceae_species.csv",
    package = "DAISIEprepExtra"
  ),
  header = TRUE
)

names <- paste(tbl$Genus, tbl$Species, sep = "_")
status_species <- tbl$Endemic

hawaii_asters <- data.frame(name = names, status_species = status_species)

status_species <- c()
for (i in seq_along(hawaii_asters$name)) {
  status_species[i] <- DAISIEprep::translate_status(
    hawaii_asters$status_species[i]
  )
}
hawaii_asters$status_species <- status_species

names(hawaii_asters) <- c("tip_labels", "tip_endemicity_status")

missing_madagascar_mammals <- tbl[, c("Genus", "Species", "In_Upham")]

missing_madagascar_mammals <- missing_species_tbl[which(missing_species_tbl$In_Upham == "No"), ]

usethis::use_data(hawaii_asters, overwrite = TRUE)
usethis::use_data(missing_madagascar_mammals, overwrite = TRUE)

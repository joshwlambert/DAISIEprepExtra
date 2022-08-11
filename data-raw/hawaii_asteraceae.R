## code to prepare `hawaii_asteraceae` dataset goes here

tbl <- utils::read.csv(
  file = system.file(
    "extdata", "hawaii_asteraceae_checklist.csv",
    package = "DAISIEprepExtra"
  ),
  header = TRUE
)

names <- paste(tbl$Genus, tbl$Species, sep = "_")
status_species <- tbl$DAISIE_status

hawaii_asters <- data.frame(name = names, status_species = status_species)

status_species <- c()
for (i in seq_along(hawaii_asters$name)) {
  status_species[i] <- DAISIEprep::translate_status(
    hawaii_asters$status_species[i]
  )
}
hawaii_asters$status_species <- status_species

names(hawaii_asters) <- c("tip_labels", "tip_endemicity_status")

usethis::use_data(hawaii_asters, overwrite = TRUE)


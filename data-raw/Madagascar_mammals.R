## code to prepare `Madagascar_mammals` dataset goes here

tbl <- read.csv(
  file = system.file(
    "inst/extdata/madagascar_mammal_species.csv",
    package = "DAISIEprepExtra"
  ),
  header = TRUE,
  sep = ";"
)

madagascar_mammals <- tbl[, c("Name_in_Upham", "Status_Species")]
madagascar_mammals <- na.omit(madagascar_mammals)
name_in_upham <- gsub(pattern = " ", replacement = "_", x = madagascar_mammals$Name_in_Upham)
madagascar_mammals$Name_in_Upham <- name_in_upham

status_species <- c()
for (i in seq_along(madagascar_mammals$Status_Species)) {
  status_species[i] <- DAISIEprep::translate_status(madagascar_mammals$Status_Species)
}
madagascar_mammals$Status_Species <- status_species

names(madagascar_mammals) <- c("tip_labels", "tip_endemicity_status")



missing_madagascar_mammals <- tbl[, c("Genus", "Species", "In_Upham")]

missing_madagascar_mammals <- missing_species_tbl[which(missing_species_tbl$In_Upham == "No"), ]

usethis::use_data(madagascar_mammals, overwrite = TRUE)
usethis::use_data(missing_madagascar_mammals, overwrite = TRUE)

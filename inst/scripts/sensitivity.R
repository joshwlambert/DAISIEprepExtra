# script to run the sensitivity analysis

# load madagascar mammals species data table
data("madagascar_mammals", package = "DAISIEprepExtra")

# load the DNA only and complete posterior distribution of trees
phylos_dna <- ape::read.nexus(
  file = system.file(
    "extdata", "Upham_dna_posterior_100.nex",
    package = "DAISIEprepExtra"
  )
)
phylos_complete <- ape::read.nexus(
  file = system.file(
    "extdata", "Upham_complete_posterior_100.nex",
    package = "DAISIEprepExtra"
  )
)

sensitivity_dna <- lapply(
  phylos_dna,
  DAISIEprep::sensitivity,
  island_species = madagascar_mammals,
  extraction_method = c("min", "asr"),
  asr_method = c("parsimony", "mk"),
  tie_preference = c("island", "mainland"),
  island_age = c(5),
  num_mainland_species = c(1000)
)

sensitivity_complete <- lapply(
  phylos_complete,
  DAISIEprep::sensitivity,
  island_species = madagascar_mammals,
  extraction_method = c("min", "asr"),
  asr_method = c("parsimony", "mk"),
  tie_preference = c("island", "mainland"),
  island_age = c(5),
  num_mainland_species = c(1000)
)

sensitivity <- list(
  sensitivity_dna = sensitivity_dna,
  sensitivity_complete = sensitivity = complete
)

output_name <- "sensitivity.rds"

output_folder <- file.path("results")

output_file_path <- file.path(output_folder, output_name)

saveRDS(object = sensitivity, file = output_file_path)

message("Finished")


# script to run the sensitivity analysis in parallel

# get argument passed from command line and store in tree_index
args <- commandArgs(TRUE)
args <- as.numeric(args)
tree_index <- args[1]

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

sensitivity_dna <- DAISIEprep::sensitivity(
  phylo = phylos_dna[[tree_index]],
  island_species = madagascar_mammals,
  extraction_method = c("min", "asr"),
  asr_method = c("parsimony", "mk"),
  tie_preference = c("mainland"),
  island_age = 88,
  num_mainland_species = c(1000),
  verbose = TRUE
)

sensitivity_complete <- DAISIEprep::sensitivity(
  phylo = phylos_dna[[tree_index]],
  island_species = madagascar_mammals,
  extraction_method = c("min", "asr"),
  asr_method = c("parsimony", "mk"),
  tie_preference = c("mainland"),
  island_age = 88,
  num_mainland_species = c(1000),
  verbose = TRUE
)

sensitivity <- list(
  sensitivity_dna = sensitivity_dna,
  sensitivity_complete = sensitivity_complete
)

output_name <- paste0("sensitivity_tree_index_", tree_index, ".rds")

output_folder <- file.path("results")

output_file_path <- file.path(output_folder, output_name)

saveRDS(object = sensitivity, file = output_file_path)

message("Finished")


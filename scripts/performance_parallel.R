parameter_space <- expand.grid(
  tree_size = c(10, 50, 100, 500, 1000, 5000, 10000),
  prob_on_island = c(0.2, 0.5),
  prob_endemic = c(0.2, 0.8)
)

set.seed(
  2639688,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

times_list <- list()
for (i in seq_len(nrow(parameter_space))) {

  message("Parameter set: ", i)

  median_times_min <- c()
  median_times_asr <- c()
  for (j in 1:10) {

    # simulate phylogeny
    phylo <- ape::rcoal(n = parameter_space$tree_size[i])

    # generate a set of unique tip labels that conform to standard
    tip_labels <- expand.grid(letters, letters, letters)
    tip_labels <- do.call(paste0, tip_labels)
    tip_labels <- tip_labels[1:parameter_space$tree_size[i]]
    tip_labels <- paste("bird", tip_labels, sep = "_")
    phylo$tip.label <- tip_labels

    prob_not_present <- 1 - parameter_space$prob_on_island[i]
    prob_endemic <-
      parameter_space$prob_endemic[i] * parameter_space$prob_on_island[i]
    prob_nonendemic <-
      (1 - parameter_space$prob_endemic[i]) * parameter_space$prob_on_island[i]


    empty_island <- TRUE
    while (empty_island) {
      # generate tip states under uniform sampling
      endemicity_status <- sample(
        x = c("endemic", "nonendemic", "not_present"),
        size = parameter_space$tree_size[i],
        replace = TRUE,
        prob = c(prob_endemic, prob_nonendemic, prob_not_present)
      )
      if (any(endemicity_status != "not_present")) {
        empty_island <- FALSE
      }
    }

    # add not present outgroup
    phylo <- DAISIEprep::add_outgroup(phylo)
    endemicity_status <- c("not_present", endemicity_status)

    # format data for DAISIEprep
    phylod <- phylobase::phylo4d(phylo, as.data.frame(endemicity_status))
    phylod <- DAISIEprep::add_asr_node_states(
      phylod = phylod,
      asr_method = "parsimony",
      tie_preference = "mainland",
      earliest_col = FALSE
    )


    # run extraction
    min_time <- microbenchmark::microbenchmark(
      island_tbl_min <- DAISIEprep::extract_island_species(
        phylod = phylod,
        extraction_method = "min",
        island_tbl = NULL,
        include_not_present = FALSE
      ),
      times = 10L
    )

    asr_time <- microbenchmark::microbenchmark(
      island_tbl_asr <- DAISIEprep::extract_island_species(
        phylod = phylod,
        extraction_method = "asr",
        island_tbl = NULL,
        include_not_present = FALSE
      ),
      times = 10L
    )

    median_times_min[j] <- median(min_time$time)
    median_times_asr[j] <- median(asr_time$time)
  }
  times_list[[i]] <- list(median_times_min, median_times_asr)
}

# convert list to data frame
results <- data.frame(
  parameter_space = rep(1:nrow(parameter_space), each = 2),
  tree_size = rep(parameter_space$tree_size, each = 2),
  prob_on_island = rep(parameter_space$prob_on_island, each = 2),
  prob_endemic = rep(parameter_space$prob_endemic, each = 2),
  extraction_method = rep(c("min", "asr"), nrow(parameter_space))
)

times <- unlist(lapply(times_list, function(x) {lapply(x, FUN =  median)}))

#convert from nanoseconds to seconds
times <- times / 1e9

results$median_time <- times

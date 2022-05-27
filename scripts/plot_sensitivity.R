# load all the files from the path
sensitivity_data_files <- list.files(
  path = system.file("sensitivity_data", package = "DAISIEprepExtra")
)
if (length(sensitivity_data_files) == 0) {
  stop("No results are in the results directory")
} else {
  file_paths <- as.list(paste0(
    system.file("sensitivity_data", package = "DAISIEprepExtra"),
    "/",
    sensitivity_data_files
  ))
  sensitivity_data <- lapply(file_paths, readRDS)
}

dna <- lapply(sensitivity_data, function(x) {
  lapply(x, "[[", 1)
})

complete <- lapply(sensitivity_data, function(x) {
  lapply(x, "[[", 2)
})

parameters <- lapply(sensitivity_data, function(x) {
  lapply(x, "[[", 3)
})

dna <- unlist(dna, recursive = FALSE)
dna_clado <- sapply(dna, "[[", "lambda_c")
dna_ext <- sapply(dna, "[[", "mu")
dna_k <- sapply(dna, "[[", "K")
dna_immig <- sapply(dna, "[[", "gamma")
dna_ana <- sapply(dna, "[[", "lambda_a")

complete <- unlist(complete, recursive = FALSE)
complete_clado <- sapply(complete, "[[", "lambda_c")
complete_ext <- sapply(complete, "[[", "mu")
complete_k <- sapply(complete, "[[", "K")
complete_immig <- sapply(complete, "[[", "gamma")
complete_ana <- sapply(complete, "[[", "lambda_a")

parameters <- unlist(parameters, recursive = FALSE)
extraction_method <- sapply(parameters, "[[", "extraction_method")
asr_method <- sapply(parameters, "[[", "asr_method")
tie_preference <- sapply(parameters, "[[", "tie_preference")

plotting_data_dna <- data.frame(
  extraction_method = extraction_method,
  asr_method = asr_method,
  tie_preference = tie_preference,
  dna_clado = dna_clado,
  dna_ext = dna_ext,
  #dna_k = dna_k,
  dna_immig = dna_immig,
  dna_ana = dna_ana
)

plotting_data_dna <- tidyr::unite(
  data = plotting_data_dna,
  col = extraction,
  extraction_method:tie_preference
)

plotting_data_dna <- tidyr::pivot_longer(
  data = plotting_data_dna,
  names_to = "parameter",
  dna_clado:dna_ana,
  values_to = "rates"
)

sensitivity <- ggplot2::ggplot(plotting_data_dna) +
  geom_density(
    mapping = ggplot2::aes(rates, fill = extraction),
    alpha = 0.5
  ) +
  ggplot2::scale_x_continuous(name = "Rate") +
  ggplot2::scale_y_continuous(name = "Density") +
  ggplot2::scale_fill_discrete(
    name = "Extraction Method",
    type = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"),
    labels = c("ASR (Mk island)", "ASR (Mk mainland)",
               "ASR (parsimony island)", "ASR (parsimony mainland)",
               "min")) +
  ggplot2::theme_classic() +
  facet_wrap(
    facets = "parameter",
    scales = "free",
    labeller = ggplot2::as_labeller(c(
      dna_ana = "Anagenesis",
      dna_clado = "Cladogenesis",
      dna_ext = "Extinction",
      dna_immig =  "Colonisation")
    )
  ) +
  ggplot2::theme(strip.background = ggplot2::element_blank(),
                 strip.text = ggtext::element_markdown())

ggplot2::ggsave(
  plot = sensitivity,
  filename = file.path("plots", "sensitivity.png"),
  device = "png",
  width = 300,
  height = 200,
  units = "mm",
  dpi = 600
)


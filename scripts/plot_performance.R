# load all the files from the path
performance_data_files <- list.files(
  path = system.file("performance_data", package = "DAISIEprepExtra")
)
if (length(sensitivity_data_files) == 0) {
  stop("No results are in the results directory")
} else {
  file_paths <- as.list(paste0(
    system.file("performance_data", package = "DAISIEprepExtra"),
    "/",
    performance_data_files
  ))
  performance_data <- lapply(file_paths, readRDS)
}

min <- lapply(performance_data, "[[", "min")

asr <- lapply(performance_data, "[[", "asr")

parameter <- lapply(performance_data, "[[", "parameter_index")

parameter_space <- expand.grid(
  tree_size = c(10, 50, 100, 500, 1000, 5000, 10000),
  prob_on_island = c(0.2, 0.5),
  prob_endemic = c(0.2, 0.8)
)

parameter <- lapply(parameter, function(x) {
  parameter_space[x, ]
})

tree_size <- sapply(parameter, "[[", "tree_size")
prob_on_island <- sapply(parameter, "[[", "prob_on_island")
prob_endemic <- sapply(parameter, "[[", "prob_endemic")

results <- data.frame(
  tree_size = rep(tree_size, each = 10),
  prob_on_island = rep(prob_on_island, each = 10),
  prob_endemic = rep(prob_endemic, each = 10),
  time_min = unlist(min),
  time_asr = unlist(asr)
)

results <- tidyr::pivot_longer(
  data = results,
  names_to = "parameter",
  time_min:time_asr,
  values_to = "rates"
)

grouped_results <- dplyr::group_by(results, tree_size)
mean_results <- dplyr::summarise(grouped_results, Mean = mean(rates))
sd_results <- dplyr::summarise(grouped_results, SD = sd(rates))

summary_results <- dplyr::right_join(mean_results, sd_results)

ggplot2::ggplot(data = summary_results) +
  ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = tree_size,
      y = Mean
      )
  ) +
  ggplot2::geom_errorbar(
    mapping = ggplot2::aes(
      x = tree_size,
      y = Mean,
      ymin = Mean-SD,
      ymax = Mean+SD
    ),
    width = 0.2
  ) +
  ggplot2::theme_classic() +
  ggplot2::scale_x_continuous(
    name = "Tree Size (num. tips)",
    trans = "log",
    breaks = scales::breaks_log()
  ) +
  ggplot2::scale_y_continuous(
    name = "Mean run time (seconds)",
    trans = "log",
    breaks = scales::breaks_log(),
    labels = scales::comma_format(accuracy = 0.1)
  )



#ggplot2::geom_point(
#    mapping = ggplot2::aes(
#      x = tree_size,
#      y = median_time,
#      colour = extraction_method,
#      shape = as.factor(prob_on_island),
#      size = as.factor(prob_endemic)),
#    alpha = 0.75) +

#ggplot2::guides(colour = ggplot2::guide_legend(title = "Extraction method"),
#                  shape = ggplot2::guide_legend(title = "Probability on island"),
#                  size = ggplot2::guide_legend(title = "Probability endemic")) +



ggplot2::ggsave(
  plot = performance,
  filename = file.path("plots", "performance.png"),
  device = "png",
  width = 150,
  height = 100,
  units = "mm",
  dpi = 600
)

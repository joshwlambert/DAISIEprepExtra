performance <- ggplot2::ggplot(data = results) +
  ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = tree_size,
      y = median_time,
      colour = extraction_method,
      shape = as.factor(prob_on_island),
      size = as.factor(prob_endemic)),
    alpha = 0.75) +
  ggplot2::theme_classic() +
  ggplot2::guides(colour = ggplot2::guide_legend(title = "Extraction method"),
                  shape = ggplot2::guide_legend(title = "Probability on island"),
                  size = ggplot2::guide_legend(title = "Probability endemic")) +
  ggplot2::scale_colour_brewer(palette="Set2") +
  ggplot2::scale_x_continuous(name = "Tree Size (num. tips)") +
  ggplot2::scale_y_continuous(name = "Median run time (seconds)")


ggplot2::ggsave(
  plot = performance,
  filename = file.path("plots", "performance.png"),
  device = "png",
  width = 150,
  height = 100,
  units = "mm",
  dpi = 600
)

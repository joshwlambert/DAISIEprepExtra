set.seed(
  8,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)
phylo <- ape::rphylo(n = 5, birth = 0.1, death = 0)
phylo$tip.label <- c("species_a", "species_b", "species_c", "species_d", "species_e")
phylo <- phylobase::phylo4(phylo)
endemicity_status <- c("endemic", "nonendemic", "endemic", "not_present",
                     "not_present")
phylod <- phylobase::phylo4d(phylo, as.data.frame(endemicity_status))
min_plot <- DAISIEprep::plot_phylod(phylod = phylod)
phylod <- DAISIEprep::add_asr_node_states(phylod = phylod, asr_method = "parsimony", tie_preference = "mainland")
asr_plot <- DAISIEprep::plot_phylod(phylod = phylod)

min_plot <- min_plot +
  ggplot2::ggtitle("min") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      hjust = 0.5,
      face = "bold",
      size = 16,
      family = "mono"
    )
  )
asr_plot <- asr_plot +
  ggplot2::ggtitle("asr") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      hjust = 0.5,
      face = "bold",
      size = 16,
      family = "mono"
    )
  )

min_plot <- min_plot +
  ggplot2::scale_colour_manual(
    values = c("#FF7F50", "#20B0B0", "#000000"),
    labels = c("Endemic", "Non-endemic", "Not present")
  )

asr_plot <- asr_plot +
  ggplot2::scale_colour_manual(
    values = c("#FF7F50", "#20B0B0", "#000000"),
    labels = c("Endemic", "Non-endemic", "Not present")
  )

min_plot <- min_plot +
  ggplot2::annotate(
    "text",
    x = -5,
    y = 3.0,
    label = "gamma",
    parse = TRUE,
    size = 5
  ) +
  ggplot2::annotate(
    "text",
    x = -2,
    y = 5.0,
    label = "gamma",
    parse = TRUE,
    size = 5
  ) +
  ggplot2::annotate(
    "text",
    x = -2,
    y = 4.0,
    label = "gamma",
    parse = TRUE,
    size = 5
  ) +
  ggplot2::annotate(
    "rect",
    xmin = -2,
    xmax = 0.5,
    ymin = 4.75,
    ymax = 5.25,
    alpha = 0.2,
    fill = "#FF7F50"
  ) +
  ggplot2::annotate(
    "rect",
    xmin = -2,
    xmax = 0.5,
    ymin = 4.25,
    ymax = 3.75,
    alpha = 0.2,
    fill = "#20B0B0"
  ) +
  ggplot2::annotate(
    "rect",
    xmin = -4.75,
    xmax = 0.5,
    ymin = 3.25,
    ymax = 2.75,
    alpha = 0.2,
    fill = "#FF7F50"
  )

asr_plot <- asr_plot +
  ggplot2::annotate(
    "text",
    x = -17.0,
    y = 4.0,
    label = "gamma",
    parse = TRUE,
    size = 5
  ) +
  ggplot2::annotate(
    "rect",
    xmin = -17.5,
    xmax = 0.5,
    ymin = 2.75,
    ymax = 5.25,
    alpha = 0.2,
    fill = "#FF7F50"
  )

prow <- cowplot::plot_grid(
  min_plot + ggplot2::theme(legend.position="none"),
  asr_plot + ggplot2::theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1
)


legend <- cowplot::get_legend(
  min_plot +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(x = 1, units = "cm"),
      legend.text = ggplot2::element_text(size = 12))
)

# add the legend underneath the row we made earlier. Give it 10%
# of the height of one plot (via rel_heights).
algo_plot <- cowplot::plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .2))

ggplot2::ggsave(
  plot = algo_plot,
  filename = file.path("plots", "algo_plot.png"),
  device = "png",
  width = 300,
  height = 150,
  units = "mm",
  dpi = 600
)

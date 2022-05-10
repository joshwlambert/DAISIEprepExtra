set.seed(
  1,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

phylo <- ape::rcoal(25)


phylo$tip.label <- c("species_a", "species_b", "species_c", "species_d",
                     "species_e", "species_f", "species_g", "species_h",
                     "species_i", "species_j", "species_k", "species_l",
                     "species_m", "species_n", "species_o", "species_p",
                     "species_q", "species_r", "species_s", "species_t",
                     "species_u", "species_v", "species_w", "species_x",
                     "species_y")

phylo <- phylobase::phylo4(phylo)

endemicity_status <- c("not_present", "not_present", "endemic",
                       "endemic", "endemic", "endemic",
                       "not_present", "not_present", "not_present",
                       "not_present", "not_present", "not_present",
                       "nonendemic", "nonendemic", "not_present",
                       "not_present", "not_present", "endemic",
                       "endemic", "endemic", "not_present",
                       "not_present", "not_present", "nonendemic",
                       "not_present")

phylod <- phylobase::phylo4d(phylo, as.data.frame(endemicity_status))

phylod <- DAISIEprep::add_asr_node_states(
  phylod = phylod,
  asr_method = "mk",
  tie_preference = "mainland"
)

  # generate plot
  # suppress Scale for 'y' is already present.
phylo <- ggtree::ggtree(phylod) +
  ggtree::theme_tree2() +
  ggtree::geom_tiplab(as_ylab = TRUE)

# suppress Scale for 'x' is already present.
phylo <- ggtree::revts(treeview = phylo) +
  ggplot2::scale_x_continuous(labels = abs) +
  ggplot2::xlab("Time (Million years ago)")

phylo <- phylo +
  ggtree::geom_tippoint(
    ggplot2::aes(colour = endemicity_status),
    size = 2.5,
    alpha = 0.5
  ) +
  ggplot2::labs(colour = "Endemicity status") +
  ggtree::geom_nodepoint(
    ggplot2::aes(colour = island_status),
    size = 2.5,
    alpha = 0.5
  ) +
  ggplot2::scale_colour_manual(
    values = c("#FF7F50", "#20B0B0", "#000000"),
    labels = c("Endemic", "Non-endemic", "Not present")
  )

island_tbl <- DAISIEprep::extract_island_species(
  phylod = phylod,
  extraction_method = "min"
)

cols <- DAISIEprep::plot_colonisation(
  island_tbl = island_tbl,
  island_age = 0.5
)

utility_plots <- cowplot::plot_grid(
  phylo,
  cols,
  nrow = 1,
  align = "h",
  labels = "AUTO",
  rel_widths = c(1, 0.75)
)

ggplot2::ggsave(
  plot = utility_plots,
  filename = file.path("plots", "utility_plots.png"),
  device = "png",
  width = 300,
  height = 100,
  units = "mm",
  dpi = 600
)

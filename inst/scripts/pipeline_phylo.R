set.seed(
  1,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)
phylo <- ape::rcoal(10)

phylo$tip.label <- c("bird_a", "bird_b", "bird_c", "bird_d", "bird_e", "bird_f",
                     "bird_g", "bird_h", "bird_i", "bird_j")

phylo <- phylobase::phylo4(phylo)

endemicity_status <- sample(
  x = c("not_present", "endemic", "nonendemic"),
  size = length(phylobase::tipLabels(phylo)),
  replace = TRUE,
  prob = c(0.6, 0.2, 0.2)
)

phylod <- phylobase::phylo4d(phylo, as.data.frame(endemicity_status))

phylod <- DAISIEprep::add_asr_node_states(phylod = phylod, asr_method = "parsimony")

phylo <- ggtree::ggtree(phylod) +
  ggtree::geom_tippoint(
    ggplot2::aes(colour = endemicity_status),
    size = 3,
    alpha = 0.75
  ) +
  ggplot2::labs(colour = "Island Species Endemicity") +
  ggtree::geom_nodepoint(
    ggplot2::aes(colour = island_status),
    size = 3,
    alpha = 0.75
  ) +
  ggplot2::scale_colour_manual(
    values = c("#FF7F50", "#20B0B0", "#000000"),
    labels = c("Endemic", "Non-endemic", "Not present")
  )

ggplot2::ggsave(
  plot = phylo,
  filename = file.path("plots", "pipeline_phylo.png"),
  device = "png",
  width = 150,
  height = 100,
  units = "mm",
  dpi = 600
)

island_tbl <- DAISIEprep::extract_island_species(
  phylod = phylod,
  extraction_method = "min"
)

island_tbl <- DAISIEprep::extract_island_species(
  phylod = phylod,
  extraction_method = "asr"
)


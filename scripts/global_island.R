library(ggplot2)
library(ggspatial)
library(rnaturalearth)
library(cowplot)

world <- ne_countries(scale = "large", returnclass = "sf")

world_map <- ggplot(data = world) +
  geom_sf(color = "gray50", fill = "gray50") +
  geom_rect(xmin = -166, xmax = -149, ymin = 16, ymax = 25,
            fill = NA, colour = "forestgreen", size = 1) +
  theme_bw() +
  theme(plot.margin = margin(0, 0, 0, 0))

hawaii_map <- ggplot(data = world) +
  geom_sf(color = "gray50", fill = "gray50") +
  coord_sf(xlim = c(-165, -150), ylim = c(18, 23), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "bl", which_north ="true",
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotate(geom = "text", x = -163, y = 23, label = "Hawaiian Archipelago",
           fontface = "bold", color = "grey22", size = 4) +
  annotate(geom = "text", x = -161, y = 19, label = "Pacific Ocean",
           fontface = "italic", color = "grey22", size = 4) +
  annotate(geom = "text", x = -153, y = 22, label = "Pacific Ocean",
           fontface = "italic", color = "grey22", size = 4) +
  ylab("Latitude") +
  xlab("Longitude") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "forestgreen", size = 2),
        plot.margin = margin(0, 0, 0, 0))


facet_map <- plot_grid(
  world_map,
  hawaii_map,
  nrow = 2,
  align = "v",
  axis = "lr",
  labels = c("A", "C")
)

library(TreeSim)
library(phylobase)
library(ggplot2)
library(ggtree)

set.seed(2)
phylo <- TreeSim::sim.bd.taxa(n = 25, numbsim = 1, lambda = 1, mu = 0)[[1]]
phylo <- phylo4(phylo)
tip_data <- data.frame(
  colour = c("black", "#008080", "black", "black", "#FF7F50", "black", "#FF7F50",
             "#FF7F50", "#FF7F50", "#FF7F50", "black", "black", "black", "black",
             "#FF7F50", "#FF7F50", "#FF7F50", "black", "black", "black", "black",
             "black", "black", "black", "black"),
  row.names = nodeId(phylo, "tip"))
node_data <- data.frame(
  colour = c("black", "black", "black", "black", "black", "#FF7F50", "black",
             "black", "#FF7F50", "black", "black", "black", "black", "black",
             "black", "black", "#FF7F50", "black", "#FF7F50", "#FF7F50",
             "black", "black", "black", "black"),
  row.names = nodeId(phylo, "internal"))

phylod <- phylo4d(x = phylo, tip.data = tip_data, node.data = node_data)

global_phylo <- ggtree(
  phylod, aes(
    color = I(colour)
  ),
  size = 0.6) +
  geom_tippoint(
    aes(colour = colour, shape = colour),
    alpha = 0.5,
    size = 1.5
  ) +
  geom_nodepoint(aes(colour = colour), alpha = 0.5, size = 1.5) +
  geom_highlight(
    data = data.frame(
      node = c(2, 17, 31, 44),
      type = c("Non-endemic", "Endemic", "Endemic", "Endemic")
    ),
    aes(node = node, fill = type),
    fill = c("#008080", "#FF7F50", "#FF7F50", "#FF7F50"),
    type = "roundrect",
    alpha = 0.3
  ) +
  ggplot2::scale_shape_manual(values = c(15, 16, 20)) +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    legend.position = "none"
  )


horizontal <- data.frame(
  species_id = as.factor(c(1, 2, 4, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)),
  species_type = c("Non-endemic", rep("Endemic", 13)),
  col_t = c(0.29, 1.37, 0.54, NA, NA, NA, 1.85, NA, NA, NA, NA, NA, NA, NA),
  spec_origin_t = c(0.29, 1.37, 0.54, 0.21, 0.21, 1.16, 1.85, 0.76, 1.16, 0.30,
                    0.76, 0.19, 0.3, 0.19),
  spec_ex_t = c(0, 0, 0.21, 0, 0, 0, 1.16, 0, 0.76, 0, 0.3, 0, 0.19, 0)
)

vertical <- data.frame(
  species_type = c(rep("Endemic", 5)),
  x = c(0.21, 1.16, 0.76, 0.3, 0.19),
  y = as.factor(c(3, 6, 8, 10, 12)),
  yend = as.factor(c(5, 9, 11, 13, 14))
)

violin <- data.frame(
  col_t = c(
    rnorm(n = 1000, mean = 0.29, sd = 0.01),
    rnorm(n = 1000, mean = 1.37, sd = 0.02),
    rnorm(n = 1000, mean = 0.54, sd = 0.02),
    rnorm(n = 1000, mean = 1.85, sd = 0.01)
  ),
  species_id = c(
    rep(1, 1000),
    rep(2, 1000),
    rep(4, 1000),
    rep(7, 1000)
  ),
  species_type = c(
    rep("Non-endemic", 1000),
    rep("Endemic", 1000),
    rep("Endemic", 1000),
    rep("Endemic", 1000)
  )
)

island_phylo <- ggplot2::ggplot(data = horizontal) +
  ggplot2::geom_segment(
    data = horizontal,
    ggplot2::aes(
      x = spec_origin_t,
      xend = spec_ex_t,
      y = species_id,
      yend = species_id,
      colour = species_type
    ),
  ) +
  ggplot2::geom_segment(
    data = vertical,
    ggplot2::aes(
      x = x,
      xend = x,
      y = y,
      yend = yend,
      colour = species_type
    )
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = col_t,
      y = species_id,
      color = species_type,
      shape = species_type
    ), size = 2) +
  ggplot2::geom_violin(
    data = violin,
    ggplot2::aes(x = col_t, y = as.factor(species_id), colour = species_type),
    adjust = 3,
    alpha = 0.1
  ) +
  ggplot2::scale_shape_manual(values = c(16, 15)) +
  ggplot2::scale_colour_manual(values = c("#FF7F50", "#008080")) +
  ggplot2::scale_x_continuous(name = "Time before present (Myr)",
                              trans = "reverse", limits = c(2, 0)) +
  ggplot2::guides(colour = guide_legend(title = "Island Species Endemicity"),
                  shape = guide_legend(title = "Island Species Endemicity")) +
  ggplot2::geom_vline(xintercept = 2, linetype = "dashed") +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.ticks.y = ggplot2::element_blank(),
    axis.line.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

prow <- plot_grid(
  global_phylo + theme(legend.position="none"),
  island_phylo + theme(legend.position="none"),
  align = 'vh',
  labels = c("B", "D"),
  hjust = -1,
  nrow = 2,
  rel_heights = c(1, 1)
)
legend <- get_legend(
  # create some space to the left of the legend
  island_phylo + theme(legend.box.margin = margin(0, 0, 0, 12))
)
facet_phylo <- plot_grid(prow, legend, rel_widths = c(3, .8))

global_island <- cowplot::plot_grid(facet_map, facet_phylo, nrow = 1)

ggplot2::ggsave(
  plot = global_island,
  filename = file.path("plots", "global_island.png"),
  device = "png",
  width = 300,
  height = 180,
  units = "mm",
  dpi = 600
)



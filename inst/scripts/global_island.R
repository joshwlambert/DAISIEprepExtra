library(ggplot2)
library(ggspatial)
library(rnaturalearth)
library(cowplot)
library(sf)
library(TreeSim)
library(phylobase)
library(ggtree)

load("data/world_border.rda")
world <- ne_countries(scale = "medium", returnclass = "sf")

PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

world <- st_transform(world, PROJ)

world_map <- ggplot(data = world) +
  geom_sf(color = "gray50", fill = "gray50") +
  geom_polygon(
    data=world_border,
    aes(x=long, y=lat),
    colour="black",
    fill="transparent",
    size = 0.25
  ) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0),
        panel.grid.major = element_line(colour = "transparent"))

world <- ne_countries(scale = "large", returnclass = "sf")

hawaii_map <- ggplot(data = world) +
  geom_sf(color = "gray50", fill = "gray50") +
  coord_sf(xlim = c(-162.5, -152.4), ylim = c(18.5, 22.5), expand = TRUE) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0),
        panel.border = ggplot2::element_rect(
          fill = NA,
          colour = "black",
          size = 2
        )
  )

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

global <- cowplot::ggdraw(global_phylo) +
  cowplot::draw_plot(world_map, 0, 0.6, 0.4, 0.4) +
  cowplot::draw_text(text = "1", x = 0.5, y = 0.28, size = 10) +
  cowplot::draw_text(text = "2", x = 0.74, y = 0.53, size = 10) +
  cowplot::draw_text(text = "3", x = 0.88, y = 0.72, size = 10) +
  cowplot::draw_text(text = "4", x = 0.86, y = 0.90, size = 10)

horizontal <- data.frame(
  species_id = as.factor(c(1, 2, 4, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)),
  species_type = c(rep("Endemic", 10), "Non-endemic", rep("Endemic", 3)),
  col_t = c(NA, 1.85, NA, NA, NA, NA, NA, NA, NA, 1.37, 0.29, NA, 0.54, NA),
  spec_origin_t = c(1.16, 1.85, 1.16, 0.76, 0.30, 0.76, 0.19, 0.3, 0.19, 1.37,
                    0.29, 0.21, 0.54, 0.21),
  spec_ex_t = c(0, 1.16, 0.76, 0, 0, 0.3, 0, 0.19, 0, 0, 0, 0, 0.21, 0)
)

vertical <- data.frame(
  species_type = c(rep("Endemic", 5)),
  x = c(0.21, 1.16, 0.76, 0.3, 0.19),
  y = as.factor(c(12, 1, 3, 5, 7)),
  yend = as.factor(c(14, 4, 6, 8, 9))
)

violin <- data.frame(
  col_t = c(
    rnorm(n = 1000, mean = 1.85, sd = 0.016),
    rnorm(n = 1000, mean = 1.37, sd = 0.02),
    rnorm(n = 1000, mean = 0.29, sd = 0.02),
    rnorm(n = 1000, mean = 0.54, sd = 0.014)
  ),
  species_id = c(
    rep(2, 1000),
    rep(10, 1000),
    rep(11, 1000),
    rep(13, 1000)
  ),
  species_type = c(
    rep("Endemic", 1000),
    rep("Endemic", 1000),
    rep("Non-endemic", 1000),
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
    alpha = 0.1,
    width = 1.5
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
  ) +
  cowplot::draw_text(text = "1", x = 1.5, y = 2.5, size = 10) +
  cowplot::draw_text(text = "2", x = 1.2, y = 10.5, size = 10) +
  cowplot::draw_text(text = "3", x = 0.4, y = 11.5, size = 10) +
  cowplot::draw_text(text = "4", x = 0.45, y = 13.5, size = 10)

legend <- get_legend(island_phylo)

island_phylo <- island_phylo + theme(legend.position = "none")

island <- cowplot::ggdraw(island_phylo) +
  cowplot::draw_plot(hawaii_map, 0.1, 0.65, 0.4, 0.4)

prow <- plot_grid(
  global + theme(legend.position="none"),
  island + theme(legend.position="none"),
  align = 'vh',
  labels = "AUTO",
  label_x = -0.01,
  hjust = -1,
  nrow = 1,
  rel_heights = c(1, 1)
)

global_island <- plot_grid(prow, legend, rel_widths = c(4, 1))

ggplot2::ggsave(
  plot = global_island,
  filename = file.path("inst", "plots", "global_island.png"),
  device = "png",
  width = 300,
  height = 150,
  units = "mm",
  dpi = 600
)


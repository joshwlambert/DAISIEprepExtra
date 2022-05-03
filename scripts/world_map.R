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
  theme(panel.background = element_rect(fill = "aliceblue"))

hawaii_map <- ggplot(data = world) +
  geom_sf(color = "gray50", fill = "gray50") +
  coord_sf(xlim = c(-165, -150), ylim = c(18, 23), expand = TRUE) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "bl", which_north ="true",
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotate(geom = "text", x = -163, y = 23, label = "Hawaiian Archipelago",
           fontface = "bold", color = "grey22", size = 5) +
  annotate(geom = "text", x = -161, y = 19, label = "Pacific Ocean",
           fontface = "italic", color = "grey22", size = 5) +
  annotate(geom = "text", x = -153, y = 22, label = "Pacific Ocean",
           fontface = "italic", color = "grey22", size = 5) +
  ylab("Latitude") +
  xlab("Longitude") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "aliceblue"),
        panel.border = element_rect(colour = "forestgreen", size = 2))


facet_map <- plot_grid(world_map, hawaii_map, nrow = 2, align = "v", axis = "lr")

ggsave(facet_map)



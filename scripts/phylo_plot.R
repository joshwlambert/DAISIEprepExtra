library(TreeSim)
library(phylobase)
library(ggplot2)
library(ggtree)

set.seed(2)
phylo <- TreeSim::sim.bd.taxa(n = 25, numbsim = 1, lambda = 1, mu = 0)[[1]]
phylo <- phylo4(phylo)
tip_data <- data.frame(
  colour = c("black", "blue", "black", "black", "red", "black", "red",
             "red", "red", "red", "black", "black", "black", "black",
             "red", "red", "red", "black", "black", "black", "black",
             "black", "black", "black", "black"),
  row.names = nodeId(phylo, "tip"))
node_data <- data.frame(
  colour = c("black", "black", "black", "black", "black", "red", "black",
             "black", "red", "black", "black", "black", "black", "black",
             "black", "black", "red", "black", "red", "red",
             "black", "black", "black", "black"),
  row.names = nodeId(phylo, "internal"))

phylod <- phylo4d(x = phylo, tip.data = tip_data, node.data = node_data)

ggtree(phylod, aes(color = I(colour)), layout = "circular", size = 0.6) +
  geom_tippoint(aes(colour = colour), alpha = 0.5, size = 2, stroke = 1) +
  geom_nodepoint(aes(colour = colour), alpha = 0.5, size = 2, stroke = 1) +
  geom_highlight(
    data = data.frame(
      node = c(2, 17, 31, 44),
      type = c("Non-endemic", "Endemic", "Endemic", "Endemic")
    ),
    aes(node = node, fill = type),
    type = "rect"
  ) +
  labs(fill =  "Island Species Endemicity")


island <- data.frame(
  species_id = as.factor(c(rep(1, 10), 2, 3, 4, 5, 6, 7, 8, 9, 10)),
  clade_id = c(rep(1, 10), 2, 2, 2, 2, 2, 3, 3, 3, 4),
  species_type = c(rep("endemic", 10), "endemic", "endemic", "endemic", "endemic",
                   "endemic", "endemic", "endemic", "endemic", "nonendemic"),
  mean_col_t = c(rep(1.26, 10), 0.07, NA, NA, NA, NA, 1.57, NA, NA, 0.56),
  spec_origin_t = c(rep(1.26, 10), 0.07, 1.04, 1.04, 1.62, 1.62, 1.57, 1.75, 1.75, 0.56),
  spec_ex_t = c(rep(2, 10), 1.03, 2, 1.62, 2, 2, 1.75, 2, 2, 2),
  branch_code = c(rep("A", 10), "A", "AA", "AB", "ABA", "ABB", "A", "AA", "AB", "A"),
  col_t = c(rnorm(n = 10, mean = 1.26, sd = 0.2), 0.07, NA, NA, NA, NA, 1.57, NA, NA, 0.56)
)

# x1 = x = = spec_origin_t
# x2 = xend = spec_ex_t (always conveniently stopped at time 'total_time')
# y1 = y = unique_species_id                                                  # nolint this is no commented code
# y2 = yend = unique_species_id                                               # nolint this is no commented code
island$y <- Vectorize(DAISIEmainland::branch_code_to_y)(
  island$branch_code
)

# Number all species of all clades individually
# Keep the clade ID first; it is assumed the branch code is at the end:
# by removing the last character, the ancestor is found
island$unique_species_id <- paste0(
  island$clade_id, "-", island$branch_code
)

# Do not make a factor, as we need to work on the string
# t_mainland$unique_species_id <- as.factor(t_mainland$unique_species_id)     # nolint yup, this is code

# Create a table for the vertical lines
# x1 = spec_ex_t (of ancestor)
# x2 = spec_origin_t (of derived)
# y1 = y of branch_code ancestor
# y2 = y of branch_code of derived species
#
# Work backwards
t_ancestors <- data.frame(
  ancestor_branch_code = island$branch_code,
  ancestor_unique_species_id = island$unique_species_id,
  ancestor_spec_ex_t = island$spec_ex_t,
  clade_id = island$clade_id
)
t_offspring <- data.frame(
  offspring_branch_code = island$branch_code,
  offspring_unique_species_id = island$unique_species_id,
  ancestor_unique_species_id = strtrim(
    island$unique_species_id,
    nchar(island$unique_species_id) - 1
  ),
  offspring_spec_origin_t = island$spec_origin_t
)
t_vertical <- merge(t_ancestors, t_offspring)
t_vertical$ancestor_y <- Vectorize(DAISIEmainland::branch_code_to_y)(
  t_vertical$ancestor_branch_code
)
t_vertical$offspring_y <- Vectorize(DAISIEmainland::branch_code_to_y)(
  t_vertical$offspring_branch_code
)

# Here, we reverse the time axis,
# from time after the island came into existance,
# to time before present

# aka the island age. There is always a species at the present
total_time <- max(island$spec_ex_t)
island$spec_origin_t <- total_time - island$spec_origin_t
island$spec_ex_t <- total_time - island$spec_ex_t
island$mean_col_t <- total_time - island$mean_col_t
t_vertical$ancestor_spec_ex_t <- total_time - t_vertical$ancestor_spec_ex_t
t_vertical$offspring_spec_origin_t <- total_time - t_vertical$offspring_spec_origin_t # nolint indeed a long line
t_vertical$species_type <- rep("endemic", nrow(t_vertical))

p <- ggplot2::ggplot(data = island) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = spec_origin_t,
      xend = spec_ex_t,
      y = y,
      yend = y,
      colour = species_type
    )
  ) + ggplot2::geom_segment(
    data = t_vertical,
    ggplot2::aes(
      x = ancestor_spec_ex_t,
      xend = offspring_spec_origin_t,
      y = ancestor_y,
      yend = offspring_y,
      colour = species_type
    )
  ) +
  geom_point(
    mapping = aes(
      x = mean_col_t,
      y = y,
      shape = species_type,
      color = species_type),
    size = 3
  ) +
  scale_shape_manual(values = c(16, 15))

p <- p + ggplot2::scale_x_reverse(
  name = "Time before present (Myr)",
  limits = c(2, 0)
) +
  ggplot2::facet_grid(
    clade_id ~ .,
    scales = "free",
    space = "free"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.line.y = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_blank(),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 2, linetype = "dashed")

ggplot2::ggplot(data = island, aes(species_id, col_t)) +
  geom_violin()
p <- ggplot(mtcars, aes(factor(cyl), mpg))
p + geom_violin()

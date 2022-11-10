library(ggplot2)
library(cowplot)
library(patchwork)
library(ggh4x)

source("scripts/sensitivity/functions.R")

theme_set(theme_cowplot() + background_grid())

res_dir_base <- "res_special"
plot_dir_base <- paste(res_dir_base, "plots", sep = "/")

if (!dir.exists(plot_dir_base)) {
  dir.create(plot_dir_base, recursive = TRUE)
}

params <- expand.grid(
  c(2, 1000),
  c(2, 1000),
  c(0, 25),
  c(0, 1),
  -5:5 %>% as.numeric
) %>%
  as.data.table %>%
  setnames(c("j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag", "gamma", "pe"))

# we can luckily use this function without it overwriting our specified pe
# values, as within the function it does params[, pe := pe], which just sets
# them to themselves, if they have already been assigned (if they haven't, the
# pe on the RHS references the global variable)
add_constants_and_gammas_to_param_table(params)

read_spatial_data_at_times <- function(times) {
  data <- lapply(
    times, function(t_inf) read_spatial_data(params, t_inf, res_dir_base)
  )

  data <- rbindlist(data)
  # Round the time to a few decimal places as a hack to ensure that all times
  # we consider "equal" when plotting (e.g. 30.000000 and 30.000020) are
  # treated as such by facet_grid etc.
  data[, time_inf := round(time_inf, 3)]
  melt(data, measure.vars = c("C_b", "phi_C_b"))
}

# Combined spatial plots
# ----------------------

spatial_plot_subset_combined <- function(plot_times, pes, lag) {
  data_subset_long <- read_spatial_data_at_times(plot_times)

  # Manually set the factor labels so they appear correctly in the facet strips
  # without having to use a complex labeller
  data_subset_long[, j_phi_i_i_factor := factor(
    j_phi_i_i_factor,
    labels = c(
      expression(paste(J[phi[i]]^I, " factor = 2")),
      expression(paste(J[phi[i]]^I, " factor = 1000"))
    )
  )]

  data_subset_long[, m_i_factor := factor(
    m_i_factor,
    labels = c(
      expression(paste(M^I, " factor = 2")),
      expression(paste(M^I, " factor = 1000"))
    )
  )]

  data_subset_long[, time_inf := factor(
    time_inf,
    labels = sapply(plot_times, function(t) paste0("t[inf] == ", t / 10))
  )]

  data_subset_long[, variable := factor(
    variable,
    labels = c(
      variable_labels["C_b"],
      variable_labels["phi_C_b"]
    )
  )]

  ggplot(
    data_subset_long[t_j_phi_i_lag == lag & pe %in% pes],
    aes(
      x = x,
      y = value,
      group = rep,
      colour = factor(pe),
      linetype = factor(gamma)
    )
  ) +
    facet_nested(
      rows = vars(j_phi_i_i_factor, time_inf),
      cols = vars(m_i_factor, variable),
      scales = "free_y",
      independent = "y",
      labeller = label_parsed,
      switch = "y"
    ) +
    geom_line(size = 1) +
    scale_linetype_discrete(limits = rev) +
    theme(
      strip.placement = "outside",
      strip.text = element_text(size = rel(1.25)),
      strip.background = element_rect(fill = "#eaeaea")
    ) +
    labs(
      x = expression(x),
      y = NULL,
      colour = param_labels["pe"],
      linetype = param_labels["gamma"]
    )
}

spatial_plot_subset <- function(plot_time, pes, lag) {
  data_subset_long <- read_spatial_data_at_times(plot_time)

  # Manually set the factor labels so they appear correctly in the facet strips
  # without having to use a complex labeller
  data_subset_long[, j_phi_i_i_factor := factor(
    j_phi_i_i_factor,
    labels = c(
      expression(paste(J[phi[i]]^I, " factor = 2")),
      expression(paste(J[phi[i]]^I, " factor = 1000"))
    )
  )]

  data_subset_long[, m_i_factor := factor(
    m_i_factor,
    labels = c(
      expression(paste(M^I, " factor = 2")),
      expression(paste(M^I, " factor = 1000"))
    )
  )]

  data_subset_long[, variable := factor(
    variable,
    labels = c(
      variable_labels["C_b"],
      variable_labels["phi_C_b"]
    )
  )]

  ggplot(
    data_subset_long[
      pe %in% pes &
        # NOTE: dividing the output number by 10 to get the time is a bit of a
        # hack
        time_inf == plot_time / 10 &
        t_j_phi_i_lag == lag
    ],
    aes(
      x = x,
      y = value,
      colour = factor(pe),
      linetype = factor(gamma),
      group = rep
    )
  ) +
    geom_line(size = 1) +
    facet_nested(
      rows = vars(j_phi_i_i_factor),
      cols = , vars(m_i_factor, variable),
      scales = "free_y",
      independent = "y",
      labeller = label_parsed,
      switch = "y"
    ) +
    scale_linetype_discrete(limits = rev) +
    theme(
      strip.placement = "outside",
      #strip.text = element_text(size = rel(1.25)),
      strip.background = element_rect(fill = "#eaeaea")
    ) +
    labs(
      x = "Spatial distance",
      #x = expression(x),
      y = NULL,
      colour = param_labels["pe"],
      linetype = param_labels["gamma"]
    )
}

spatial_plot_homeostasis <- function(pes) {
  data_subset_hom <- read_spatial_data_at_times(0)
  data_subset_hom <- data_subset_hom[
    pe %in% pes &
      m_i_factor == 2 &
      j_phi_i_i_factor == 2 &
      t_j_phi_i_lag == 0
  ]

  data_subset_hom[, variable := factor(
    variable,
    labels = c(
      variable_labels["C_b"],
      variable_labels["phi_C_b"]
    )
  )]

  ggplot(
    data_subset_hom,
    aes(
      x = x,
      y = value,
      colour = factor(pe),
      linetype = factor(gamma)
    )
  ) +
    geom_line(size = 1) +
    facet_wrap(
      vars(variable),
      scales = "free_y",
      labeller = label_parsed
    ) +
    scale_linetype_discrete(limits = rev) +
    theme(
      strip.placement = "outside",
      #strip.text = element_text(size = rel(1.25)),
      strip.background = element_rect(fill = "#eaeaea")
    ) +
    labs(
      x = "Spatial distance",
      #x = expression(x),
      y = NULL,
      colour = param_labels["pe"],
      linetype = param_labels["gamma"]
    )
}

pes <- c(-5, -3, -1, 1, 3, 5)

# Spatial plots for a range of time points

# No lag
plot_times <- c(0, 50, 100, 150, 250, 500)
p_spatial_lag_0 <- spatial_plot_subset_combined(plot_times, pes, 0)

ggsave_with_defaults(
  plot = p_spatial_lag_0,
  paste(plot_dir_base, "spatial_lag_0.pdf", sep = "/"),
  width = 14,
  height = 20
)

# Lag 25
plot_times <- c(0, 250, 300, 350, 400, 500)
p_spatial_lag_25 <- spatial_plot_subset_combined(plot_times, pes, 25)

ggsave_with_defaults(
  plot = p_spatial_lag_25,
  paste(plot_dir_base, "spatial_lag_25.pdf", sep = "/"),
  width = 14,
  height = 20
)

#########

p_spatial_hom <- spatial_plot_homeostasis(pes)

p_spatial_hom +
  #guides(colour = guide_legend(nrow = 1)) +
  theme(
    legend.position = "bottom",
    legend.justification = "centre",
    legend.box.just = "top",
    legend.box = "horizontal"
  )

p_spatial_t_50_lag_0 <- spatial_plot_subset(500, pes, 0)

# Aligns well, collects guides, but unequal panel sizes
(p_spatial_hom / p_spatial_t_50_lag_0) +
  plot_layout(guides = "collect")

# Aligns well, collects guides, equal panel sizes, but completely broken axes
# (because they are aligned before the panels are resized)
(p_spatial_hom / p_spatial_t_50_lag_0) +
  plot_layout(guides = "collect") &
  force_panelsizes(rows = unit("0.50", "npc"), cols = unit("0.21", "npc"))

# Force the sizes of the panels to a fixed value, extract the shared legend to
# a separate grob thing, and store each plot without the legends (but with a
# title)
forced_sizes <- force_panelsizes(rows = unit("0.25", "npc"), cols = unit("0.2", "npc"))
leg <- get_legend(p_spatial_hom)
p1 <- p_spatial_hom +
  forced_sizes +
  theme(legend.position = "none") +
  labs(title = expression(t[inf] == 0), tag = "A")
p2 <- p_spatial_t_50_lag_0 +
  forced_sizes +
  theme(legend.position = "none") +
  labs(title = expression(t[inf] == 50), tag = "B")

# Save the two plots without legends to pdf
ggsave_with_defaults("p1.pdf", plot = p1, width = 12)
ggsave_with_defaults("p2.pdf", plot = p2, width = 12)

# Save the shared legend to pdf
ggsave_with_defaults("legend.pdf", plot = leg, width = 12)

# Puts the plots one above the other, retaining the fixed sizes, but leaves too
# much of a gap. Presumably it's assigning 1/2 of the height to each plot, but
# the top one is smaller so it pads it. If you specify say `rel_widths = c(1,
# 2)`, it resizes the panels. Axis alignment stuff doesn't do much (anything?)
# in this case
plot_grid(p1, p2, ncol = 1, axis = "ltrb", align = "v")

# Essentially the same result as `plot_grid`
grid.arrange(p1, p2, ncol = 1)

library(egg)

# Very similar to basic patchwork, but (afaik) can't collect guides
ggarrange(p_spatial_hom, p_spatial_t_50_lag_0)

fixed_size_panels <- lapply(
  list(
    p_spatial_hom + theme(legend.position = "none"),
    p_spatial_t_50_lag_0 + theme(legend.position = "none")
  ),
  set_panel_size,
  width = unit(0.19, "npc"),
  height = unit(0.38, "npc")
)

grid.arrange(grobs = fixed_size_panels)

plot_grid(
  fixed_size_panels[[1]],
  fixed_size_panels[[2]],
  ncol = 1
)

plot_grid(
  plot_grid(
    fixed_size_panels[[1]],
    fixed_size_panels[[2]],
    ncol = 1
  ),
  leg,
  ncol = 2,
  rel_widths = c(20, 1)
)

(
  wrap_elements(fixed_size_panels[[1]]) +
  wrap_elements(fixed_size_panels[[2]]) +
  plot_layout(ncol = 1, byrow = FALSE)
) + wrap_elements(leg) +
plot_layout(ncol = 2, widths = c(20, 1))
#+
  #plot_layout(widths = , guides = "collect")

ggsave_with_defaults(
  plot = p_asdf,
  "asdf.pdf"
)

################

a <- plot_grid(
  NULL,
  p_spatial_hom + theme(legend.position = "none"),
  NULL,
  axis = "lr",
  align = "v",
  ncol = 3,
  rel_widths = c(1, 2, 1)
)
plot_grid(a,
  p_spatial_t_50_lag_0 + theme(legend.position = "none"),
  axis = "lr",
  align = "v",
  ncol = 1,
  rel_heights = c(1, 2)
)

(plot_spacer() + p_spatial_hom + plot_spacer()) / p_spatial_t_50_lag_0 + plot_layout(guides = "collect")

(plot_spacer() + p_spatial_hom) /
  p_spatial_t_50_lag_0 +
  plot_layout(ncol = 1, nrow = 2, guides = "collect")

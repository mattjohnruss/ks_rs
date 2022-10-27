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

pes <- c(-5, -3, -1, 1, 3, 5)

plot_times <- c(0, 50, 100, 150, 250, 500)
p_spatial_lag_0 <- spatial_plot_subset_combined(plot_times, pes, 0)

ggsave_with_defaults(
  plot = p_spatial_lag_0,
  paste(plot_dir_base, "spatial_lag_0.pdf", sep = "/"),
  width = 14,
  height = 20
)

plot_times <- c(0, 250, 300, 350, 400, 500)
p_spatial_lag_25 <- spatial_plot_subset_combined(plot_times, pes, 25)

ggsave_with_defaults(
  plot = p_spatial_lag_25,
  paste(plot_dir_base, "spatial_lag_25.pdf", sep = "/"),
  width = 14,
  height = 20
)

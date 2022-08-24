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

spatial_plot_subset <- function(times,
                                params_subset,
                                colour_by,
                                linetype_by,
                                alpha_by) {
  data <- lapply(
    times, function(t_inf) read_spatial_data(params, t_inf, res_dir_base)
  )

  data <- rbindlist(data)
  # Round the time to a few decimal places as a hack to ensure that all times
  # we consider "equal" when plotting (e.g. 30.000000 and 30.000020) are
  # treated as such by facet_grid etc.
  data[, time_inf := round(time_inf, 3)]
  data_long <- melt(data, measure.vars = c("C_b", "phi_C_b"))

  p <- ggplot(
    data_long[
      j_phi_i_i_factor %in% as.numeric(params_subset$j_phi_i_i_factor) &
        m_i_factor %in% as.numeric(params_subset$m_i_factor) &
        t_j_phi_i_lag %in% as.numeric(params_subset$t_j_phi_i_lag) &
        gamma %in% as.numeric(params_subset$gamma) &
        pe %in% as.numeric(params_subset$pe)
    ],
    aes(
      x = x,
      y = value,
      group = rep,
      colour = factor({{ colour_by }})
    )
  ) +
  facet_nested(
    rows = vars(time_inf),
    cols = vars(variable),
    scales = "free",
    independent = "y",
    labeller = labeller(
      .cols = as_labeller(function(v) all_labels[v], label_parsed),
      .rows = as_labeller(
        function(t) paste0("t[inf] == ", round(as.numeric(t), 3)),
        label_parsed
      )
    ),
    switch = "y"
  )

  p <- p +
    geom_line(size = 1) +
    theme(
      #strip.background.y = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = rel(1.25)),
      #strip.text.y.left = element_text(angle = 0)
    ) +
    labs(
      x = expression(x),
      y = variable_labels[deparse(substitute(variable))],
      colour = param_labels[deparse(substitute(colour_by))]
    )

  if (!missing(linetype_by)) {
    p <- p + aes(linetype = factor({{ linetype_by }})) +
      labs(linetype = param_labels[deparse(substitute(linetype_by))])
  }

  if (!missing(alpha_by)) {
    # scale_alpha_manual only works for a parameter with at most the number of
    # values specified here
    p <- p + aes(alpha = factor({{ alpha_by }})) +
      scale_alpha_manual(values = c(1.0, 0.4)) +
      labs(alpha = param_labels[deparse(substitute(alpha_by))])
  }

  p
}

plot_times <- c(0, 50, 100, 150, 250, 500)

params_subset <- list(
  gamma = c(0, 1),
  pe = c(-5, -3, -1, 1, 3, 5)
)

# t_j_phi_i_lag = 0
# -----------------
plot_dir <- paste(plot_dir_base, "t_j_phi_i_lag=0", sep = "/")

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

params_subset$t_j_phi_i_lag <- 0

# j_phi_i_i_factor = 2, m_i_factor = 2
params_subset$j_phi_i_i_factor <- 2
params_subset$m_i_factor <- 2

p_spatial_2_2 <- spatial_plot_subset(
  plot_times,
  params_subset,
  colour_by = pe,
  linetype_by = gamma
)

ggsave_with_defaults(
  plot = p_spatial_2_2,
  paste(plot_dir, "spatial_2_2.pdf", sep = "/"),
  width = 8.5,
  height = 14
)

# j_phi_i_i_factor = 1000, m_i_factor = 2
params_subset$j_phi_i_i_factor <- 1000
params_subset$m_i_factor <- 2

p_spatial_1000_2 <- spatial_plot_subset(
  plot_times,
  params_subset,
  colour_by = pe,
  linetype_by = gamma
)

ggsave_with_defaults(
  plot = p_spatial_1000_2,
  paste(plot_dir, "spatial_1000_2.pdf", sep = "/"),
  width = 8.5,
  height = 14
)

# j_phi_i_i_factor = 2, m_i_factor = 1000
params_subset$j_phi_i_i_factor <- 2
params_subset$m_i_factor <- 1000

p_spatial_2_1000 <- spatial_plot_subset(
  plot_times,
  params_subset,
  colour_by = pe,
  linetype_by = gamma
)

ggsave_with_defaults(
  plot = p_spatial_2_1000,
  paste(plot_dir, "spatial_2_1000.pdf", sep = "/"),
  width = 8.5,
  height = 14
)

# j_phi_i_i_factor = 1000, m_i_factor = 1000
params_subset$j_phi_i_i_factor <- 1000
params_subset$m_i_factor <- 1000

p_spatial_1000_1000 <- spatial_plot_subset(
  plot_times,
  params_subset,
  colour_by = pe,
  linetype_by = gamma
)

ggsave_with_defaults(
  plot = p_spatial_1000_1000,
  paste(plot_dir, "spatial_1000_1000.pdf", sep = "/"),
  width = 8.5,
  height = 14
)

# t_j_phi_i_lag = 25
# -----------------
plot_dir <- paste(plot_dir_base, "t_j_phi_i_lag=25", sep = "/")

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

params_subset$t_j_phi_i_lag <- 25

# j_phi_i_i_factor = 2, m_i_factor = 2
params_subset$j_phi_i_i_factor <- 2
params_subset$m_i_factor <- 2

p_spatial_2_2 <- spatial_plot_subset(
  plot_times,
  params_subset,
  colour_by = pe,
  linetype_by = gamma
)

ggsave_with_defaults(
  plot = p_spatial_2_2,
  paste(plot_dir, "spatial_2_2.pdf", sep = "/"),
  width = 8.5,
  height = 14
)

# j_phi_i_i_factor = 1000, m_i_factor = 2
params_subset$j_phi_i_i_factor <- 1000
params_subset$m_i_factor <- 2

p_spatial_1000_2 <- spatial_plot_subset(
  plot_times,
  params_subset,
  colour_by = pe,
  linetype_by = gamma
)

ggsave_with_defaults(
  plot = p_spatial_1000_2,
  paste(plot_dir, "spatial_1000_2.pdf", sep = "/"),
  width = 8.5,
  height = 14
)

# j_phi_i_i_factor = 2, m_i_factor = 1000
params_subset$j_phi_i_i_factor <- 2
params_subset$m_i_factor <- 1000

p_spatial_2_1000 <- spatial_plot_subset(
  plot_times,
  params_subset,
  colour_by = pe,
  linetype_by = gamma
)

ggsave_with_defaults(
  plot = p_spatial_2_1000,
  paste(plot_dir, "spatial_2_1000.pdf", sep = "/"),
  width = 8.5,
  height = 14
)

# j_phi_i_i_factor = 1000, m_i_factor = 1000
params_subset$j_phi_i_i_factor <- 1000
params_subset$m_i_factor <- 1000

p_spatial_1000_1000 <- spatial_plot_subset(
  plot_times,
  params_subset,
  colour_by = pe,
  linetype_by = gamma
)

ggsave_with_defaults(
  plot = p_spatial_1000_1000,
  paste(plot_dir, "spatial_1000_1000.pdf", sep = "/"),
  width = 8.5,
  height = 14
)

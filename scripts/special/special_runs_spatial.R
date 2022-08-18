library(ggplot2)
library(cowplot)
library(patchwork)

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

spatial_plot_subset <- function(time,
                                params_subset,
                                colour_by,
                                linetype_by,
                                alpha_by) {
  data <- read_spatial_data(params, time, res_dir_base)
  data_long <- melt(
    data,
    measure.vars = c("C_b", "phi_C_b")
  )
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

  p <- p +
    geom_line(size = 1) +
    facet_wrap(
      vars(variable),
      scales = "free",
      labeller = as_labeller(function(v) all_labels[v], label_parsed),
      strip.position = "left"
    ) +
    theme(strip.background = element_blank(), strip.placement = "outside") +
    labs(
      x = expression(x),
      y = NULL,
      title = paste("Time since inflammation: ", data[, time_inf] %>% unique),
      colour = param_labels[deparse(substitute(colour_by))]
    )
  p
}

spatial_plot_subset_2 <- function(time,
                                  params_subset,
                                  colour_by,
                                  linetype_by,
                                  alpha_by) {
  data <- read_spatial_data(params, time, res_dir_base)

  spatial_plot_panel <- function(variable) {
    p <- ggplot(
      data[
        j_phi_i_i_factor %in% as.numeric(params_subset$j_phi_i_i_factor) &
          m_i_factor %in% as.numeric(params_subset$m_i_factor) &
          t_j_phi_i_lag %in% as.numeric(params_subset$t_j_phi_i_lag) &
          gamma %in% as.numeric(params_subset$gamma) &
          pe %in% as.numeric(params_subset$pe)
      ],
      aes(
        x = x,
        y = {{ variable }},
        group = rep,
        colour = factor({{ colour_by }})
      )
    )

    p <- p +
      geom_line(size = 1) +
      theme(strip.background = element_blank(), strip.placement = "outside") +
      labs(
        x = expression(x),
        y = variable_labels[deparse(substitute(variable))]
      )
    p
  }

  plots <- list()

  plots[[1]] <- spatial_plot_panel(C_b) +
    labs(
      title = paste(
        "Time since inflammation: ",
        round(data[, time_inf] %>% unique, 3)
      ),
      colour = param_labels[deparse(substitute(colour_by))]
    )
  plots[[2]] <- spatial_plot_panel(phi_C_b) +
    labs(colour = param_labels[deparse(substitute(colour_by))])

  if (!missing(linetype_by)) {
    for (i in 1:length(plots)) {
      plots[[i]] <- plots[[i]] + aes(linetype = factor({{ linetype_by }})) +
        labs(linetype = param_labels[deparse(substitute(linetype_by))])
    }
  }

  if (!missing(alpha_by)) {
    # scale_alpha_manual only works for a parameter with at most the number of
    # values specified here
    for (i in 1:length(plots)) {
      plots[[i]] <- plots[[i]] + aes(alpha = factor({{ alpha_by }})) +
        scale_alpha_manual(values = c(1.0, 0.4)) +
        labs(alpha = param_labels[deparse(substitute(alpha_by))])
    }
  }

  plots
}

# general

spatial_plot_general <- function(params_subset) {
  p <- lapply(c(0, 50, 100, 150, 250, 500), function(t_inf) {
    spatial_plot_subset_2(
      t_inf,
      params_subset,
      colour_by = pe,
      linetype_by = gamma
    )
  })

  p <- wrap_plots(
    rlang::flatten(p),
    ncol = 2,
    guides = "collect"
    ) +
  plot_annotation(tag_levels = "A")

  p
}

spatial_plot_general(
  list(
    j_phi_i_i_factor = 1000,
    m_i_factor = 2,
    t_j_phi_i_lag = 25,
    gamma = c(0, 1),
    pe = c(-5, -3, -1, 1, 3, 5)
  )
)

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

p_spatial_2_2 <- spatial_plot_general(params_subset)

ggsave_with_defaults(
  plot = p_spatial_2_2,
  paste(plot_dir, "spatial_2_2.pdf", sep = "/"),
  width = 8.5,
  height = 14
)

# j_phi_i_i_factor = 1000, m_i_factor = 2
params_subset$j_phi_i_i_factor <- 1000
params_subset$m_i_factor <- 2

p_spatial_1000_2 <- spatial_plot_general(params_subset)

ggsave_with_defaults(
  plot = p_spatial_1000_2,
  paste(plot_dir, "spatial_1000_2.pdf", sep = "/"),
  width = 8.5,
  height = 14
)

# j_phi_i_i_factor = 2, m_i_factor = 1000
params_subset$j_phi_i_i_factor <- 2
params_subset$m_i_factor <- 1000

p_spatial_2_1000 <- spatial_plot_general(params_subset)

ggsave_with_defaults(
  plot = p_spatial_2_1000,
  paste(plot_dir, "spatial_2_1000.pdf", sep = "/"),
  width = 8.5,
  height = 14
)

# j_phi_i_i_factor = 1000, m_i_factor = 1000
params_subset$j_phi_i_i_factor <- 1000
params_subset$m_i_factor <- 1000

p_spatial_1000_1000 <- spatial_plot_general(params_subset)

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

p_spatial_2_2 <- spatial_plot_general(params_subset)

ggsave_with_defaults(
  plot = p_spatial_2_2,
  paste(plot_dir, "spatial_2_2.pdf", sep = "/"),
  width = 8.5,
  height = 14
)

# j_phi_i_i_factor = 1000, m_i_factor = 2
params_subset$j_phi_i_i_factor <- 1000
params_subset$m_i_factor <- 2

p_spatial_1000_2 <- spatial_plot_general(params_subset)

ggsave_with_defaults(
  plot = p_spatial_1000_2,
  paste(plot_dir, "spatial_1000_2.pdf", sep = "/"),
  width = 8.5,
  height = 14
)

# j_phi_i_i_factor = 2, m_i_factor = 1000
params_subset$j_phi_i_i_factor <- 2
params_subset$m_i_factor <- 1000

p_spatial_2_1000 <- spatial_plot_general(params_subset)

ggsave_with_defaults(
  plot = p_spatial_2_1000,
  paste(plot_dir, "spatial_2_1000.pdf", sep = "/"),
  width = 8.5,
  height = 14
)

# j_phi_i_i_factor = 1000, m_i_factor = 1000
params_subset$j_phi_i_i_factor <- 1000
params_subset$m_i_factor <- 1000

p_spatial_1000_1000 <- spatial_plot_general(params_subset)

ggsave_with_defaults(
  plot = p_spatial_1000_1000,
  paste(plot_dir, "spatial_1000_1000.pdf", sep = "/"),
  width = 8.5,
  height = 14
)

library(sensitivity)
library(ggplot2)
library(cowplot)
library(patchwork)
library(GGally)
library(ggnewscale)

source("scripts/sensitivity/functions.R")

# Set the random seed to make the parameter sample reproducible. This will be
# useful for setting off some runs before we've fully decided what quantities
# to extract from the simulations. Once we've decided, we can then re-run this
# script, skipping running the simulations and just calculate the Sobol indices
# etc.
set.seed(12345L)

res_dir_base <- "res_sensitivity_3_neg_pe"

# Set and create the plot directory if it doesn't exist
plot_dir <- paste(res_dir_base, "plots", sep = "/")

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# Parameter and variable labels as `expression`s for nicer formatting in plots
param_labels <- c(
  "j_phi_i_i_factor" = expression(paste(J[phi[i]]^I, " factor")),
  "m_i_factor" = expression(paste(M^I, " factor")),
  "t_j_phi_i_lag" = expression(paste(J[phi[i]]^I, " delay")),
  "gamma" = expression(gamma)
)

# Add any other variables as needed
variable_labels <- c(
  "C_b^{tot}" = expression(paste("Total ", C[b])),
  "phi_{C_b}^{tot}" = expression(paste("Total ", phi[C[b]])),
  "-F_{phi_{C_b}}(x=0)" = expression(paste("Cell flux at l.v."))
)

all_labels <- c(param_labels, variable_labels)

# Min/max for the parameters we're varying
j_phi_i_i_factor_min <- 1.0
j_phi_i_i_factor_max <- 1000.0
m_i_factor_min <- 1.0
m_i_factor_max <- 1000.0
t_j_phi_i_lag_min <- 0.0
t_j_phi_i_lag_max <- 25.0
gamma_min <- 0.0
gamma_max <- 10.0

names <- c("j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag", "gamma")
mins <- c(j_phi_i_i_factor_min, m_i_factor_min, t_j_phi_i_lag_min, gamma_min)
maxs <- c(j_phi_i_i_factor_max, m_i_factor_max, t_j_phi_i_lag_max, gamma_max)

x_1 <- gen_param_sample(100, names, mins, maxs)
x_2 <- gen_param_sample(100, names, mins, maxs)

#x <- soboljansen(model = NULL, X1 = x_1, X2 = x_2, nboot = 100)
x <- sobolmartinez(model = NULL, X1 = x_1, X2 = x_2, nboot = 100)

add_constants_and_gammas_to_param_table(x$X)

#write_json_params_files(x$X, res_dir_base)
#fwrite(x$X, paste(res_dir_base, "d_m.csv", sep = "/"), sep = " ")

#simulate_all(x$X, res_dir_base)

# TODO: handle the fact that we will sometimes only want to run the simulations
# and not do the sensitivity analysis in the same run, and vice-versa. In the
# latter we'll need to skip writing the parameters to file, skip running the
# sims and just read the combined outputs from the file. Currently
# `run_all.bash` is responsible for running sims and combining outputs, so
# maybe we should split this into two scripts to make it easier to call them as
# needed from R.

trace_data_full <- read_trace_data(x$X, res_dir_base)
trace_data <- trace_data_full[`t_{inf}` >= 0]

trace_data_long <- melt(
  trace_data,
  measure.vars = c(
    "C_u^{tot}",
    "C_b^{tot}",
    "C_s^{tot}",
    "phi_i^{tot}",
    "phi_m^{tot}",
    "phi_{C_u}^{tot}",
    "phi_{C_b}^{tot}",
    "-F_{phi_i}(x=1)",
    "-F_{phi_{C_b}}(x=0)"
  )
)

integrated_fluxes <- calculate_integrated_fluxes(trace_data)
integrated_fluxes[, net_change := cells_in - cells_out]

cells_vs_params <- x$X[
  ,
  .(j_phi_i_i_factor,
    m_i_factor,
    t_j_phi_i_lag,
    gamma,
    cells_in = integrated_fluxes[, cells_in],
    cells_out = integrated_fluxes[, cells_out],
    net_change = integrated_fluxes[, net_change]
  )
]

cells_vs_params_long <- cells_vs_params %>%
  melt(
    .,
    measure.vars = c(
      "j_phi_i_i_factor",
      "m_i_factor",
      "t_j_phi_i_lag",
      "gamma"
    ),
    variable.name = "param",
    value.name = "param_value"
  ) %>%
  melt(
    .,
    measure.vars = c("cells_in", "cells_out", "net_change"),
    variable.name = "variable",
    value.name = "cells"
  )

sobol_indices_cells_plot <- function(x, title) {
  rn <- rownames(x$T)

  data <- rbind(
    cbind(
      parameter = factor(rn, levels = rn),
      effect = "main",
      x$S %>% as.data.table()
    ),
    cbind(
      parameter = factor(rn, levels = rn),
      effect = "total",
      x$T %>% as.data.table()
    )
  )

  ggplot(
    data,
    aes(
      parameter,
      y = original,
      ymin = `min. c.i.`,
      ymax = `max. c.i.`,
      group = effect,
      shape = effect,
      colour = effect)
    ) +
    geom_point(size = 3, position = position_dodge(width = 0.3)) +
    geom_errorbar(width = 0.2, position = position_dodge(width = 0.3)) +
    scale_x_discrete(labels = param_labels) +
    scale_shape_discrete(labels = c("main" = "Main", "total" = "Total")) +
    scale_colour_discrete(labels = c("main" = "Main", "total" = "Total")) +
    coord_cartesian(ylim = c(0.0, 1.0)) +
    labs(
      x = NULL,
      y = "Sobol index",
      shape = "Effect",
      colour = "Effect",
      title = title
    ) +
    theme_cowplot() +
    theme(plot.title = element_text(hjust = 0.5))
}

# Sobol indices for cells in
# --------------------------
y <- integrated_fluxes$cells_in

tell(x, y)
print(x)

p_sobol_indices_cells_in <- sobol_indices_cells_plot(x, "Cells in")

# Sobol indices for cells out
# ---------------------------
y <- integrated_fluxes$cells_out

tell(x, y)
print(x)

p_sobol_indices_cells_out <- sobol_indices_cells_plot(x, "Cells out")

# Combined plot for cells in and out
#-----------------------------------

p_sobol_indices_cells <- p_sobol_indices_cells_in +
  p_sobol_indices_cells_out +
  plot_layout(guides = "collect")

ggsave_with_defaults(
  plot = p_sobol_indices_cells,
  paste(plot_dir, "sobol_indices_cells.pdf", sep = "/")
)

# Sobol indices for net change in number of cells
# -----------------------------------------------
y <- integrated_fluxes$net_change

tell(x, y)
print(x)
sobol_indices_cells_plot(x, "Net change in cells")

# Other plots
# -----------

# Plot the integrated cell numbers against parameter values
p_cells_vs_params <- ggplot(cells_vs_params_long) +
  geom_point(aes(
    x = param_value,
    y = cells,
    colour = variable,
    group = variable,
  ), size = 1.5) +
  facet_wrap(
    vars(param),
    scales = "free",
    labeller = as_labeller(function(v) param_labels[v], label_parsed),
    strip.position = "bottom"
  ) +
  scale_color_discrete(
    labels = c(
      "cells_in" = "Cells in",
      "cells_out" = "Cells out",
      "net_change" = "Net change"
    )
  ) +
  theme_cowplot() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  labs(colour = NULL, x = NULL, y = "Cells (dimensionless)")

ggsave_with_defaults(
  plot = p_cells_vs_params,
  paste(plot_dir, "cells_vs_params.pdf", sep = "/")
)

# Individual cells in/out plots
#ggplot(cells_vs_params_long[variable == "cells_in"]) +
  #geom_point(
    #aes(x = param_value, y = cells, colour = param, group = variable)
  #) +
  #facet_wrap(vars(param), scales = "free") +
  #theme_cowplot() +
  #xlab("Parameter value") +
  #ylab("Cells (dimensionless)")

#ggplot(cells_vs_params_long[variable == "cells_out"]) +
  #geom_point(
    #aes(x = param_value, y = cells, colour = param, group = variable)
  #) +
  #facet_wrap(vars(param), scales = "free") +
  #theme_cowplot() +
  #xlab("Parameter value") +
  #ylab("Cells (dimensionless)")

colour_bar_gradient_ordered <- function(high, order) {
  scale_colour_gradient(
    low = "black",
    high = high,
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      order = order
    )
  )
}
blue <- colour_bar_gradient_ordered("blue", 1)
red <- colour_bar_gradient_ordered("red", 2)
green <- colour_bar_gradient_ordered("green", 3)
orange <- colour_bar_gradient_ordered("orange", 4)
colour_scales <- c(
  "j_phi_i_i_factor" = blue,
  "m_i_factor" = red,
  "t_j_phi_i_lag" = green,
  "gamma" = orange
)

#grid_panel <- function(y_var, colour_by, colour_style, axis_style = NULL) {
  #p <- ggplot(
      #trace_data[rep %in% 0:9],
      #aes(
        #x = `t_{inf}`,
        #y = {{ y_var }},
        #group = rep,
        #colour = {{ colour_by }}
      #)
    #) +
    #geom_line(alpha = 0.2) +
    #labs(x = "Time since inflammation") +
    #theme_cowplot() +
    #theme(legend.position = "top") +
    #colour_style
  #if (!is.null(axis_style)) {
    #p <- p + axis_style
  #}
  #p
#}

trace_data_longer <- trace_data[
  ,
  .(
    rep, `t_{inf}`, `C_b^{tot}`, `phi_{C_b}^{tot}`, `-F_{phi_{C_b}}(x=0)`,
    j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma
  )
  ] %>%
  melt(
    .,
    measure.vars = c(
      "C_b^{tot}", "phi_{C_b}^{tot}", "-F_{phi_{C_b}}(x=0)"
    )
  ) %>%
  melt(
    .,
    measure.vars = c(
      "j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag", "gamma"
    ),
    variable.name = "param",
    value.name = "param_value"
  )

p_trace_grid <- ggplot(
  trace_data_longer,
  aes(
    x = `t_{inf}`,
    y = value,
  )
)
for (param_name in trace_data_longer$param %>% levels) {
  p_trace_grid <- p_trace_grid +
    geom_line(
      aes(
        group = rep,
        colour = param_value
      ),
      alpha = 0.2,
      data = trace_data_longer[param == param_name]
    ) +
    colour_scales[param_name] +
    labs(colour = param_labels[param_name]) +
    new_scale_colour()
}
p_trace_grid <- p_trace_grid +
  labs(x = "Time since inflammation", y = NULL) +
  facet_grid(
    rows = vars(variable),
    cols = vars(param),
    scales = "free",
    switch = "y",
    labeller = as_labeller(function(v) all_labels[v], label_parsed)
  ) +
  theme_cowplot() +
  theme(
    plot.background = element_rect("white"),
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.box.just = "bottom",
    legend.position = "top",
    legend.justification = "centre",
    # no idea why 1/24 seems to be correct here given that "npc" units mean
    # that surely it should be 0.25, but who cares
    legend.key.width = unit(1.0 / 24.0, "npc"),
    legend.key.height = unit(0.3, "cm")
  )

ggsave_with_defaults(
  paste(plot_dir, "trace_grid.png", sep = "/"),
  plot = p_trace_grid
)

############################

# Plot some specific repeats from interesting parts of paramteter space

interesting_runs <- c(138, 347, 68, 419)

ggplot(trace_data_long[rep %in% interesting_runs]) +
  geom_line(aes(x = `t_{inf}`, y = value, colour = factor(rep)), size = 2) +
  facet_wrap(vars(variable), scales = "free") +
  theme_cowplot()

grid_panel(`phi_i^{tot}` + `phi_m^{tot}` + `phi_{C_u}^{tot}` + `phi_{C_b}^{tot}`, j_phi_i_i_factor, blue)

############################

# Sobol indices as a function of time - general methods

t_inf_max <- trace_data %>%
    .[, .(t_stop = max(`t_{inf}`)), by = rep] %>%
    .[, min(t_stop)]

output_inf_max <- trace_data %>%
  .[, .(output_stop = max(output_inf)), by = rep] %>%
  .[, min(output_stop)]

get_t_inf <- function(i) {
  trace_data[rep == i & `t_{inf}` <= t_inf_max, .(`t_{inf}`)]
}



sobol_at_time <- function(data, output, variable) {
  y <- data[output_inf == output, ..variable] %>% unlist
  #y <- y - mean(y)
  tell(x, y)
  list(S = x$S, T = x$T)
}

# Sobol indices of cell outflux as a function of time

flux_first_order_sobol_indices <- list()
flux_total_order_sobol_indices <- list()

for (i in 1:output_inf_max) {
  print(i)
  st <- sobol_at_time(trace_data, i, "-F_{phi_{C_b}}(x=0)")
  s <- st$S
  t <- st$T

  s <- cbind(rownames(s), s %>% as.data.table)
  setnames(s, "V1", "variable")
  t <- cbind(rownames(t), t %>% as.data.table)
  setnames(t, "V1", "variable")
  t_inf <- trace_data[output_inf == i & rep == 1, `t_{inf}`]
  s[, `t_{inf}` := t_inf]
  t[, `t_{inf}` := t_inf]

  flux_first_order_sobol_indices[[i]] <- s
  flux_total_order_sobol_indices[[i]] <- t
}

flux_first_order_sobol_indices <- rbindlist(flux_first_order_sobol_indices)
flux_total_order_sobol_indices <- rbindlist(flux_total_order_sobol_indices)

p_flux_first_order <- ggplot(
  flux_first_order_sobol_indices,
  aes(x = `t_{inf}`, y = original, group = variable, colour = variable)
) +
  geom_ribbon(
    aes(ymin = `min. c.i.`, ymax = `max. c.i.`, fill = variable, colour = NULL),
    alpha = 0.2
  ) +
  geom_line() +
  labs(
    x = "Time since inflammation",
    y = "First-order Sobol index",
    colour = "Parameter",
    fill = "Parameter"
  ) +
  scale_colour_discrete(labels = param_labels) +
  scale_fill_discrete(labels = param_labels) +
  coord_cartesian(ylim = c(-0.1, 1.0)) +
  theme_cowplot()

p_flux_first_order

ggsave_with_defaults(
  plot = p_flux_first_order,
  paste(plot_dir, "flux_first_order_sobol_indices.pdf", sep = "/"),
)

p_flux_total_order <- ggplot(
  flux_total_order_sobol_indices,
  aes(x = `t_{inf}`, y = original, group = variable, colour = variable)
) +
  geom_ribbon(
    aes(ymin = `min. c.i.`, ymax = `max. c.i.`, fill = variable, colour = NULL),
    alpha = 0.2
  ) +
  geom_line() +
  labs(x = "Time since inflammation", y = "Total-order Sobol index") +
  coord_cartesian(ylim = c(0, 1.2)) +
  theme_cowplot()

p_flux_total_order

ggsave_with_defaults(
  plot = p_flux_total_order,
  paste(plot_dir, "flux_total_order_sobol_indices.pdf", sep = "/")
)

################

# Sobol indices of maximum gradient of C_b as a function of time

# This section currently requires data from
# `scripts/extract_max_cell_density.R` to be loaded into the session
# TODO: refactor

gradient_first_order_sobol_indices <- list()
gradient_total_order_sobol_indices <- list()

for (i in 1:output_inf_max) {
  print(i)
  st <- sobol_at_time(max_dc_b_dx_all_and_params, i, "dC_b_dx")
  s <- st$S
  t <- st$T

  s <- cbind(rownames(s), s %>% as.data.table)
  setnames(s, "V1", "variable")
  t <- cbind(rownames(t), t %>% as.data.table)
  setnames(t, "V1", "variable")
  t_inf <- trace_data[output_inf == i & rep == 1, `t_{inf}`]
  s[, `t_{inf}` := t_inf]
  t[, `t_{inf}` := t_inf]

  gradient_first_order_sobol_indices[[i]] <- s
  gradient_total_order_sobol_indices[[i]] <- t
}

gradient_first_order_sobol_indices <- rbindlist(gradient_first_order_sobol_indices)
gradient_total_order_sobol_indices <- rbindlist(gradient_total_order_sobol_indices)

p_gradient_first_order <- ggplot(
  gradient_first_order_sobol_indices,
  aes(x = `t_{inf}`, y = original, group = variable, colour = variable)
) +
  geom_ribbon(
    aes(ymin = `min. c.i.`, ymax = `max. c.i.`, fill = variable, colour = NULL),
    alpha = 0.2
  ) +
  geom_line() +
  labs(x = "Time since inflammation", y = "First-order Sobol index") +
  coord_cartesian(ylim = c(-0.1, 1.0)) +
  theme_cowplot()

p_gradient_first_order

ggsave_with_defaults(
  plot = p_gradient_first_order,
  paste(plot_dir, "gradient_first_order_sobol_indices.pdf", sep = "/")
)

p_gradient_total_order <- ggplot(
  gradient_total_order_sobol_indices,
  aes(x = `t_{inf}`, y = original, group = variable, colour = variable)
) +
  geom_ribbon(
    aes(ymin = `min. c.i.`, ymax = `max. c.i.`, fill = variable, colour = NULL),
    alpha = 0.2
  ) +
  geom_line() +
  labs(x = "Time since inflammation", y = "Total-order Sobol index") +
  coord_cartesian(ylim = c(0, 1.2)) +
  theme_cowplot()

p_gradient_total_order

ggsave_with_defaults(
  plot = p_gradient_total_order,
  paste(plot_dir, "gradient_total_order_sobol_indices.pdf", sep = "/")
)

library(sensitivity)
library(ggplot2)
library(cowplot)
library(patchwork)
library(GGally)
library(ggnewscale)
library(ggh4x)
library(ggbreak)

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

# add dimensional time so we can plot against it
trace_data[, t_inf_h := 0.6944444444 * `t_{inf}`]

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
    scale_shape_discrete(labels = c("main" = "First", "total" = "Total")) +
    scale_colour_discrete(labels = c("main" = "First", "total" = "Total")) +
    coord_cartesian(ylim = c(0.0, 1.0)) +
    labs(
      x = NULL,
      y = "Sobol index",
      shape = "Order",
      colour = "Order",
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
p_sobol_indices_cells_net_change <- sobol_indices_cells_plot(x, "Net change in cells")

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
  geom_hline(yintercept = 0, linetype = "dashed") +
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

## Plot some specific repeats from interesting parts of paramteter space

#interesting_runs <- c(138, 347, 68, 419)

#ggplot(trace_data_long[rep %in% interesting_runs]) +
  #geom_line(aes(x = `t_{inf}`, y = value, colour = factor(rep)), size = 2) +
  #facet_wrap(vars(variable), scales = "free") +
  #theme_cowplot()

#grid_panel(`phi_i^{tot}` + `phi_m^{tot}` + `phi_{C_u}^{tot}` + `phi_{C_b}^{tot}`, j_phi_i_i_factor, blue)

############################

# Sobol indices as a function of time - general methods

t_inf_max <- trace_data %>%
    .[, .(t_stop = max(`t_{inf}`)), by = rep] %>%
    .[, min(t_stop)]

output_inf_max <- trace_data %>%
  .[, .(output_stop = max(output_inf)), by = rep] %>%
  .[, min(output_stop)]

# checks that the t_inf values for different reps are essentially equal
check_t_infs <- function() {
  get_t_inf <- function(i) {
    trace_data[rep == i & `t_{inf}` <= t_inf_max, .(`t_{inf}`)]
  }

  compare_t_infs <- function(i, j) {
    (abs(get_t_inf(i) - get_t_inf(j)) <= 1e-4) %>% all
  }

  for (i in 2:599) {
    if (!compare_t_infs(1, i)) {
      print(paste0("1 and ", i, " t_infs differ more than 1e-4"))
    }
  }
}

# check_t_infs()

sobol_at_time <- function(data, output, variable) {
  y <- data[output_inf == output, ..variable] %>% unlist
  #y <- y - mean(y)
  tell(x, y)
  list(S = x$S, T = x$T)
}

p_sobol_vs_time <- function(sobol_data, title, ylim) {
  ggplot(
    sobol_data,
    aes(x = `t_{inf}`, y = original, group = variable, colour = variable)
  ) +
    geom_ribbon(
      aes(
        ymin = `min. c.i.`,
        ymax = `max. c.i.`,
        fill = variable,
        colour = NULL
      ),
      alpha = 0.2
    ) +
    geom_line() +
    labs(
      x = "Time since inflammation",
      y = title,
      colour = "Parameter",
      fill = "Parameter"
    ) +
    scale_colour_discrete(labels = param_labels) +
    scale_fill_discrete(labels = param_labels) +
    coord_cartesian(ylim = ylim) +
    theme_cowplot()
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

# First order

p_flux_first_order <- p_sobol_vs_time(
  flux_first_order_sobol_indices,
  "First-order Sobol index",
  c(-0.1, 1.0)
)

ggsave_with_defaults(
  plot = p_flux_first_order,
  paste(plot_dir, "flux_first_order_sobol_indices.pdf", sep = "/"),
)

# Total order

p_flux_total_order <- p_sobol_vs_time(
  flux_total_order_sobol_indices,
  "Total-order Sobol index",
  c(0, 1.2)
)

ggsave_with_defaults(
  plot = p_flux_total_order,
  paste(plot_dir, "flux_total_order_sobol_indices.pdf", sep = "/")
)

################

# Sobol indices of maximum gradient of C_b as a function of time

source("scripts/extract_max_cell_density.R")

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

# First order

p_gradient_first_order <- p_sobol_vs_time(
  gradient_first_order_sobol_indices,
  "First-order Sobol index",
  c(-0.1, 1.0)
)

ggsave_with_defaults(
  plot = p_gradient_first_order,
  paste(plot_dir, "gradient_first_order_sobol_indices.pdf", sep = "/")
)

# Total order

p_gradient_total_order <- p_sobol_vs_time(
  gradient_total_order_sobol_indices,
  "Total-order Sobol index",
  c(0.0, 1.2)
)

ggsave_with_defaults(
  plot = p_gradient_total_order,
  paste(plot_dir, "gradient_total_order_sobol_indices.pdf", sep = "/")
)

# Max value and location of max value plots for phi_c_b and dc_b_dx
# -----------------------------------------------------------------

# Some plotting styles and helper functions

p_location_panel <- function(data, colour_by, colour_style) {
  print(deparse(substitute(colour_by)))
  print(param_labels[deparse(substitute(colour_by))])

  ggplot(
    data,
    aes(time_inf, x, group = rep, colour = {{ colour_by }})
    ) +
  geom_path(alpha = 0.1) +
  coord_cartesian(ylim = c(0, 1)) +
  colour_style +
  labs(
    x = "Time since inflammation",
    y = expression(x),
    colour = param_labels[deparse(substitute(colour_by))]
  )
}

p_value_panel <- function(data, colour_by, colour_style, ylabel) {
  print(deparse(substitute(colour_by)))
  print(param_labels[deparse(substitute(colour_by))])

  ggplot(
    data,
    aes(time_inf, phi_C_b, group = rep, colour = {{ colour_by }})
    ) +
  geom_path(alpha = 0.1) +
  colour_style +
  labs(
    x = "Time since inflammation",
    y = ylabel,
    colour = param_labels[deparse(substitute(colour_by))]
  )
}

# Location of max phi_c_b

p_j_phi_i_i_factor <- p_location_panel(max_phi_c_b_all_and_params, j_phi_i_i_factor, blue)
p_m_i_factor <- p_location_panel(max_phi_c_b_all_and_params, m_i_factor, red)
p_t_j_phi_i_lag <- p_location_panel(max_phi_c_b_all_and_params, t_j_phi_i_lag, green)
p_gamma <- p_location_panel(max_phi_c_b_all_and_params, gamma, orange)

p_phi_c_b_location <- p_j_phi_i_i_factor + p_m_i_factor + p_t_j_phi_i_lag + p_gamma +
  plot_annotation(
    title = expression(
      paste("Location of maxmimum ", phi[C[b]], " concentration")
    )
  ) &
  theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect("white")
  )

ggsave_with_defaults(
  plot = p_phi_c_b_location,
  paste(plot_dir, "max_phi_c_b_location.png", sep = "/")
)

# Value of max phi_c_b

p_phi_c_b_value_panel <- function(...) {
  p_value_panel(..., expression(phi[C[b]]))
}

p_j_phi_i_i_factor <- p_phi_c_b_value_panel(max_phi_c_b_all_and_params, j_phi_i_i_factor, blue)
p_m_i_factor <- p_phi_c_b_value_panel(max_phi_c_b_all_and_params, m_i_factor, red)
p_t_j_phi_i_lag <- p_phi_c_b_value_panel(max_phi_c_b_all_and_params, t_j_phi_i_lag, green)
p_gamma <- p_phi_c_b_value_panel(max_phi_c_b_all_and_params, gamma, orange)

p_phi_c_b_value <- p_j_phi_i_i_factor + p_m_i_factor + p_t_j_phi_i_lag + p_gamma +
  plot_annotation(
    title = expression(
      paste("Value of maxmimum ", phi[C[b]], " concentration")
    )
  ) &
  theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect("white")
  )

ggsave_with_defaults(
  plot = p_phi_c_b_value,
  paste(plot_dir, "max_phi_c_b_value.png", sep = "/")
)

# Location of max gradient of c_b

p_j_phi_i_i_factor <- p_location_panel(max_dc_b_dx_all_and_params, j_phi_i_i_factor, blue)
p_m_i_factor <- p_location_panel(max_dc_b_dx_all_and_params, m_i_factor, red)
p_t_j_phi_i_lag <- p_location_panel(max_dc_b_dx_all_and_params, t_j_phi_i_lag, green)
p_gamma <- p_location_panel(max_dc_b_dx_all_and_params, gamma, orange)

p_dc_b_dx_location <- p_j_phi_i_i_factor + p_m_i_factor + p_t_j_phi_i_lag + p_gamma +
  plot_annotation(
    title = expression(
      paste("Location of maxmimum ", C[b], " gradient")
    )
  ) &
  theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect("white")
  )

ggsave_with_defaults(
  plot = p_dc_b_dx_location,
  paste(plot_dir, "max_dc_b_dx_location.png", sep = "/")
)

# Value of max gradient of c_b

p_dc_b_dx_value_panel <- function(...) {
  p_value_panel(..., expression(dC[b]/dx))
}

p_j_phi_i_i_factor <- p_dc_b_dx_value_panel(max_dc_b_dx_all_and_params, j_phi_i_i_factor, blue)
p_m_i_factor <- p_dc_b_dx_value_panel(max_dc_b_dx_all_and_params, m_i_factor, red)
p_t_j_phi_i_lag <- p_dc_b_dx_value_panel(max_dc_b_dx_all_and_params, t_j_phi_i_lag, green)
p_gamma <- p_dc_b_dx_value_panel(max_dc_b_dx_all_and_params, gamma, orange)

p_dc_b_dx_value <- p_j_phi_i_i_factor + p_m_i_factor + p_t_j_phi_i_lag + p_gamma +
  plot_annotation(
    title = expression(
      paste("Value of maxmimum ", C[b], " gradient")
    )
  ) &
  theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect("white")
  )

ggsave_with_defaults(
  plot = p_dc_b_dx_value,
  paste(plot_dir, "max_dc_b_dx_value.png", sep = "/")
)

################

# Get the rows corresponding to the (global) maximum cell flux across lymphatic
# vessel
trace_max_flux <- trace_data[, .SD[which.max(`-F_{phi_{C_b}}(x=0)`)], by = rep]

trace_max_flux_long <- melt(
  trace_max_flux,
  measure.vars = c(
    "j_phi_i_i_factor",
    "m_i_factor",
    "t_j_phi_i_lag",
    "gamma"),
  variable.name = "param",
  value.name = "param_value"
)

# Get the first local maximum assuming it's before t_inf = 50
trace_max_flux_first <-
  trace_data[`t_{inf}` < 49, .SD[which.max(`-F_{phi_{C_b}}(x=0)`)], by = rep]

trace_max_flux_first_long <- melt(
  trace_max_flux_first,
  measure.vars = c(
    "j_phi_i_i_factor",
    "m_i_factor",
    "t_j_phi_i_lag",
    "gamma"),
  variable.name = "param",
  value.name = "param_value"
)

# Get the second local maximum assuming it's after t_inf = 50
trace_max_flux_second <-
  trace_data[`t_{inf}` >= 49, .SD[which.max(`-F_{phi_{C_b}}(x=0)`)], by = rep]

trace_max_flux_second_long <- melt(
  trace_max_flux_second,
  measure.vars = c(
    "j_phi_i_i_factor",
    "m_i_factor",
    "t_j_phi_i_lag",
    "gamma"),
  variable.name = "param",
  value.name = "param_value"
)

plot_max_flux <- function(max_data) {
  p <- ggplot(
    max_data,
    aes(x = `t_{inf}`, y = `-F_{phi_{C_b}}(x=0)`)
  )

  panel <- function(param_str, new_scale_after = TRUE) {
    p <<- p +
      geom_point(
        aes(colour = param_value),
        data = max_data[param == param_str]
      ) +
      colour_scales[param_str] +
      labs(colour = param_labels[param_str])

    if (new_scale_after) {
      p <<- p + new_scale_colour()
    }
  }
  panel("j_phi_i_i_factor")
  panel("m_i_factor")
  panel("t_j_phi_i_lag")
  panel("gamma", new_scale_after = FALSE)

  p +
    #scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 5)) +
    facet_wrap(
      vars(param),
      labeller = as_labeller(function(v) all_labels[v], label_parsed)
    ) +
    #scale_x_break(c(10, 48)) +
    theme_cowplot() +
    background_grid() +
    #theme(
      #strip.background = element_blank(),
      #strip.text = element_blank()
    #) +
    labs(x = "Time since inflammation", y = "Max. cell flux at l.v.")
}

#plot_max_flux(trace_max_flux_second_long)
p_max_flux <- plot_max_flux(trace_max_flux_long)

ggsave_with_defaults(
  plot = p_max_flux,
  paste(plot_dir, "max_flux.pdf", sep = "/")
)

#################
# other way to do plots

max_flux_panel <- function(param) {
  ggplot(
    trace_max_flux,
    aes(x = `t_{inf}`, y = `-F_{phi_{C_b}}(x=0)`, colour = {{ param }})
  ) +
    geom_point(size = 2) +
    scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 5)) +
    scale_x_break(c(10, 48)) +
    theme_cowplot() +
    background_grid()
}

p_1 <- max_flux_panel(j_phi_i_i_factor) + blue
p_2 <- max_flux_panel(m_i_factor) + red
p_3 <- max_flux_panel(t_j_phi_i_lag) + green
p_4 <- max_flux_panel(gamma) + orange

#((p_1 | p_2) / (p_3 | p_4)) + plot_layout(guides = "collect")

aplot::plot_list(p_1, p_2, p_3, p_4, guides = "collect")

trace_max_flux_full_traj <-
  trace_data[rep %in% trace_max_flux[`t_{inf}` < 20, rep]]

ggplot(
  trace_max_flux_full_traj,
  aes(x = `t_{inf}`, y = `-F_{phi_{C_b}}(x=0)`, group = rep)
) +
  geom_line()

#################################

# Find all local maxima and minima of fluxes

# Helper functions
find_local_maxima_indices <- function(values) {
  max_ind <- NULL

  for (i in 2:(length(values) - 1)) {
    if (values[i] > values[i - 1] && values[i] > values[i + 1]) {
      max_ind <- c(max_ind, i)
    }
  }
  max_ind
}

find_local_minima_indices <- function(values) {
  min_ind <- NULL

  for (i in 2:(length(values) - 1)) {
    if (values[i] < values[i - 1] && values[i] < values[i + 1]) {
      min_ind <- c(min_ind, i)
    }
  }
  min_ind
}

flux_trace_plot_subset <- function(subset_rep_ids) {
  ggplot(
    trace_data[rep %in% subset_rep_ids],
    aes(
      x = t_inf_h,
      y = `-F_{phi_{C_b}}(x=0)`,
      group = rep
    )
  ) +
  geom_line() +
  theme_cowplot() +
  background_grid() +
  labs(
    x = "Time since inflammation (h)",
    y = "Cell flux at l.v."
  )
}

parameter_pairs_plot <- function(data) {
  ggpairs(
    data,
    columns = c("j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag", "gamma"),
    upper = "blank",
    diag = list(continuous = "barDiag"),
    progress = FALSE,
    labeller = as_labeller(function(v) param_labels[v], label_parsed)
  ) +
  theme_cowplot() +
  background_grid()
}

# Find the local maxima of the flux for each rep
flux_local_maxima <- trace_data[,
  .SD[find_local_maxima_indices(`-F_{phi_{C_b}}(x=0)`)],
  by = rep
] %>%
  .[, extrema_id := 1:.N, by = rep] %>%
  .[, extrema_type := factor("maximum")]

# Find the local minima of the flux for each rep
flux_local_minima <- trace_data[,
  .SD[find_local_minima_indices(`-F_{phi_{C_b}}(x=0)`)],
  by = rep
] %>%
  .[, extrema_id := 1:.N, by = rep] %>%
  .[, extrema_type := factor("minimum")]

flux_local_extrema <- rbind(flux_local_maxima, flux_local_minima)

# Plot all extrema at once
ggplot(
  flux_local_extrema[extrema_id == 1],
  aes(
    x = t_inf_h,
    y = `-F_{phi_{C_b}}(x=0)`,
    colour = extrema_type,
    #shape = extrema_type,
    group = rep
  )
) +
  geom_point(size = 1) +
  geom_line(alpha = 0.2) +
  scale_x_continuous(limits = c(0, 39)) +
  scale_shape_manual(values = c(19, 6)) +
  theme_cowplot() +
  background_grid() +
  labs(
    x = "Time since inflammation (h)",
    y = "Cell flux at l.v.",
    colour = "Count",
    shape = "Type"
  ) #+
  #geom_point(data = flux_local_extrema[extrema_type == "minimum"])

# Plot all minima
ggplot(
  flux_local_minima,
  aes(x = `t_{inf}`, y = `-F_{phi_{C_b}}(x=0)`, colour = factor(extrema_id))
) +
  geom_point() +
  scale_x_continuous(limits = c(0, 55)) +
  theme_cowplot() +
  background_grid()

# Plot all maxima
ggplot(
  flux_local_maxima,
  aes(x = `t_{inf}`, y = `-F_{phi_{C_b}}(x=0)`, colour = factor(extrema_id))
) +
  geom_point() +
  geom_vline(xintercept = 8 - 2) +
  geom_vline(xintercept = 8 + 2) +
  scale_x_continuous(limits = c(0, 55)) +
  theme_cowplot() +
  background_grid()

# Find the reps that have only one local maximum
flux_single_max <- flux_local_maxima %>%
  .[, .SD[which.max(extrema_id)], by = rep] %>%
  .[extrema_id == 1]

# Find the reps that have only one local minimum
flux_single_min <- flux_local_minima %>%
  .[, .SD[which.max(extrema_id)], by = rep] %>%
  .[extrema_id == 1]

# Bindi's suggested max/min thresholds

# (i) Solutions that have only 1 max at early time, t1 (t1 = 20?)
flux_early_max_ids <- flux_single_max[`t_{inf}` < 20, rep]

# (ii) Solutions that have only 1 min at t2, where t1 < t2 < t3
# TODO what are sensible t1 and t3 here? in the previous case, the choice of t1
# makes no difference, as long as it's > 5...
flux_middle_min_ids <- flux_single_min[`t_{inf}` > 0, rep]

# (iii) Solutions that have only 1 max at the end of inflammation, t3 (t3 ~ 50)
flux_end_inf_max_ids <- flux_single_max[abs(`t_{inf}` - 50) < 1e-4, rep]

# (iv) Solutions that have only 1 max after the end of inflammation, "t4" (t_inf > 50 (+ a
# bit due to output times not strictly being evenly spaced)
flux_after_inf_max_ids <- flux_single_max[`t_{inf}` > 50.00002, rep]

flux_trace_plot_subset(flux_early_max_ids) +
  geom_point(
    data = flux_single_max[`t_{inf}` < 20],
    aes(x = `t_{inf}`, y = `-F_{phi_{C_b}}(x=0)`),
    colour = "red",
    shape = "cross"
  )

flux_trace_plot_subset(flux_middle_min_ids) +
  geom_point(
    data = flux_single_min,
    aes(x = `t_{inf}`, y = `-F_{phi_{C_b}}(x=0)`),
    colour = "red",
    shape = "cross"
  )

flux_trace_plot_subset(flux_end_inf_max_ids) +
  geom_point(
    data = flux_single_max[abs(`t_{inf}` - 50) < 1e-4],
    aes(x = `t_{inf}`, y = `-F_{phi_{C_b}}(x=0)`),
    colour = "red",
    shape = "cross"
  )

flux_trace_plot_subset(flux_after_inf_max_ids) +
  geom_point(
    data = flux_single_max[`t_{inf}` > 50.00002],
    aes(x = `t_{inf}`, y = `-F_{phi_{C_b}}(x=0)`),
    colour = "red",
    shape = "cross"
  )

parameter_pairs_plot(flux_single_min)
parameter_pairs_plot(trace_data[output_inf == 1])

#####

# Get the ids for the reps that have at least 2 maxima
reps_w_al_2_maxima <- flux_local_maxima[extrema_id == 2, rep]

# Get the actual rows from the flux table up to the 2nd maximum
flux_local_al_2_maxima <- flux_local_maxima[
  rep %in% reps_w_al_2_maxima & extrema_id <= 2
]

# Only keep the columns we need here. This makes transforming into wide form
# below easier because we don't need to account for columns that vary between
# rows for one rep...
flux_local_al_2_maxima <- flux_local_al_2_maxima[
  ,
  .(
    rep,
    `t_{inf}`,
    t_inf_h,
    j_phi_i_i_factor,
    m_i_factor,
    t_j_phi_i_lag,
    gamma,
    extrema_id,
    `-F_{phi_{C_b}}(x=0)`
  )
]

# dcast the two maxima into their own columns
flux_local_al_2_maxima_wide <- dcast(
  flux_local_al_2_maxima,
  ... ~ extrema_id,
  value.var = c("t_{inf}", "t_inf_h", "-F_{phi_{C_b}}(x=0)")
)

# calculate the ratio max_2/max_1
flux_local_al_2_maxima_wide[
  ,
  ratio_21 := `-F_{phi_{C_b}}(x=0)_2` / `-F_{phi_{C_b}}(x=0)_1`,
  by = rep
]

# "join" the ratio_21 column with the trace_data
trace_data_with_max_21_ratio <-
  trace_data[rep %in% flux_local_al_2_maxima_wide[, rep]][
  flux_local_al_2_maxima_wide[, .(rep, ratio_21)],
  on = "rep"
]

## convert to long form
#trace_data_with_ratios_long <- melt(
  #trace_data_with_max_21_ratio,
  #measure.vars = c(
    #"C_u^{tot}",
    #"C_b^{tot}",
    #"C_s^{tot}",
    #"phi_i^{tot}",
    #"phi_m^{tot}",
    #"phi_{C_u}^{tot}",
    #"phi_{C_b}^{tot}",
    #"-F_{phi_i}(x=1)",
    #"-F_{phi_{C_b}}(x=0)"
  #)
#)

ggplot(
  flux_local_al_2_maxima_wide,
  aes(x = ratio_21, colour = j_phi_i_i_factor)
) +
  geom_histogram(bins = 50)
  #geom_density()
  #geom_point(position = position_jitter(seed = 1))

ggplot(
  trace_data_with_max_21_ratio,
  aes(x = t_inf_h, y = `-F_{phi_{C_b}}(x=0)`, colour = ratio_21, group = rep)
) +
  geom_line(alpha = 0.5, size = 1) +
  scale_colour_distiller(palette = "Spectral") +
  labs(
    x = "Time since inflammation (h)",
    y = "Cell flux at l.v.",
    colour = "Max 2 / Max 1"
  ) +
  theme_cowplot() +
  background_grid()

###########

# Find the ids of all reps that have no local minima

# First get the reps that do have minima
reps_w_some_minima <- flux_local_minima[, .(rep)] %>% unique

# Then take the complement. There's definitely a better way to do this but this
# works...
reps_w_no_minima <-
  data.table(rep = 0:trace_data[, max(rep)])[!reps_w_some_minima, rep, on = "rep"]

flux_trace_plot_subset(reps_w_no_minima)

# We want the subset of the above reps that have a maximum at (or maybe soon
# after?) t_inf = 50, regardless of the total number of maxima

reps_w_no_min_but_max_50 <-
  flux_local_maxima[rep %in% reps_w_no_minima & `t_{inf}` >= 50, rep]

flux_trace_plot_subset(reps_w_no_min_but_max_50)

parameter_pairs_plot(trace_data[rep %in% reps_w_no_min_but_max_50])

###########

# Work out, for reps that have a minimum before t_inf == 50 (i.e. somewhere
# during inflammation), the ratio of the minimum to the first maximum. The
# first condition finds all reps that have any kind of extrema before t = 50,
# including some reps that have say a maximum but no minimum. We filter these
# out using na.omit to leave just the reps we actually want.
flux_w_max_min_ratio <- flux_local_extrema[
  `t_{inf}` < 50 & extrema_id == 1,
  .(rep, `-F_{phi_{C_b}}(x=0)`, extrema_type)
] %>%
  dcast(., ... ~ extrema_type, value.var = c("-F_{phi_{C_b}}(x=0)")) %>%
  .[, ratio_max_min := maximum / minimum] %>%
  na.omit

ggplot(flux_w_max_min_ratio, aes(x = ratio_max_min)) +
  #geom_density() +
  geom_histogram(bins = 50) +
  coord_cartesian(xlim = c(0, NA))

reps_w_min_and_max_50 <- flux_w_max_min_ratio[, rep]

# "join" the ratio_max_min column with the trace_data
trace_data_with_max_min_ratio <-
  trace_data[rep %in% reps_w_min_and_max_50][
  flux_w_max_min_ratio[, .(rep, ratio_max_min)],
  on = "rep"
]

# Plot the flux vs time, coloured by the ratio of max to min
ggplot(
  trace_data_with_max_min_ratio,
  aes(
    x = `t_{inf}`,
    y = `-F_{phi_{C_b}}(x=0)`,
    group = rep,
    colour = ratio_max_min
  )
) +
  geom_line(alpha = 0.5, size = 1) +
  scale_colour_distiller(palette = "Spectral") + theme_cowplot() +
  background_grid()

# Plot the flux vs time for reps that have ratio below a threshold
ggplot(
  trace_data_with_max_min_ratio[ratio_max_min <= 1.1],
  aes(
    x = `t_{inf}`,
    y = `-F_{phi_{C_b}}(x=0)`,
    group = rep,
    colour = ratio_max_min
  )
) +
  geom_line(alpha = 0.5, size = 1) +
  scale_colour_distiller(palette = "Spectral") + theme_cowplot() +
  background_grid()

# UPDATE: these aren't actually the relevant reps at all (I misunderstood...)
# Combine the above subset of reps with those that have no minimum but still
# have a max at 50. These should be the relevant ones???
reps_relevant <- c(
  flux_w_max_min_ratio[ratio_max_min <= 1.1, rep],
  reps_w_no_min_but_max_50
)

ggplot(
  trace_data[rep %in% reps_relevant],
  aes(
    x = `t_{inf}`,
    y = `-F_{phi_{C_b}}(x=0)`,
    group = rep,
    #colour = ratio_max_min
  )
) +
  geom_line(alpha = 0.5, size = 1) +
  scale_colour_distiller(palette = "Spectral") + theme_cowplot() +
  background_grid()

parameter_pairs_plot(trace_data[rep %in% reps_relevant])

###########

# Find the actually relevant reps - those that have a first peak at around 6
# hours (need to decide size of range around 6h to include), then have a second
# peak with height not too much lower than the first.

first_max_range_mean <- 6

# The time either side of 6h we allow
first_max_range_width <- 2.78
ratio_21_threshold <- 0.75
relevant_reps <- flux_local_al_2_maxima_wide[
  first_max_range_mean - first_max_range_width <= `t_inf_h_1` &
    `t_inf_h_1` <= first_max_range_mean + first_max_range_width &
    ratio_21 >= ratio_21_threshold,
  rep
]
flux_trace_plot_subset(relevant_reps) +
  geom_rect(
    data = ~ head(.x, 1),
    xmin = first_max_range_mean - first_max_range_width,
    xmax = first_max_range_mean + first_max_range_width,
    ymin = -1,
    ymax = 1,
    colour = NA,
    fill = "black",
    alpha = 0.1
  ) +
  geom_point(
    data = flux_local_al_2_maxima_wide[rep %in% relevant_reps],
    aes(x = `t_inf_h_1`, y = `-F_{phi_{C_b}}(x=0)_1`),
    colour = "red"
  ) +
  geom_point(
    data = flux_local_al_2_maxima_wide[rep %in% relevant_reps],
    aes(x = `t_inf_h_2`, y = `-F_{phi_{C_b}}(x=0)_2`),
    colour = "blue"
  ) +
  geom_vline(xintercept = 50 * 0.6944444444, linetype = "dashed")

#ggsave(
  #plot = last_plot(),
  #"simulation_subset.png",
  #width = 9,
  #height = 5
#)

parameter_pairs_plot(
  flux_local_al_2_maxima_wide[rep %in% relevant_reps]
)

# TODO the next couple of bits were copied from earlier - refactor into
# functions if we keep it
trace_data_relevant_longer <- trace_data[
  rep %in% relevant_reps,
  .(
    rep, `t_{inf}`, t_inf_h, `C_b^{tot}`, `phi_{C_b}^{tot}`,
    `-F_{phi_{C_b}}(x=0)`, j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma
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

p_relevant_trace_grid <- ggplot(
  trace_data_relevant_longer,
  aes(
    x = t_inf_h,
    y = value,
  )
)
for (param_name in trace_data_relevant_longer$param %>% levels) {
  p_relevant_trace_grid <- p_relevant_trace_grid +
    geom_line(
      aes(
        group = rep,
        colour = param_value
      ),
      alpha = 0.5,
      data = trace_data_relevant_longer[param == param_name]
    ) +
    colour_scales[param_name] +
    labs(colour = param_labels[param_name]) +
    new_scale_colour()
}
p_relevant_trace_grid <- p_relevant_trace_grid +
  labs(x = "Time since inflammation (h)", y = NULL) +
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
p_relevant_trace_grid

# Comparing with all reps (with >= 2 max) with m_i_factor approx in the range
# indicated by the above subset

flux_trace_plot_subset(
  flux_local_al_2_maxima_wide[m_i_factor %between% c(50, 200), rep]
)

setdiff(
  relevant_reps,
  flux_local_al_2_maxima_wide[m_i_factor %between% c(50, 200), rep]
)

# Plot ratio_21 against parameters

flux_local_al_2_maxima_wide_long_params <- melt(
  flux_local_al_2_maxima_wide,
  measure.vars = c("j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag", "gamma"),
  variable.name = "param",
  value.name = "param_value"
)

quantity_vs_param_plot <- function(quantity, colour_by = NULL) {
  p <- ggplot(
    # this join copies the un-melted parameter columns into the melted table so
    # we can use their values for e.g. colouring
    flux_local_al_2_maxima_wide_long_params[
      flux_local_al_2_maxima_wide[,
        .(rep, j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma)
      ],
      on = "rep"
    #],
    ][`t_{inf}_2` > 45], # & m_i_factor <= 200],
    aes(x = param_value, y = {{ quantity }}, colour = {{ colour_by }})
  ) +
    geom_point(size = 4, alpha = 0.6) +
    #geom_hline(yintercept = 1) +
    facet_wrap(vars(param), scales = "free") +
    scale_colour_distiller(palette = "Spectral") +
    scale_x_continuous(trans = "log10",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_continuous(trans = "log10",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    theme_cowplot() +
    background_grid()

  q_str <- gsub("[{}]", "", deparse(substitute(quantity)))
  c_str <- gsub("[{}]", "", deparse(substitute(colour_by)))
  plot_name <- paste0(q_str, "_by_", c_str, ".pdf")

  ggsave_with_defaults(
    plot = p,
    paste(plot_dir, "quantity_vs_param", plot_name, sep = "/")
  )

  p
}

quantity_vs_param_plot(t_inf_h_1, j_phi_i_i_factor)
quantity_vs_param_plot(t_inf_h_1, m_i_factor)
quantity_vs_param_plot(t_inf_h_1, t_j_phi_i_lag)
quantity_vs_param_plot(t_inf_h_1, gamma)
# => Time of first peak almost entirely controlled by m_i_factor. gamma appears
# to have a much smaller influence, with some low gamma reps having slightly higher t_inf_h_1

quantity_vs_param_plot(t_inf_h_2, j_phi_i_i_factor)
quantity_vs_param_plot(t_inf_h_2, m_i_factor)
quantity_vs_param_plot(t_inf_h_2, t_j_phi_i_lag)
quantity_vs_param_plot(t_inf_h_2, gamma)

quantity_vs_param_plot(`-F_{phi_{C_b}}(x=0)_1`, j_phi_i_i_factor)
quantity_vs_param_plot(`-F_{phi_{C_b}}(x=0)_1`, m_i_factor)
quantity_vs_param_plot(`-F_{phi_{C_b}}(x=0)_1`, t_j_phi_i_lag)
quantity_vs_param_plot(`-F_{phi_{C_b}}(x=0)_1`, gamma)

quantity_vs_param_plot(`-F_{phi_{C_b}}(x=0)_2`, j_phi_i_i_factor)
quantity_vs_param_plot(`-F_{phi_{C_b}}(x=0)_2`, m_i_factor)
quantity_vs_param_plot(`-F_{phi_{C_b}}(x=0)_2`, t_j_phi_i_lag)
quantity_vs_param_plot(`-F_{phi_{C_b}}(x=0)_2`, gamma)

interp_irregular_data <- function(x, y, z) {
  grid <- interp::interp(x, y, z, duplicate = "drop")
  subset(
    data.table::data.table(
      x = rep(grid$x, nrow(grid$z)),
      y = rep(grid$y, each = ncol(grid$z)),
      z = as.numeric(grid$z)),
    !is.na(z)
  )
}

plot_irregular_data <- function(x, y, z) {
  ggplot(
    data.table::data.table(x = x, y = y, z = z) %>% unique,
    aes(x = x, y = y, colour = z)
  ) +
    geom_point()
}

with(
  flux_local_al_2_maxima_wide,
  #interp_irregular_data(j_phi_i_i_factor, `t_{inf}_1`, gamma)
  plot_irregular_data(j_phi_i_i_factor, `t_{inf}_1`, gamma)
)

grid <- with(origdata, interp::interp(x, y, z))

griddf <- subset(
  data.frame(
    x = rep(grid$x, nrow(grid$z)),
    y = rep(grid$y, each = ncol(grid$z)),
    z = as.numeric(grid$z)
  ),
  !is.na(z)
)

ggplot(griddf, aes(x, y, z = z)) +
  geom_contour_filled()

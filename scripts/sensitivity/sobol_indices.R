library(sensitivity)
library(ggplot2)
library(cowplot)
#library(patchwork)
library(GGally)

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

# Calculates all indices up to `order` - i.e. if order = 2, it calculates all
# combined 2nd order indices, which can be expensive for many parameters. Seems
# like this can suffer from poor conditioning like `sobol2002`, `sobol2007`,
# despite that it's not mentioned in the docs. Workaround is to centre the
# outputs by subtracting the mean before calling `tell` below.
#x <- sobol(model = NULL, X1 = x_1, X2 = x_2, order = 1, nboot = 100)

# These estimators directly calculate the regular first-order indices and total
# indices simulataneously without calculating the combined second-order indices
# individually. I.e. 2p indices rather than an amount quadratic in p. We'll
# probably want to use one of these (or similar) most of the time.
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

# Sobol indices for cell influx
# -----------------------------
y <- integrated_fluxes$cells_in

# Some of the estimators in `sensitivity` apparently suffer from a
# "conditioning problem", such as `sobol2002` and `sobol2007.` A workaround is
# to centre the model outputs before running estimators or to use alternative
# estimators like `soboljansen`, `sobolmartinez` etc. I'm assuming this is also
# the case for `sobol` (which we initially used), but its docs don't mention
# it. Generally, we'll use on the above alternatives, as they calculate
# first-order and total indices simulataneously, not requiring a number of runs
# quadratic in dimension of parameter space.
#tell(x, y - mean(y))
tell(x, y)
print(x)

# TODO: this simple plotting helper function has a typo ("effet" instead of
# "effect") in the legend - report upstream and/or make our plots from x$S and
# x$T
ggplot(x)

ggsave(
  paste(plot_dir, "sobol_indices_cells_in.pdf", sep = "/"),
  width = 13,
  height = 7
)

# Sobol indices for cell outflux
# ------------------------------
y <- integrated_fluxes$cells_out

tell(x, y)
print(x)
ggplot(x)

ggsave(
  paste(plot_dir, "sobol_indices_cells_out.pdf", sep = "/"),
  width = 13,
  height = 7
)

# Sobol indices for net change in number of cells
# -----------------------------------------------
y <- integrated_fluxes$net_change

tell(x, y)
print(x)
ggplot(x)

# Other plots
# -----------

# Plot the integrated cell numbers against parameter values
ggplot(cells_vs_params_long) +
  geom_point(aes(
    x = param_value,
    y = cells,
    colour = variable,
    group = variable,
  ), size = 2) +
  facet_wrap(vars(param), scales = "free") +
  theme_cowplot() +
  xlab("Parameter value") +
  ylab("Cells (dimensionless)")

ggsave(
  paste(plot_dir, "cells_vs_params.pdf", sep = "/"),
  width = 13,
  height = 7
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

trace_data_longer <- trace_data[
  ,
  .(
    rep, `t_{inf}`, `C_b^{tot}`, `phi_{C_b}^{tot}`, `-F_{phi_i}(x=1)`,
    `-F_{phi_{C_b}}(x=0)`, j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma
  )
  ] %>%
  melt(
    .,
    measure.vars = c(
      "C_b^{tot}", "phi_{C_b}^{tot}", "-F_{phi_i}(x=1)", "-F_{phi_{C_b}}(x=0)"
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

#ggplot(trace_data_longer[variable != "-F_{phi_i}(x=1)"]) +
  #geom_line(
    #aes(x = `t_{inf}`, y = value, colour = param_value, group = rep),
    #alpha = 0.2
  #) +
  #theme_cowplot() +
  #facet_wrap(vars(param, variable), nrow = 4, ncol = 3, scales = "free") +
  #xlab("Time since start of inflammation")

#ggplot(trace_data_longer[variable != "-F_{phi_i}(x=1)"]) +
  #geom_line(
    #aes(x = `t_{inf}`, y = value, colour = param_value, group = rep),
    #alpha = 0.2
  #) +
  #theme_cowplot() +
  #facet_grid(rows = vars(param), cols = vars(variable), scales = "free") +
  #xlab("Time since start of inflammation")

#no_x_axis <- theme(
  #axis.title.x = element_blank(),
  #axis.line.x = element_blank(),
  #axis.text.x = element_blank(),
  #axis.ticks.x = element_blank()
#)

#no_y_axis <- theme(
  #axis.title.y = element_blank(),
  #axis.line.y = element_blank(),
  #axis.text.y = element_blank(),
  #axis.ticks.y = element_blank()
#)

blue <- scale_colour_gradient(low = "black", high = "blue")
red <- scale_colour_gradient(low = "black", high = "red")
green <- scale_colour_gradient(low = "black", high = "green")
orange <- scale_colour_gradient(low = "black", high = "orange")

grid_panel <- function(y_var, colour_by, colour_style, axis_style=NULL) {
  p <- ggplot(
      trace_data,
      aes(
        x = `t_{inf}`,
        y = {{ y_var }},
        group = rep,
        colour = {{ colour_by }}
      )
    ) +
    geom_line(alpha = 0.2) +
    theme_cowplot() +
    colour_style
  if (!is.null(axis_style)) {
    p <- p + axis_style
  }
  p
}

p1 <- grid_panel(`C_b^{tot}`, j_phi_i_i_factor, blue)
p2 <- grid_panel(`phi_{C_b}^{tot}`, j_phi_i_i_factor, blue)
p3 <- grid_panel(`-F_{phi_{C_b}}(x=0)`, j_phi_i_i_factor, blue)

p4 <- grid_panel(`C_b^{tot}`, m_i_factor, red)
p5 <- grid_panel(`phi_{C_b}^{tot}`, m_i_factor, red)
p6 <- grid_panel(`-F_{phi_{C_b}}(x=0)`, m_i_factor, red)

p7 <- grid_panel(`C_b^{tot}`, t_j_phi_i_lag, green)
p8 <- grid_panel(`phi_{C_b}^{tot}`, t_j_phi_i_lag, green)
p9 <- grid_panel(`-F_{phi_{C_b}}(x=0)`, t_j_phi_i_lag, green)

p10 <- grid_panel(`C_b^{tot}`, gamma, orange)
p11 <- grid_panel(`phi_{C_b}^{tot}`, gamma, orange)
p12 <- grid_panel(`-F_{phi_{C_b}}(x=0)`, gamma, orange)

plot_list <- list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12)

# Use the GGally package to make the plot grid - seems to be the only easy way
# that does everything we need, including row/column labels

trace_grid <- ggmatrix(
  plot_list,
  3, 4,
  c("j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag", "gamma"),
  c("C_b^{tot}", "phi_{C_b^{tot}}", "-F_{phi_{C_b}}(x=0)"),
  xlab = "t_{inf}",
  byrow = FALSE
)
trace_grid

# output a png - a pdf of this plot is huge (due to the enormous number of
# lines) and renders very slowly
ggsave(
  plot = trace_grid,
  paste(plot_dir, "trace_coloured_by_params.png", sep = "/"),
  width = 13,
  height = 7
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
  labs(x = "Time since inflammation", y = "First-order Sobol index") +
  coord_cartesian(ylim = c(-0.1, 1.0)) +
  theme_cowplot()

p_flux_first_order

ggsave(
  plot = p_flux_first_order,
  paste(plot_dir, "flux_first_order_sobol_indices.pdf", sep = "/"),
  width = 13,
  height = 7
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

ggsave(
  plot = p_flux_total_order,
  paste(plot_dir, "flux_total_order_sobol_indices.pdf", sep = "/"),
  width = 13,
  height = 7
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

ggsave(
  plot = p_gradient_first_order,
  paste(plot_dir, "gradient_first_order_sobol_indices.pdf", sep = "/"),
  width = 13,
  height = 7
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

ggsave(
  plot = p_gradient_total_order,
  paste(plot_dir, "gradient_total_order_sobol_indices.pdf", sep = "/"),
  width = 13,
  height = 7
)

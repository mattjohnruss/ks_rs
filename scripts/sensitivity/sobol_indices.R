library(conflicted)
library(magrittr)
library(sensitivity)
library(ggplot2)
library(data.table)
library(jsonlite)
library(cowplot)
#library(patchwork)
library(GGally)

# Set the random seed to make the parameter sample reproducible. This will be
# useful for setting off some runs before we've fully decided what quantities
# to extract from the simulations. Once we've decided, we can then re-run this
# script, skipping running the simulations and just calculate the Sobol indices
# etc.
set.seed(12345L)

res_dir_base <- "res_sensitivity_3_neg_pe"

# Constant parameter values
phi_bar_over_c_bar <- 1.40179e-6
phi_bar_over_phi_max <- 0.1
c_bar_over_e <- 0.1
pe <- 5.0
alpha_plus <- 23.25
alpha_minus <- 0.3
beta_plus <- 0.0000385493
beta_minus <- 12.5
n_ccr7 <- 30000.0
q_u <- 0.0
q_b <- 0.0
q_s <- 0.0
d_c_s <- 0.01
d_phi_i <- 0.01
d_phi_m <- 0.01
d_phi_c_u <- 0.01
d_phi_c_b <- 0.01
chi_u <- 0.0
chi_b <- 0.004
chi_s <- 0.0
r <- 0.0
m_h <- 0.00113865
t_inflammation <- 50.0
j_phi_i_h <- 0.000455461
phi_i_init <- 0.9
phi_m_init <- 0.0

# Min/max for the parameters we're varying
j_phi_i_i_factor_min <- 1.0
j_phi_i_i_factor_max <- 1000.0
m_i_factor_min <- 1.0
m_i_factor_max <- 1000.0
t_j_phi_i_lag_min <- 0.0
t_j_phi_i_lag_max <- 25.0
gamma_min <- 0.0
gamma_max <- 10.0

add_constants_and_gammas_to_param_table <- function(params) {
  params[, phi_bar_over_c_bar := phi_bar_over_c_bar]
  params[, phi_bar_over_phi_max := phi_bar_over_phi_max]
  params[, c_bar_over_e := c_bar_over_e]
  params[, pe := pe]
  params[, alpha_plus := alpha_plus]
  params[, alpha_minus := alpha_minus]
  params[, beta_plus := beta_plus]
  params[, beta_minus := beta_minus]
  params[, n_ccr7 := n_ccr7]
  params[, q_u := q_u]
  params[, q_b := q_b]
  params[, q_s := q_s]
  params[, d_c_s := d_c_s]
  params[, d_phi_i := d_phi_i]
  params[, d_phi_m := d_phi_m]
  params[, d_phi_c_u := d_phi_c_u]
  params[, d_phi_c_b := d_phi_c_b]
  params[, chi_u := chi_u]
  params[, chi_b := chi_b]
  params[, chi_s := chi_s]
  params[, r := r]
  params[, m_h := m_h]
  params[, t_inflammation := t_inflammation]
  params[, j_phi_i_h := j_phi_i_h]
  params[, phi_i_init := phi_i_init]
  params[, phi_m_init := phi_m_init]
  params[, gamma_ui := gamma]
  params[, gamma_um := gamma]
  params[, gamma_bi := gamma]
  params[, gamma_bm := gamma]
}

# Write a config file for each set of parameter values to the directory
# corresponding to the index of this run, e.g. `res/32/config.json`. Skips the
# `gamma` column as this has already been duplicated into the various `gamma_*`
# columns that are actually used in the sims.
write_json_params_files <- function(all_params) {
  for (i in 1:nrow(all_params)) {
    j <- jsonlite::toJSON(jsonlite::unbox(all_params[i, !"gamma"]), digits = 16)
    res_dir <- paste(res_dir_base, i - 1, sep = "/")

    if (!dir.exists(res_dir)) {
      dir.create(res_dir, recursive = TRUE)
    }

    config_path <- paste(res_dir, "config.json", sep = "/")
    writeLines(j, config_path)
  }
}

# Run all simulations. First creates the config files and then it calls a shell
# script that runs GNU Parallel with the appropriate arguments. It finally
# reads all simulation outputs from file back into R for processing later.
simulate_all <- function(all_params) {
  n_params <- nrow(all_params)
  cmd <- paste("bash", "scripts/sensitivity/run_all.bash", n_params)
  system(cmd)
}

read_trace_data <- function(all_params) {
  n_runs <- nrow(all_params)

  rep <- 0
  trace_data <- data.table::fread(
    paste(res_dir_base, rep, "trace.csv", sep = "/")
  )
  trace_data <- trace_data[`t_{inf}` >= 0]
  trace_data[, rep := rep]
  trace_data[, output_inf := .I]
  trace_data <- cbind(trace_data, all_params[rep + 1])

  for (rep in 1:(n_runs - 1)) {
    filename <- paste(res_dir_base, rep, "trace.csv", sep = "/")
    print(filename)
    trace_data_new <- data.table::fread(filename)
    trace_data_new <- trace_data_new[`t_{inf}` >= 0]
    trace_data_new[, rep := rep]
    trace_data_new[, output_inf := .I]
    trace_data_new <- cbind(trace_data_new, all_params[rep + 1])
    trace_data <- rbind(trace_data, trace_data_new)
  }

  return(trace_data)
}

calculate_integrated_fluxes <- function(trace_data) {
  # Time-integration of fluxes using the trapezium rule

  # Find the stopping time of the simulation that finishes first - this is where
  # we'll integrate up to
  t_max <- trace_data %>%
    .[, .(t_stop = max(`t_{inf}`)), by = rep] %>%
    .[, min(t_stop)]

  # Calculate quantites for trapezium rule: time intervals between outputs, and
  # pairwise sum of fluxes between consective output times; then computes the
  # integral by constructing the appropriate quotients and summing
  integrated_fluxes <- trace_data[`t_{inf}` <= t_max,
    .(
      influx_sum = frollsum(`-F_{phi_i}(x=1)`, 2),
      outflux_sum = frollsum(`-F_{phi_{C_b}}(x=0)`, 2),
      dt = frollapply(`t_{inf}`, 2, function(x) x[2] - x[1])
    ), by = rep] %>%
    .[, .(q_in = influx_sum / dt, q_out = outflux_sum / dt), by = rep] %>%
    .[, .(
      cells_in = sum(q_in, na.rm = TRUE),
      cells_out = sum(q_out, na.rm = TRUE)
      ), by = rep]

  return(integrated_fluxes)
}

gen_param_sample <- function() {
  data.table(
    j_phi_i_i_factor = runif(
      n_rep,
      min = j_phi_i_i_factor_min, max = j_phi_i_i_factor_max
    ),
    m_i_factor = runif(n_rep, min = m_i_factor_min, max = m_i_factor_max),
    t_j_phi_i_lag = runif(
      n_rep,
      min = t_j_phi_i_lag_min, max = t_j_phi_i_lag_max
    ),
    gamma = runif(n_rep, min = gamma_min, max = gamma_max)
  )
}

n_rep <- 100

x_1 <- gen_param_sample()
x_2 <- gen_param_sample()

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
x <- soboljansen(model = NULL, X1 = x_1, X2 = x_2, nboot = 100)
#x <- sobolmartinez(model = NULL, X1 = x_1, X2 = x_2, nboot = 100)

add_constants_and_gammas_to_param_table(x$X)

#write_json_params_files(x$X)
#fwrite(x$X, paste(res_dir_base, "d_m.csv", sep = "/"), sep = " ")

#simulate_all(x$X)

# TODO: handle the fact that we will sometimes only want to run the simulations
# and not do the sensitivity analysis in the same run, and vice-versa. In the
# latter we'll need to skip writing the parameters to file, skip running the
# sims and just read the combined outputs from the file. Currently
# `run_all.bash` is responsible for running sims and combining outputs, so
# maybe we should split this into two scripts to make it easier to call them as
# needed from R.

trace_data <- read_trace_data(x$X)

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

cells_vs_params <- x$X[
  ,
  .(j_phi_i_i_factor,
    m_i_factor,
    t_j_phi_i_lag,
    gamma,
    cells_in = integrated_fluxes[, cells_in],
    cells_out = integrated_fluxes[, cells_out]
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
    measure.vars = c("cells_in", "cells_out"),
    variable.name = "variable",
    value.name = "cells"
  )

plot_dir <- paste(res_dir_base, "plots", sep = "/")

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

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

# Other plots
# -----------

# Plot the integrated cell numbers against parameter values
ggplot(cells_vs_params_long) +
  geom_point(aes(
    x = param_value,
    y = cells,
    colour = variable,
    group = variable
  )) +
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

no_x_axis <- theme()
#no_x_axis <- theme(
  #axis.title.x = element_blank(),
  #axis.line.x = element_blank(),
  #axis.text.x = element_blank(),
  #axis.ticks.x = element_blank()
#)

no_y_axis <- theme(
  axis.title.y = element_blank(),
  axis.line.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)

blue <- scale_colour_gradient(low = "black", high = "blue")
red <- scale_colour_gradient(low = "black", high = "red")
green <- scale_colour_gradient(low = "black", high = "green")
orange <- scale_colour_gradient(low = "black", high = "orange")

p1 <- ggplot(
    trace_data,
    aes(x = `t_{inf}`, y = `C_b^{tot}`, group = rep, colour = j_phi_i_i_factor)
  ) +
  geom_line(alpha = 0.2) +
  theme_cowplot() +
  blue +
  no_x_axis

p2 <- ggplot(
    trace_data,
    aes(
      x = `t_{inf}`,
      y = `phi_{C_b}^{tot}`,
      group = rep,
      colour = j_phi_i_i_factor
    )
  ) +
  geom_line(alpha = 0.2) +
  theme_cowplot() +
  blue +
  no_x_axis

p3 <- ggplot(
    trace_data,
    aes(
      x = `t_{inf}`,
      y = `-F_{phi_{C_b}}(x=0)`,
      group = rep,
      colour = j_phi_i_i_factor
    )
  ) +
  geom_line(alpha = 0.2) +
  theme_cowplot() +
  blue +
  no_x_axis

p4 <- ggplot(
    trace_data,
    aes(
      x = `t_{inf}`,
      y = `C_b^{tot}`,
      group = rep,
      colour = m_i_factor
      )
    ) +
  geom_line(alpha = 0.2) +
  theme_cowplot() +
  red +
  no_x_axis

p5 <- ggplot(
    trace_data,
    aes(
      x = `t_{inf}`,
      y = `phi_{C_b}^{tot}`,
      group = rep,
      colour = m_i_factor
      )
    ) +
  geom_line(alpha = 0.2) +
  theme_cowplot() +
  red +
  no_x_axis

p6 <- ggplot(
    trace_data,
    aes(
      x = `t_{inf}`,
      y = `-F_{phi_{C_b}}(x=0)`,
      group = rep,
      colour = m_i_factor
      )
    ) +
  geom_line(alpha = 0.2) +
  theme_cowplot() +
  red +
  no_x_axis

p7 <- ggplot(
    trace_data,
    aes(
      x = `t_{inf}`,
      y = `C_b^{tot}`,
      group = rep,
      colour = t_j_phi_i_lag
      )
    ) +
  geom_line(alpha = 0.2) +
  theme_cowplot() +
  green +
  no_x_axis

p8 <- ggplot(
    trace_data,
    aes(
      x = `t_{inf}`,
      y = `phi_{C_b}^{tot}`,
      group = rep,
      colour = t_j_phi_i_lag
      )
    ) +
  geom_line(alpha = 0.2) +
  theme_cowplot() +
  green +
  no_x_axis

p9 <- ggplot(
    trace_data,
    aes(
      x = `t_{inf}`,
      y = `-F_{phi_{C_b}}(x=0)`,
      group = rep,
      colour = t_j_phi_i_lag
      )
    ) +
  geom_line(alpha = 0.2) +
  theme_cowplot() +
  green +
  no_x_axis

p10 <- ggplot(
    trace_data,
    aes(
      x = `t_{inf}`,
      y = `C_b^{tot}`,
      group = rep,
      colour = gamma
      )
    ) +
  geom_line(alpha = 0.2) +
  theme_cowplot() +
  orange

p11 <- ggplot(
    trace_data,
    aes(
      x = `t_{inf}`,
      y = `phi_{C_b}^{tot}`,
      group = rep,
      colour = gamma
      )
    ) +
  geom_line(alpha = 0.2) +
  theme_cowplot() +
  orange

p12 <- ggplot(
    trace_data,
    aes(
      x = `t_{inf}`,
      y = `-F_{phi_{C_b}}(x=0)`,
      group = rep,
      colour = gamma
      )
    ) +
  geom_line(alpha = 0.2) +
  theme_cowplot() +
  orange

#p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + p11 + p12 +
  #plot_layout(guides = "collect", ncol = 3)

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

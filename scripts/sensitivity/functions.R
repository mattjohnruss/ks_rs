library(data.table)
library(magrittr)
library(jsonlite)
library(stringr)

# the ratio of the timescales used in the 1D and cleaving models - used to convert parameters inferred from cleaving data into the appropriate scale for the 1D model
t_1d_over_t_cleaving <- 2500.0 / 3600.0

# Constant parameter values
const_params <- list(
  phi_bar_over_c_bar = 1.40179e-6,
  phi_bar_over_phi_max = 0.1,
  c_bar_over_e = 0.1,
  # IDEA: get rid of pe from here and do e.g. `x$X[, pe := -5]` in the script,
  pe = -5.0,
  alpha_plus = 23.25,
  alpha_minus = 0.3,
  beta_plus = 0.0000385493,
  beta_minus = 12.5,
  n_ccr7 = 30000.0,
  n = 2.0,
  a_bar = 5.566645e+00 * t_1d_over_t_cleaving,
  gamma = 4.140805e+00 * t_1d_over_t_cleaving,
  q_u = 0.0,
  q_b = 0.0,
  q_s = 0.0,
  d_f = 2.793900e-02 * t_1d_over_t_cleaving,
  d_t = 9.196955e-03 * t_1d_over_t_cleaving,
  d_c_s = 10.0,
  d_phi_i = 0.01,
  d_phi_m = 0.01,
  d_phi_c_u = 0.01,
  d_phi_c_b = 0.01,
  d_phi_c_s = 0.01,
  d_j = 1.0,
  chi_u = 0.004,
  chi_b = 0.004,
  chi_s = 0.004,
  mu_m = 2.532745e-05 * t_1d_over_t_cleaving,
  m_h = 0.00113865,
  t_inflammation = 50.0,
  j_phi_i_h = 0.000455461,
  phi_i_init = 0.9,
  phi_m_init = 0.0
)

add_constants_to_param_table <- function(params, const_params) {
  params[, names(const_params) := const_params]
}

# Parameter and variable labels as `expression`s for nicer formatting in plots
param_labels <- c(
  "j_phi_i_i_factor" = expression(paste(J[phi[i]]^I, " factor")),
  "m_i_factor" = expression(paste(M^I, " factor")),
  "t_j_phi_i_lag" = expression(paste(J[phi[i]]^I, " delay")),
  "gamma" = expression(gamma),
  "pe" = expression(Pe)
)

param_labels_words <- c(
  "j_phi_i_i_factor" = "Ingress ratio",
  "m_i_factor" = "Maturation ratio",
  "t_j_phi_i_lag" = "Lag",
  "gamma" = "Cleavage\n rate",
  "pe" = expression(Pe)
)

param_labels_words_no_breaks <- c(
  "j_phi_i_i_factor" = "Ingress ratio",
  "m_i_factor" = "Maturation ratio",
  "t_j_phi_i_lag" = "Lag",
  "gamma" = "Cleavage rate",
  "pe" = expression(Pe)
)

# Add any other variables as needed
variable_labels <- c(
  "C_u" = expression(C[u]),
  "C_b" = expression(C[b]),
  "C_s" = expression(C[s]),
  "phi_i" = expression(phi[i]),
  "phi_m" = expression(phi[m]),
  "phi_C_u" = expression(phi[C[u]]),
  "phi_C_b" = expression(phi[C[b]]),
  "phi_C_s" = expression(phi[C[s]]),
  "C_u^{tot}" = expression(paste("Total ", C[u])),
  "C_b^{tot}" = expression(paste("Total ", C[b])),
  "C_s^{tot}" = expression(paste("Total ", C[s])),
  "phi_{C_u}^{tot}" = expression(paste("Total ", phi[C[u]])),
  "phi_{C_b}^{tot}" = expression(paste("Total ", phi[C[b]])),
  "phi_{C_s}^{tot}" = expression(paste("Total ", phi[C[s]])),
  "J^{tot}" = expression(paste("Total ", J)),
  "cell_outflux" = expression(paste("Cell LV flux"))
)

all_labels <- c(param_labels, variable_labels)

# Min/max for the parameters we're varying
j_phi_i_i_factor_min <- 1.0
j_phi_i_i_factor_max <- 1000.0
m_i_factor_min <- 1.0
m_i_factor_max <- 1000.0
t_j_phi_i_lag_min <- 0.0
t_j_phi_i_lag_max <- 25.0

param_names <- c(
  "j_phi_i_i_factor",
  "m_i_factor",
  "t_j_phi_i_lag"
)
param_mins <- c(
  j_phi_i_i_factor_min,
  m_i_factor_min,
  t_j_phi_i_lag_min
)
param_maxs <- c(
  j_phi_i_i_factor_max,
  m_i_factor_max,
  t_j_phi_i_lag_max
)

# Write a config file for each set of parameter values to the directory
# corresponding to the index of this run, e.g. `res/32/config.json`.
write_json_params_files <- function(all_params, res_dir_base) {
  n_runs <- nrow(all_params)
  for (i in 1:n_runs) {
    j <- jsonlite::toJSON(jsonlite::unbox(all_params[i]), digits = 16)
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
simulate_all <- function(all_params, res_dir_base) {
  n_runs <- nrow(all_params)
  cmd <- paste(
    "bash",
    "scripts/sensitivity/run_all.bash",
    n_runs,
    res_dir_base
  )
  system(cmd)
}

read_trace_data <- function(all_params, res_dir_base) {
  n_runs <- nrow(all_params)
  trace_data <- list()

  for (rep in 1:n_runs) {
    filename <- paste(res_dir_base, rep - 1, "trace.csv", sep = "/")
    print(filename)
    trace_data[[rep]] <- data.table::fread(filename)
    trace_data[[rep]][, rep := rep - 1]
    trace_data[[rep]][`t_{inf}` >= 0, output_inf := .I]
    trace_data[[rep]] <- cbind(trace_data[[rep]], all_params[rep])
  }

  return(rbindlist(trace_data))
}

# Reads spatial data for all reps at a particular time
read_spatial_data <- function(all_params, time, res_dir_base) {
  n_runs <- nrow(all_params)
  spatial_data <- list()

  for (rep in 1:n_runs) {
    filename <- paste(
      res_dir_base,
      rep - 1,
      sprintf("inflammation/output_%05i.csv", time),
      sep = "/"
    )
    if (file.exists(filename)) {
      #print(filename)
      spatial_data[[rep]] <- data.table::fread(
        filename,
        blank.lines.skip = TRUE
      )
      spatial_data[[rep]][, rep := rep - 1]
      spatial_data[[rep]][, output_inf := time]
      spatial_data[[rep]] <- cbind(spatial_data[[rep]], all_params[rep])
    } else {
      print(paste0("skipping missing file ", filename))
    }
  }
  spatial_data <- rbindlist(spatial_data)
  setnames(spatial_data, str_remove_all(names(spatial_data), "[\\\\${}]"))
  return(spatial_data)
}

# Time-integration of fluxes using the trapezium rule
calculate_integrated_fluxes <- function(trace_data) {
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
      outflux_sum = frollsum(`cell_outflux`, 2),
      dt = frollapply(`t_{inf}`, 2, function(x) x[2] - x[1])
    ), by = rep] %>%
    .[, .(
      q_in = influx_sum * 0.5 * dt,
      q_out = outflux_sum * 0.5 * dt
      ), by = rep] %>%
    .[, .(
      cells_in = sum(q_in, na.rm = TRUE),
      cells_out = sum(q_out, na.rm = TRUE)
      ), by = rep]

  return(integrated_fluxes)
}

# Uses independent uniform samples for each param
gen_param_sample_unif <- function(n_rep, names, mins, maxs) {
  d <- list()
  for (i in seq_len(length(names))) {
    d[[i]] <- data.table(
      id = 1:n_rep,
      param = rep(names[i], n_rep),
      value = runif(n_rep, mins[i], maxs[i])
    )
  }
  d <- d %>%
    data.table::rbindlist(.) %>%
    dcast(., id ~ param) %>%
    .[, id := NULL]
  data.table::setcolorder(d, names)
  return(d)
}

# Uses Sobol sequence for all params together
gen_param_sample_sobol <- function(n_rep, min_max, init = TRUE) {
  d <- randtoolbox::sobol(n_rep, nrow(min_max), init) %>%
    as.data.table
  names(d) <- min_max$param
  d <- melt(d, measure.vars = min_max$param, variable.name = "param") %>%
    .[min_max, on = "param"] %>%
    .[, value := min + value * (max - min)] %>%
    .[, min := NULL] %>%
    .[, max := NULL] %>%
    .[, id := seq_len(.N), by = param] %>%
    dcast(., id ~ param, value.var = "value") %>%
    .[, id := NULL]
}

ggsave_with_defaults <- function(
  filename,
  width = 10,
  height = 7,
  device = NULL,
  bg = NULL,
  ...) {
  if (tools::file_ext(filename) == "pdf") {
    device <- cairo_pdf
  }

  if (tools::file_ext(filename) == "png") {
    bg <- "white"
  }

  ggplot2::ggsave(
    filename,
    width = width,
    height = height,
    device = device,
    bg = bg,
    ...
  )
}

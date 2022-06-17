library(data.table)
library(magrittr)
library(jsonlite)

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
write_json_params_files <- function(all_params, res_dir_base) {
  n_runs <- nrow(all_params)
  for (i in 1:n_runs) {
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
      outflux_sum = frollsum(`-F_{phi_{C_b}}(x=0)`, 2),
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

# This works but I suspect it could be written in a nicer way
gen_param_sample <- function(n_rep, names, mins, maxs) {
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

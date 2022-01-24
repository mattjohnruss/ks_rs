library(conflicted)
library(magrittr)
library(sensitivity)
library(ggplot2)
library(data.table)
library(jsonlite)

# Set the random seed to make the parameter sample reproducible. This will be
# useful for setting off some runs before we've fully decided what quantities
# to extract from the simulations. Once we've decided, we can then re-run this
# script, skipping running the simulations and just calculate the Sobol indices
# etc.
set.seed(12345L)

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
    j <- toJSON(unbox(all_params[i, !"gamma"]), digits = 16)
    res_dir <- paste0("res_sensitivity/", i - 1)

    if (!dir.exists(res_dir)) {
      dir.create(res_dir, recursive = TRUE)
    }

    config_path <- paste0(res_dir, "/config.json")
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

read_combined_output <- function(all_params) {
  data <- all_params

  results <- fread("res_sensitivity/all_outputs.csv")
  #print(results)

  print(system.time(data[, value := results[, '-F_{phi_{C_b}}(x=0)']]))
}

gen_param_sample <- function() {
  data.table(
    j_phi_i_i_factor = runif(n_rep, min = j_phi_i_i_factor_min, max = j_phi_i_i_factor_max),
    m_i_factor = runif(n_rep, min = m_i_factor_min, max = m_i_factor_max),
    t_j_phi_i_lag = runif(n_rep, min = t_j_phi_i_lag_min, max = t_j_phi_i_lag_max),
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

#simulate_all(x$X)

# TODO: handle the fact that we will sometimes only want to run the simulations
# and not do the sensitivity analysis in the same run, and vice-versa. In the
# latter we'll need to skip writing the parameters to file, skip running the
# sims and just read the combined outputs from the file. Currently
# `run_all.bash` is responsible for running sims and combining outputs, so
# maybe we should split this into two scripts to make it easier to call them as
# needed from R.

# TODO: reenable/make flexible - see above
read_combined_output(x$X)

y <- x$X$value

# TODO: Why on earth does subtracting the mean from y here make the Sobol
# indices much more believable?
# Possible answer: Some of the estimators in `sensitivity` apparently suffer
# from a "conditioning problem", such as `sobol2002` and `sobol2007.` A
# workaround is to centre the model outputs before running estimators or to use
# alternative estimators like `soboljansen`, `sobolmartinez` etc. I'm assuming
# this is also the case for `sobol` (which we initially used), but its docs
# don't mention it.
tell(x, y - mean(y))

print(x)
plot(x)

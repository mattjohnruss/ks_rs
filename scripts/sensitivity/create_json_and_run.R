library(sensobol)

source("scripts/sensitivity/functions.R")

# Set the random seed to make the parameter sample reproducible
set.seed(12345L)

res_dir_base <- "res_sensitivity_1000_pe_5"

n_param_sample <- 1000

# Create the parameter sample and Sobol matrices
param_sample <- sobol_matrices(
  N = n_param_sample,
  params = param_names
)

# Rescale the values from U(0, 1) -> U(min, max)
for (i in seq_along(param_names)) {
  min <- param_mins[i]
  max <- param_maxs[i]
  name <- param_names[i]
  param_sample[, name] <- qunif(param_sample[, name], min, max)
}

param_sample_dt <- data.table(param_sample)

add_constants_and_gammas_to_param_table(param_sample_dt, const_params)

write_json_params_files(param_sample_dt, res_dir_base)

# Save the parameter sample (matrix with just the varied params, and full
# data.table)
saveRDS(
  list(
    n_param_sample = n_param_sample,
    param_sample = param_sample,
    param_sample_dt = param_sample_dt
  ),
  paste(res_dir_base, "param_sample.rds", sep = "/")
)

simulate_all(param_sample_dt, res_dir_base)

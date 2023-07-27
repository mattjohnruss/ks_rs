library(sensitivity)

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

#x_1 <- gen_param_sample_unif(100, names, mins, maxs)
#x_2 <- gen_param_sample_unif(100, names, mins, maxs)

param_min_max <- data.table(param = names, min = mins, max = maxs)
x_1 <- gen_param_sample_sobol(100, param_min_max)
x_2 <- gen_param_sample_sobol(100, param_min_max, init = FALSE)

#x <- soboljansen(model = NULL, X1 = x_1, X2 = x_2, nboot = 100)
x <- sobolmartinez(model = NULL, X1 = x_1, X2 = x_2, nboot = 100)

add_constants_and_gammas_to_param_table(x$X, const_params)

write_json_params_files(x$X, res_dir_base)

# Save the parameter sample
fwrite(x$X1, paste(res_dir_base, "x_1.csv", sep = "/"), sep = " ")
fwrite(x$X2, paste(res_dir_base, "x_2.csv", sep = "/"), sep = " ")
fwrite(x$X, paste(res_dir_base, "d_m.csv", sep = "/"), sep = " ")

#simulate_all(x$X, res_dir_base)

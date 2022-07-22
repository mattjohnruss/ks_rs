library(data.table)
library(magrittr)
library(stringr)

params <- fread(paste(res_dir_base, "d_m.csv", sep = "/"))

n_rep <- 600

get_max_single <- function(rep_id) {
  rep <- rep_id - 1

  filenames <- list.files(
    paste(res_dir_base, rep, "inflammation", sep = "/"),
    pattern = "output_[0-9]{5}\\.csv",
    full.names = TRUE
  ) %>% sort()

  n_files <- length(filenames)

  data <- list()

  for (i in 1:n_files) {
    filename <- sprintf(
      paste(res_dir_base, rep, "inflammation", "output_%05d.csv", sep = "/"),
      i - 1
    )
    data[[i]] <- fread(filename, blank.lines.skip = TRUE)
    data[[i]][, rep := rep]
    data[[i]][, output_inf := i]
  }

  data <- rbindlist(data)

  # Remove all the awkward characters from the columns names
  setnames(data, str_remove_all(names(data), "[\\\\${}]"))

  # When there are multiple equal maxima we just take the first one. They are
  # likely to be close to each other, and it only matters when considering the
  # maximum as a function of space.

  # Find the maximum of phi_C_b at each time and save all variables at that
  # time
  phi_c_b_data <- data[data[, .I[phi_C_b == max(phi_C_b)][1], by = time_inf]$V1]

  # add cell indices
  data[, cell := (1:.N - 1) %/% 2, by = time_inf]

  # calcuate the gradient of C_b within each cell
  dc_b_dx <- data[, .(dC_b_dx = diff(C_b) / diff(x)), by = .(time_inf, cell)]

  # calculate the means of all other quantities within each cell, and add the
  # gradient to this new table
  data <- data[, lapply(.SD, mean), by = .(time_inf, cell)]
  data <- cbind(data, dc_b_dx[, .(dC_b_dx)])

  # find the maximum gradient of C_b within each (computational) cell at each
  # time and save values of all variables at that time
  dc_b_dx_data <- data[
    data[, .I[abs(dC_b_dx) == max(abs(dC_b_dx))][1], by = time_inf]$V1
  ]

  list(phi_c_b_data, dc_b_dx_data)
}

#library(parallel)

#options(mc.cores = 16)

#system.time(
  #all_max_data <- mclapply(1:n_rep, get_max_single)
#)

#max_phi_c_b_all <- lapply(1:n_rep, function(i) all_max_data[[i]][[1]])
#max_dc_b_dx_all <- lapply(1:n_rep, function(i) all_max_data[[i]][[2]])

#max_phi_c_b_all <- rbindlist(max_phi_c_b_all)
#max_dc_b_dx_all <- rbindlist(max_dc_b_dx_all)

#fwrite(
  #max_phi_c_b_all,
  #paste(res_dir_base, "max_phi_c_b_all.csv", sep = "/"),
  #sep = " "
#)
#fwrite(
  #max_dc_b_dx_all,
  #paste(res_dir_base, "max_dc_b_dx_all.csv", sep = "/"),
  #sep = " "
#)

max_phi_c_b_all <- fread(paste(res_dir_base, "max_phi_c_b_all.csv", sep = "/"))
max_dc_b_dx_all <- fread(paste(res_dir_base, "max_dc_b_dx_all.csv", sep = "/"))

max_phi_c_b_all_and_params <- cbind(
  max_phi_c_b_all,
  params[
    max_phi_c_b_all[, rep + 1],
    .(j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma)
  ]
)

max_dc_b_dx_all_and_params <- cbind(
  max_dc_b_dx_all,
  params[
    max_dc_b_dx_all[, rep + 1],
    .(j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma)
  ]
)

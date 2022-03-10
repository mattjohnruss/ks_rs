library(data.table)
library(magrittr)
library(ggplot2)
library(cowplot)
library(stringr)

params <- fread("d_m.csv")

n_rep <- 600

get_max_single <- function(rep_id) {
  rep <- rep_id - 1

  filenames <- list.files(
    paste(rep, "inflammation", sep = "/"),
    pattern = "output_[0-9]{5}\\.csv",
    full.names = TRUE
  ) %>% sort()

  n_files <- length(filenames)

  data <- fread(
    paste(rep, "inflammation", "output_00000.csv", sep = "/"),
    blank.lines.skip = TRUE
  )
  data[, rep := rep]
  data[, file := 0]

  for (i in 1:(n_files - 1)) {
    filename <- sprintf(
      paste(rep, "inflammation", "output_%05d.csv", sep = "/"),
      i
    )
    data_new <- fread(filename, blank.lines.skip = TRUE)
    data_new[, rep := rep]
    data_new[, file := i]
    data <- rbind(data, data_new)
  }

  # Remove all the awkward characters from the columns names
  setnames(data, str_remove_all(names(data), "[\\\\${}]"))

  # Find the maximum of phi_C_b at each time and return all variables at that
  # time
  data[data[, .I[phi_C_b == max(phi_C_b)], by = time_inf]$V1]
}

#max_phi_c_b_all <- list()
#for (rep_id in 1:n_rep) {
  #max_phi_c_b_all[[rep_id]] <- get_max_single(rep_id)
#}
max_phi_c_b_all <- rbindlist(max_phi_c_b_all)

#fwrite(max_phi_c_b_all, "max_phi_c_b_all.csv", sep = " ")

max_phi_c_b_all <- fread("max_phi_c_b_all.csv")

max_phi_c_b_all_and_params <- cbind(
  max_phi_c_b_all,
  params[
    max_phi_c_b_all[, rep + 1],
    .(j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma)
  ]
)

ggplot(
  max_phi_c_b_all_and_params,
  aes(time_inf, x, group = rep, colour = gamma)
) +
  geom_path(alpha = 0.1) +
  theme_cowplot()

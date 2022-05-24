library(data.table)
library(magrittr)
library(ggplot2)
library(cowplot)
library(stringr)
library(patchwork)
library(parallel)

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

# Plots
#------

# Some plotting styles and helper functions

blue <- scale_colour_gradient(
  low = "black",
  high = "blue",
  guide = guide_colorbar(title = expression(paste(J[phi[i]]^I, " factor")))
)
red <- scale_colour_gradient(
  low = "black",
  high = "red",
  guide = guide_colorbar(title = expression(paste(M^I, " factor")))
)
green <- scale_colour_gradient(
  low = "black",
  high = "green",
  guide = guide_colorbar(title = expression(paste(J[phi[i]]^I, " delay")))
)
orange <- scale_colour_gradient(
  low = "black",
  high = "orange",
  guide = guide_colorbar(
    title = expression(paste(gamma[u], ", ", gamma[u], " ..."))
  )
)

p_location_panel <- function(data, colour_by, colour_style) {
  ggplot(
    data,
    aes(time_inf, x, group = rep, colour = {{ colour_by }})
    ) +
  geom_path(alpha = 0.1) +
  coord_cartesian(ylim = c(0, 1)) +
  colour_style +
  xlab("Time since inflammation") +
  ylab(expression(x))
}

p_value_panel <- function(data, colour_by, colour_style, ylabel) {
  ggplot(
    data,
    aes(time_inf, phi_C_b, group = rep, colour = {{ colour_by }})
    ) +
  geom_path(alpha = 0.1) +
  colour_style +
  xlab("Time since inflammation") +
  ylab(ylabel)
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

ggsave(
  plot = p_phi_c_b_location,
  paste(plot_dir, "max_phi_c_b_location.png", sep = "/"),
  width = 13,
  height = 7
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

p_phi_c_b_value

ggsave(
  plot = p_phi_c_b_value,
  paste(plot_dir, "max_phi_c_b_value.png", sep = "/"),
  width = 13,
  height = 7
)

# Location of max gradient of c_b

p_j_phi_i_i_factor <- p_location_panel(max_dc_b_dx_all_and_params, j_phi_i_i_factor, blue)
p_m_i_factor <- p_location_panel(max_dc_b_dx_all_and_params, m_i_factor, red)
p_t_j_phi_i_lag <- p_location_panel(max_dc_b_dx_all_and_params, t_j_phi_i_lag, green)
p_gamma <- p_location_panel(max_dc_b_dx_all_and_params, gamma, orange)

p_dc_b_dx_location <- p_j_phi_i_i_factor + p_m_i_factor + p_t_j_phi_i_lag + p_gamma +
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

p_dc_b_dx_location

ggsave(
  plot = p_dc_b_dx_location,
  paste(plot_dir, "max_dc_b_dx_location.png", sep = "/"),
  width = 13,
  height = 7
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

p_dc_b_dx_value

ggsave(
  plot = p_dc_b_dx_value,
  paste(plot_dir, "max_dc_b_dx_value.png", sep = "/"),
  width = 13,
  height = 7
)

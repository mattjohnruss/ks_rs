library(data.table)
library(magrittr)
library(ggplot2)
library(cowplot)
library(stringr)
library(patchwork)

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

#max_phi_c_b_all <- rbindlist(max_phi_c_b_all)

#fwrite(max_phi_c_b_all, "max_phi_c_b_all.csv", sep = " ")

max_phi_c_b_all <- fread("max_phi_c_b_all.csv")

max_phi_c_b_all_and_params <- cbind(
  max_phi_c_b_all,
  params[
    max_phi_c_b_all[, rep + 1],
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
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  plot = p_phi_c_b_location,
  "plots/max_phi_c_b_location.png",
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
  theme(plot.title = element_text(hjust = 0.5))

p_phi_c_b_value

ggsave(
  plot = p_phi_c_b_value,
  "plots/max_phi_c_b_value.png",
  width = 13,
  height = 7
)

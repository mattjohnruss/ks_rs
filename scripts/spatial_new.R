library(sensobol)
library(ggplot2)
library(cowplot)
library(patchwork)
library(GGally)
library(ggnewscale)
library(ggh4x)
library(ggbreak)

source("scripts/sensitivity/functions.R")

res_dir_base <- "res_sensitivity_100_pe_-5_inhib"

# Set and create the plot directory if it doesn't exist
plot_dir <- paste(res_dir_base, "plots/spatial", sep = "/")

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# Read the param_sample object from disk. This contains the initial (Sobol
# sequence or uniformly random) parameter samples and the resulting "design
# matrix" (`x$X`)
x <- readRDS(paste(res_dir_base, "param_sample.rds", sep = "/"))
# The size of the original parameter samples, not the total number of runs
n_param_sample <- x$n_param_sample
param_sample <- x$param_sample
param_sample_dt <- x$param_sample_dt

#time <- 100

#variable_names <- c("C_u", "C_b", "C_s", "phi_i", "phi_m", "phi_C_u", "phi_C_b", "phi_C_s", "J")

# find the max of each variable across all repeats
#maxes <- list()
#i <- 1
#for (time in c(seq(0, 100, 10))) {
  #print(time)
  #spatial_data <- read_spatial_data(param_sample_dt, time, res_dir_base)
  #maxes[[i]] <- spatial_data[, lapply(.SD, max), .SDcols = variable_names]
  #i <- i + 1
#}

#maxes <- rbindlist(maxes)
#maxes <- maxes[, lapply(.SD, max), .SDcols = variable_names]

# plot spatial profiles for a range of times and colour by the 3 params we are varying
for (time in c(seq(0, 1000, 10))) {
  print(time)

  spatial_data <- read_spatial_data(param_sample_dt, time, res_dir_base)

  spatial_data_long <- spatial_data |>
    melt(measure.vars = c("C_u", "C_b", "C_s", "phi_i", "phi_m", "phi_C_u", "phi_C_b", "phi_C_s", "J"))

  #param <- "j_phi_i_i_factor"

  for (param in c("j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag")) {
    print(param)
    p_spatial <- ggplot(
      spatial_data_long,
      aes(x = x, y = value, colour = .data[[param]], group = rep)
    ) +
      geom_line(linewidth = 0.2, alpha = 0.2) +
      facet_wrap(
        vars(variable),
        scales = "free",
        labeller = as_labeller(function(v) variable_labels[v], label_parsed)
      ) +
      labs(
        title = paste0("Time: ", time),
        x = "Spatial distance",
        y = NULL,
        colour = param_labels_words[param]
      ) +
      colour_scales[param] +
      theme(legend.position = "top", legend.justification = "centre") +
      # for pe = -5
      #ggh4x::facetted_pos_scales(y = list(
        #variable == "C_u" ~ scale_y_continuous(limits = c(0, 1)),
        #variable == "C_b" ~ scale_y_continuous(limits = c(0, 9)),
        #variable == "C_s" ~ scale_y_continuous(limits = c(0, 0.03)),
        #variable == "phi_i" ~ scale_y_continuous(limits = c(0, 8)),
        #variable == "phi_m" ~ scale_y_continuous(limits = c(0, 7)),
        #variable == "phi_C_u" ~ scale_y_continuous(limits = c(0, 0.055)),
        #variable == "phi_C_b" ~ scale_y_continuous(limits = c(0, 1.3)),
        #variable == "phi_C_s" ~ scale_y_continuous(limits = c(0, 0.07)),
        #variable == "J" ~ scale_y_continuous(limits = c(0, 0.9))
      #))
      # for pe = 5
      #ggh4x::facetted_pos_scales(y = list(
        #variable == "C_u" ~ scale_y_continuous(limits = c(0, 1)),
        #variable == "C_b" ~ scale_y_continuous(limits = c(0, 9)),
        #variable == "C_s" ~ scale_y_continuous(limits = c(0, 0.1)),
        #variable == "phi_i" ~ scale_y_continuous(limits = c(0, 8.5)),
        #variable == "phi_m" ~ scale_y_continuous(limits = c(0, 4.1)),
        #variable == "phi_C_u" ~ scale_y_continuous(limits = c(0, 0.29)),
        #variable == "phi_C_b" ~ scale_y_continuous(limits = c(0, 3.9)),
        #variable == "phi_C_s" ~ scale_y_continuous(limits = c(0, 0.25)),
        #variable == "J" ~ scale_y_continuous(limits = c(0, 0.35))
      #))
    p_spatial

    name <- sprintf("spatial_t=%05d_colour=%s.png", time, param)

    ggsave_with_defaults(
      plot = p_spatial,
      paste(plot_dir, name, sep = "/")
    )
  }
}

## manually plot
#ggplot(
  #spatial_data_long,
  #aes(x = x, y = value, colour = .data[[param]], group = rep)
#) +
  #geom_line(alpha = 0.2) +
  #facet_wrap(vars(variable), scales = "free") +
  #colour_scales[param] +
  #theme_cowplot()

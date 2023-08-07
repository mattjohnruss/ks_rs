library(sensobol)
library(ggplot2)
library(cowplot)
library(patchwork)
library(GGally)

source("scripts/sensitivity/functions.R")

#theme_set(theme_cowplot(font_size = 17))
theme_set(theme_cowplot())

sobol_neg_and_pos <- function() {

  # Functions
  # ---------
  load_trace_data <- function(pos = FALSE) {
    # Set the seed manually so multiple calls generate the same param samples
    set.seed(12345L)

    # Set the (global, hence <<-) pe variable, and the data dir, to the
    # appropriate value
    if (pos) {
      pe <<- 5.0
      res_dir_base <- "res_sensitivity_1000_pe_5_traces_only"
    } else {
      pe <<- -5.0
      res_dir_base <- "res_sensitivity_1000_pe_-5_traces_only"
    }

    # Min/max for the parameters we're varying
    j_phi_i_i_factor_min <- 1.0
    j_phi_i_i_factor_max <- 1000.0
    m_i_factor_min <- 1.0
    m_i_factor_max <- 1000.0
    t_j_phi_i_lag_min <- 0.0
    t_j_phi_i_lag_max <- 25.0
    gamma_min <- 0.0
    gamma_max <- 10.0

    names <- c("j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag", "gamma")
    mins <- c(
      j_phi_i_i_factor_min,
      m_i_factor_min,
      t_j_phi_i_lag_min,
      gamma_min
    )
    maxs <- c(j_phi_i_i_factor_max,
      m_i_factor_max,
      t_j_phi_i_lag_max,
      gamma_max
    )

    x <- readRDS(paste(res_dir_base, "param_sample.rds", sep = "/"))
    n_param_sample <- x$n_param_sample
    param_sample <- x$param_sample
    param_sample_dt <- x$param_sample_dt

    trace_data_full <- read_trace_data(param_sample_dt, res_dir_base)
    trace_data <- trace_data_full[`t_{inf}` >= 0]

    list(
      n_param_sample = n_param_sample,
      param_sample = param_sample,
      param_sample_dt = param_sample_dt,
      trace_data = trace_data,
      param_names = names
    )
  }

  sobol_indices_cells_panel <- function(
    ind_neg,
    ind_pos,
    title,
    first_plot = TRUE
  ) {
    ind_neg$results[, pe := "neg"]
    ind_pos$results[, pe := "pos"]
    data <- rbind(ind_neg$results, ind_pos$results)
    data[, factor(
      parameters,
      levels = c("j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag", "gamma")
    )]

    y_label <- if (first_plot) {
      "Sobol index"
    } else {
      NULL
    }

    p <- ggplot(
      data,
      aes(
        x = parameters,
        y = original,
        ymin = low.ci,
        ymax = high.ci,
        shape = sensitivity,
        colour = pe
      )
    ) +
      geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
      geom_errorbar(width = 0.5, position = position_dodge(width = 0.6)) +
      scale_x_discrete(labels = param_labels_words_no_breaks) +
      scale_shape_discrete(
        labels = c("Si" = "First-order", "Ti" = "Total-order")
      ) +
      scale_colour_discrete(
        labels = c("pos" = expression(Pe == 5), "neg" = expression(Pe == -5))
      ) +
      coord_cartesian(ylim = c(0.0, 1.0)) +
      labs(
        x = NULL,
        y = y_label,
        shape = "Sobol index",
        colour = "Flow",
        title = title
      ) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.text.align = 0,
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
      )

    if (!first_plot) {
      p <- p + theme(
        #axis.line.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
      )
    }

    p
  }

  sobol_indices_cells_panel_2 <- function(
    ind_neg,
    ind_pos
  ) {
    ind_neg$results[, pe := "neg"]
    ind_pos$results[, pe := "pos"]
    data <- rbind(ind_neg$results, ind_pos$results)
    data[, factor(
      parameters,
      levels = c("j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag", "gamma")
    )]

    y_label <- "Sobol index"

    p <- ggplot(
      data,
      aes(
        x = parameters,
        y = original,
        ymin = low.ci,
        ymax = high.ci,
        shape = sensitivity,
        colour = pe
      )
    ) +
      geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
      geom_errorbar(width = 0.5, position = position_dodge(width = 0.6)) +
      scale_x_discrete(labels = param_labels_words_no_breaks) +
      scale_shape_discrete(
        labels = c("Si" = "First-order", "Ti" = "Total-order")
      ) +
      scale_colour_discrete(
        labels = c("pos" = expression(Pe == 5), "neg" = expression(Pe == -5))
      ) +
      coord_cartesian(ylim = c(0.0, 1.2)) +
      labs(
        x = NULL,
        y = y_label,
        shape = "Sobol index",
        colour = "Flow"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.text.align = 0,
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
      )

    p
  }

  cells_vs_params_long <- function(param_sample_dt, integrated_fluxes) {
    cells_vs_params <- param_sample_dt[
      ,
      .(j_phi_i_i_factor,
        m_i_factor,
        t_j_phi_i_lag,
        gamma,
        cells_in = integrated_fluxes[, cells_in],
        cells_out = integrated_fluxes[, cells_out],
        net_change = integrated_fluxes[, net_change]
      )
    ]

    cells_vs_params %>%
      melt(
        .,
        measure.vars = c(
          "j_phi_i_i_factor",
          "m_i_factor",
          "t_j_phi_i_lag",
          "gamma"
        ),
        variable.name = "param",
        value.name = "param_value"
      ) %>%
      melt(
        .,
        measure.vars = c("cells_in", "cells_out", "net_change"),
        variable.name = "variable",
        value.name = "cells"
      )
  }

  cells_vs_params_panel <- function(cells_vs_params, variable_str, param_str) {
    p <- ggplot(cells_vs_params[variable == variable_str & param == param_str]) +
      geom_point(aes(
        x = param_value,
        y = cells,
        colour = factor(pe),
        group = factor(pe),
        shape = factor(pe)
      ), size = 0.2, alpha = 0.5) +
      labs(
        x = param_labels_words_no_breaks[param_str],
        y = variable_labels[variable_str],
        colour = NULL,
        shape = NULL
      ) +
      scale_colour_discrete(
        labels = c(
          "5" = expression(Pe == 5),
          "-5" = expression(Pe == -5)
        )
      ) +
      scale_shape_discrete(
        labels = c(
          "5" = expression(Pe == 5),
          "-5" = expression(Pe == -5)
        )
      ) +
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.text.align = 0,
        # TODO re-enable the legend but find a better place to fit it in - it's
        # on the RHS at the moment and squeezes the plots a lot
        legend.position = "none"
      )

      if (variable_str == "net_change") {
        p <- p + geom_hline(yintercept = 0, linetype = "dashed")
      }

      p
  }

  # Plots
  # -----

  # Load trace data for neg Pe
  data_neg <- load_trace_data(pos = FALSE)
  integrated_fluxes_neg <- calculate_integrated_fluxes(data_neg$trace_data)
  integrated_fluxes_neg[, net_change := cells_in - cells_out]

  # Load trace data for pos Pe
  data_pos <- load_trace_data(pos = TRUE)
  integrated_fluxes_pos <- calculate_integrated_fluxes(data_pos$trace_data)
  integrated_fluxes_pos[, net_change := cells_in - cells_out]

  # Cells in
  ind_neg <- sobol_indices(
    Y = integrated_fluxes_neg$cells_in,
    N = data_neg$n_param_sample,
    params = data_neg$param_names,
    boot = TRUE,
    R = 100,
    parallel = "multicore",
    ncpus = 8
  )

  ind_pos <- sobol_indices(
    Y = integrated_fluxes_pos$cells_in,
    N = data_pos$n_param_sample,
    params = data_pos$param_names,
    boot = TRUE,
    R = 100,
    parallel = "multicore",
    ncpus = 8
  )

  p_sobol_indices_cells_in <- sobol_indices_cells_panel(
    ind_neg, ind_pos, "Cells in", first_plot = TRUE
  )

  p_sobol_indices_cells_in_2 <- sobol_indices_cells_panel_2(
    ind_neg, ind_pos
  )

  # Cells out
  ind_neg <- sobol_indices(
    Y = integrated_fluxes_neg$cells_out,
    N = data_neg$n_param_sample,
    params = data_neg$param_names,
    boot = TRUE,
    R = 100,
    parallel = "multicore",
    ncpus = 8
  )

  ind_pos <- sobol_indices(
    Y = integrated_fluxes_pos$cells_out,
    N = data_pos$n_param_sample,
    params = data_pos$param_names,
    boot = TRUE,
    R = 100,
    parallel = "multicore",
    ncpus = 8
  )

  p_sobol_indices_cells_out <- sobol_indices_cells_panel(
    ind_neg, ind_pos, "Cells out", first_plot = TRUE
  )

  p_sobol_indices_cells_out_2 <- sobol_indices_cells_panel_2(
    ind_neg, ind_pos
  )

  # Net change
  ind_neg <- sobol_indices(
    Y = integrated_fluxes_neg$net_change,
    N = data_neg$n_param_sample,
    params = data_neg$param_names,
    boot = TRUE,
    R = 100,
    parallel = "multicore",
    ncpus = 8
  )

  ind_pos <- sobol_indices(
    Y = integrated_fluxes_pos$net_change,
    N = data_pos$n_param_sample,
    params = data_pos$param_names,
    boot = TRUE,
    R = 100,
    parallel = "multicore",
    ncpus = 8
  )

  p_sobol_indices_net_change <- sobol_indices_cells_panel(
    ind_neg, ind_pos, "Net change", first_plot = TRUE
  )

  p_sobol_indices_net_change_2 <- sobol_indices_cells_panel_2(
    ind_neg, ind_pos
  )

  p_sobol_indices_cells <- p_sobol_indices_cells_in +
    p_sobol_indices_cells_out +
    p_sobol_indices_net_change +
    plot_layout(guides = "collect")

  # Cell vs params
  cells_vs_params_long_neg <- cells_vs_params_long(
    data_neg$param_sample_dt,
    integrated_fluxes_neg
  )
  cells_vs_params_long_pos <- cells_vs_params_long(
    data_pos$param_sample_dt,
    integrated_fluxes_pos
  )

  cells_vs_params_long_neg[, pe := -5]
  cells_vs_params_long_pos[, pe := 5]

  cells_vs_params_long_all <- rbind(
    cells_vs_params_long_neg,
    cells_vs_params_long_pos
  )

  p_cells_in_panels <- list()
  p_cells_out_panels <- list()
  p_net_change_panels <- list()

  param_strs <- cells_vs_params_long_all[, levels(param)]

  for (i in 1:length(param_strs)) {
    p_cells_in_panels[[i]] <- cells_vs_params_panel(
      cells_vs_params_long_all,
      "cells_in",
      param_strs[i]
    )

    p_cells_out_panels[[i]] <- cells_vs_params_panel(
      cells_vs_params_long_all,
      "cells_out",
      param_strs[i]
    )

    p_net_change_panels[[i]] <- cells_vs_params_panel(
      cells_vs_params_long_all,
      "net_change",
      param_strs[i]
    )
  }

  list(
    p_sobol_indices_cells,
    p_cells_in_panels,
    p_cells_out_panels,
    p_net_change_panels,
    p_sobol_indices_cells_in_2,
    p_sobol_indices_cells_out_2,
    p_sobol_indices_net_change_2
  )
}

plot_dir <- "plots_neg_and_pos"

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

plots <- sobol_neg_and_pos()

p_sobol_neg_and_pos <- plots[[1]]
p_sobol_neg_and_pos

ggsave_with_defaults(
  plot = p_sobol_neg_and_pos,
  paste(plot_dir, "sobol_indices_cells.pdf", sep = "/")
)

p_cells_in_panels <- plots[[2]]
p_cells_out_panels <- plots[[3]]
p_net_change_panels <- plots[[4]]
p_sobol_indices_cells_in_2 <- plots[[5]]
p_sobol_indices_cells_out_2 <- plots[[6]]
p_sobol_indices_net_change_2 <- plots[[7]]

# Hacks to get things to line up better
nudge_tag <- theme(plot.tag.position = c(-0.04, 1))
nudge_x_axis_label <- theme(
  axis.title.x = element_text(margin = margin(-88, 0, 0, 0))
)

p_sobol_and_cells_in <-
  p_sobol_indices_cells_in_2 +
    theme(legend.position = c(0.5, 0.80)) +
    (
      (p_cells_in_panels[[1]] | p_cells_in_panels[[2]]) /
        (p_cells_in_panels[[3]] + nudge_x_axis_label |
          p_cells_in_panels[[4]] + nudge_x_axis_label) +
        plot_layout(guides = "collect") &
        nudge_tag &
        labs(y = "Cells in")
    ) +
    plot_layout(widths = c(1.0 / 3.0, 2.0 / 3.0)) +
    plot_annotation(tag_levels = "A") &
    theme(
      plot.tag = element_text(hjust = 0, vjust = 0)
    )

ggsave_with_defaults(
  plot = p_sobol_and_cells_in,
  paste(plot_dir, "sobol_and_cells_in.pdf", sep = "/")
)

p_sobol_and_cells_out <-
  p_sobol_indices_cells_out_2 +
    theme(legend.position = c(0.5, 0.80)) +
    (
      (p_cells_out_panels[[1]] | p_cells_out_panels[[2]]) /
        (p_cells_out_panels[[3]] + nudge_x_axis_label |
          p_cells_out_panels[[4]] + nudge_x_axis_label) +
        plot_layout(guides = "collect") &
        nudge_tag &
        labs(y = "Cells out")
    ) +
    plot_layout(widths = c(1.0 / 3.0, 2.0 / 3.0)) +
    plot_annotation(tag_levels = "A") &
    theme(
      plot.tag = element_text(hjust = 0, vjust = 0)
    )

ggsave_with_defaults(
  plot = p_sobol_and_cells_out,
  paste(plot_dir, "sobol_and_cells_out.pdf", sep = "/")
)

p_sobol_and_cells_net_change <-
  p_sobol_indices_net_change_2 +
    theme(legend.position = c(0.5, 0.80)) +
    (
      (p_net_change_panels[[1]] | p_net_change_panels[[2]]) /
        (p_net_change_panels[[3]] + nudge_x_axis_label |
          p_net_change_panels[[4]] + nudge_x_axis_label) +
        plot_layout(guides = "collect") &
        nudge_tag &
        labs(y = "Net change")
    ) +
    plot_layout(widths = c(1.0 / 3.0, 2.0 / 3.0)) +
    plot_annotation(tag_levels = "A") &
    theme(
      plot.tag = element_text(hjust = 0, vjust = 0)
    )

ggsave_with_defaults(
  plot = p_sobol_and_cells_net_change,
  paste(plot_dir, "sobol_and_cells_net_change.pdf", sep = "/")
)

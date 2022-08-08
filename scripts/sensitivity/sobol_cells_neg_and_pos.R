library(sensitivity)
library(ggplot2)
library(cowplot)
library(patchwork)
library(GGally)

source("scripts/sensitivity/functions.R")

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
      res_dir_base <- "res_sensitivity_3_pos_pe"
    } else {
      pe <<- -5.0
      res_dir_base <- "res_sensitivity_3_neg_pe"
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

    x_1 <- gen_param_sample(100, names, mins, maxs)
    x_2 <- gen_param_sample(100, names, mins, maxs)

    x <- sobolmartinez(model = NULL, X1 = x_1, X2 = x_2, nboot = 100)

    add_constants_and_gammas_to_param_table(x$X)

    trace_data_full <- read_trace_data(x$X, res_dir_base)
    trace_data <- trace_data_full[`t_{inf}` >= 0]

    list(x, trace_data)
  }

  sobol_indices_cells_panel <- function(
    x_neg,
    x_pos,
    title,
    first_plot = TRUE
  ) {
    # rownames will be the same for both, so just use x_neg
    rn <- rownames(x_neg$T)

    data <- rbind(
      cbind(
        parameter = factor(rn, levels = rn),
        effect = "main",
        pe = "neg",
        x_neg$S %>% as.data.table()
        ),
      cbind(
        parameter = factor(rn, levels = rn),
        effect = "total",
        pe = "neg",
        x_neg$T %>% as.data.table()
      ),
      cbind(
        parameter = factor(rn, levels = rn),
        effect = "main",
        pe = "pos",
        x_pos$S %>% as.data.table()
        ),
      cbind(
        parameter = factor(rn, levels = rn),
        effect = "total",
        pe = "pos",
        x_pos$T %>% as.data.table()
      )
    )

    y_label <- if (first_plot) {
      "Sobol index"
    } else {
      NULL
    }

    p <- ggplot(
      data,
      aes(
        parameter,
        y = original,
        ymin = `min. c.i.`,
        ymax = `max. c.i.`,
        shape = effect,
        colour = pe)
      ) +
      geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
      geom_errorbar(width = 0.5, position = position_dodge(width = 0.6)) +
      scale_x_discrete(labels = param_labels) +
      scale_shape_discrete(
        labels = c("main" = "First-order", "total" = "Total-order")
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
      theme_cowplot() +
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
    x_neg,
    x_pos
  ) {
    # rownames will be the same for both, so just use x_neg
    rn <- rownames(x_neg$T)

    data <- rbind(
      cbind(
        parameter = factor(rn, levels = rn),
        effect = "main",
        pe = "neg",
        x_neg$S %>% as.data.table()
      ),
      cbind(
        parameter = factor(rn, levels = rn),
        effect = "total",
        pe = "neg",
        x_neg$T %>% as.data.table()
      ),
      cbind(
        parameter = factor(rn, levels = rn),
        effect = "main",
        pe = "pos",
        x_pos$S %>% as.data.table()
      ),
      cbind(
        parameter = factor(rn, levels = rn),
        effect = "total",
        pe = "pos",
        x_pos$T %>% as.data.table()
      )
    )

    y_label <- "Sobol index"

    p <- ggplot(
      data,
      aes(
        parameter,
        y = original,
        ymin = `min. c.i.`,
        ymax = `max. c.i.`,
        shape = effect,
        colour = pe)
      ) +
      geom_point(size = 2.5, position = position_dodge(width = 0.6)) +
      geom_errorbar(width = 0.5, position = position_dodge(width = 0.6)) +
      scale_x_discrete(labels = param_labels) +
      scale_shape_discrete(
        labels = c("main" = "First-order", "total" = "Total-order")
      ) +
      scale_colour_discrete(
        labels = c("pos" = expression(Pe == 5), "neg" = expression(Pe == -5))
      ) +
      coord_cartesian(ylim = c(0.0, 1.0)) +
      labs(
        x = NULL,
        y = y_label,
        shape = "Sobol index",
        colour = "Flow"
      ) +
      theme_cowplot() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.text.align = 0,
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
      )

    p
  }

  cells_vs_params_long <- function(x, integrated_fluxes) {
    cells_vs_params <- x$X[
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
    ggplot(cells_vs_params[variable == variable_str & param == param_str]) +
      geom_point(aes(
        x = param_value,
        y = cells,
        colour = factor(pe),
        group = factor(pe),
        shape = factor(pe)
      ), size = 1.25) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(
        x = param_labels[param_str],
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
      theme_cowplot() +
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.text.align = 0,
        # TODO re-enable the legend but find a better place to fit it in - it's
        # on the RHS at the moment and squeezes the plots a lot
        legend.position = "none"
      )
  }

  # Plots
  # -----

  # Load trace data for neg Pe
  data_neg <- load_trace_data(pos = FALSE)
  x_neg <- data_neg[[1]]
  trace_data_neg <- data_neg[[2]]
  integrated_fluxes_neg <- calculate_integrated_fluxes(trace_data_neg)
  integrated_fluxes_neg[, net_change := cells_in - cells_out]

  # Load trace data for pos Pe
  data_pos <- load_trace_data(pos = TRUE)
  x_pos <- data_pos[[1]]
  trace_data_pos <- data_pos[[2]]
  integrated_fluxes_pos <- calculate_integrated_fluxes(trace_data_pos)
  integrated_fluxes_pos[, net_change := cells_in - cells_out]

  # Cells in
  tell(x_neg, integrated_fluxes_neg$cells_in)
  tell(x_pos, integrated_fluxes_pos$cells_in)

  p_sobol_indices_cells_in <- sobol_indices_cells_panel(
    x_neg, x_pos, "Cells in", first_plot = TRUE
  )

  p_sobol_indices_cells_in_2 <- sobol_indices_cells_panel_2(
    x_neg, x_pos
  )

  # Cells out
  tell(x_neg, integrated_fluxes_neg$cells_out)
  tell(x_pos, integrated_fluxes_pos$cells_out)

  p_sobol_indices_cells_out <- sobol_indices_cells_panel(
    x_neg, x_pos, "Cells out", first_plot = FALSE
  )

  p_sobol_indices_cells_out_2 <- sobol_indices_cells_panel_2(
    x_neg, x_pos
  )

  # Net change
  tell(x_neg, integrated_fluxes_neg$net_change)
  tell(x_pos, integrated_fluxes_pos$net_change)

  p_sobol_indices_net_change <- sobol_indices_cells_panel(
    x_neg, x_pos, "Net change", first_plot = FALSE
  )

  p_sobol_indices_net_change_2 <- sobol_indices_cells_panel_2(
    x_neg, x_pos
  )

  p_sobol_indices_cells <- p_sobol_indices_cells_in +
    p_sobol_indices_cells_out +
    p_sobol_indices_net_change +
    plot_layout(guides = "collect")

  # Cell vs params
  cells_vs_params_long_neg <- cells_vs_params_long(x_neg, integrated_fluxes_neg)
  cells_vs_params_long_pos <- cells_vs_params_long(x_pos, integrated_fluxes_pos)

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
  paste(plot_dir, "sobol_indices_cells.pdf", sep = "/"),
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
  axis.title.x = element_text(margin = margin(-40, 0, 0, 0))
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

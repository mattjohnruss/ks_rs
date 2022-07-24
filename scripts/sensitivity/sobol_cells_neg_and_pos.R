library(sensitivity)
library(ggplot2)
library(cowplot)
library(patchwork)
library(GGally)

source("scripts/sensitivity/functions.R")

sobol_neg_and_pos <- function() {
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

  sobol_indices_cells_panel <- function(x_neg, x_pos, title, first_plot = TRUE) {
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

  data_neg <- load_trace_data(pos = FALSE)
  x_neg <- data_neg[[1]]
  trace_data_neg <- data_neg[[2]]
  integrated_fluxes_neg <- calculate_integrated_fluxes(trace_data_neg)
  integrated_fluxes_neg[, net_change := cells_in - cells_out]

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

  # Cells out
  tell(x_neg, integrated_fluxes_neg$cells_out)
  tell(x_pos, integrated_fluxes_pos$cells_out)

  p_sobol_indices_cells_out <- sobol_indices_cells_panel(
    x_neg, x_pos, "Cells out", first_plot = FALSE
  )

  # Net change
  tell(x_neg, integrated_fluxes_neg$net_change)
  tell(x_pos, integrated_fluxes_pos$net_change)

  p_sobol_indices_cells_net_change <- sobol_indices_cells_panel(
    x_neg, x_pos, "Net change", first_plot = FALSE
  )

  p_sobol_indices_cells_in +
    p_sobol_indices_cells_out +
    p_sobol_indices_cells_net_change +
    plot_layout(guides = "collect")
}

plot_dir <- "plots_neg_and_pos"

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

p_sobol_neg_and_pos <- sobol_neg_and_pos()

p_sobol_neg_and_pos

ggsave_with_defaults(
  plot = p_sobol_neg_and_pos,
  paste(plot_dir, "sobol_indices_cells.pdf", sep = "/"),
)

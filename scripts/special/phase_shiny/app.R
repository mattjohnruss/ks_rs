library(ggplot2)
library(cowplot)

# Read the trace data for the "special" runs

###############################################################################
# COPY of some of the contents of functions.R here to make deploying shiny app
# easier. Not a good way to do this...
###############################################################################

library(data.table)
library(magrittr)

# TODO: having these as globals is awful
# Constant parameter values
phi_bar_over_c_bar <- 1.40179e-6
phi_bar_over_phi_max <- 0.1
c_bar_over_e <- 0.1
pe <- -5.0
alpha_plus <- 23.25
alpha_minus <- 0.3
beta_plus <- 0.0000385493
beta_minus <- 12.5
n_ccr7 <- 30000.0
q_u <- 0.0
q_b <- 0.0
q_s <- 0.0
d_c_s <- 0.01
d_phi_i <- 0.01
d_phi_m <- 0.01
d_phi_c_u <- 0.01
d_phi_c_b <- 0.01
chi_u <- 0.0
chi_b <- 0.004
chi_s <- 0.0
r <- 0.0
m_h <- 0.00113865
t_inflammation <- 50.0
j_phi_i_h <- 0.000455461
phi_i_init <- 0.9
phi_m_init <- 0.0

# TODO: this function should really only set the parameters that haven't
# already been set. Currently the behaviour of each line depends on whether
# the corresponding column is already set: if it's not set, it is created with
# the value of the global variable; if it is set, it is set to itself (i.e. no
# change)! This kind-of does what we want but isn't good.
add_constants_and_gammas_to_param_table <- function(params) {
  params[, phi_bar_over_c_bar := phi_bar_over_c_bar]
  params[, phi_bar_over_phi_max := phi_bar_over_phi_max]
  params[, c_bar_over_e := c_bar_over_e]
  params[, pe := pe]
  params[, alpha_plus := alpha_plus]
  params[, alpha_minus := alpha_minus]
  params[, beta_plus := beta_plus]
  params[, beta_minus := beta_minus]
  params[, n_ccr7 := n_ccr7]
  params[, q_u := q_u]
  params[, q_b := q_b]
  params[, q_s := q_s]
  params[, d_c_s := d_c_s]
  params[, d_phi_i := d_phi_i]
  params[, d_phi_m := d_phi_m]
  params[, d_phi_c_u := d_phi_c_u]
  params[, d_phi_c_b := d_phi_c_b]
  params[, chi_u := chi_u]
  params[, chi_b := chi_b]
  params[, chi_s := chi_s]
  params[, r := r]
  params[, m_h := m_h]
  params[, t_inflammation := t_inflammation]
  params[, j_phi_i_h := j_phi_i_h]
  params[, phi_i_init := phi_i_init]
  params[, phi_m_init := phi_m_init]
  params[, gamma_ui := gamma]
  params[, gamma_um := gamma]
  params[, gamma_bi := gamma]
  params[, gamma_bm := gamma]
}

read_trace_data <- function(all_params, res_dir_base) {
  n_runs <- nrow(all_params)
  trace_data <- list()

  for (rep in 1:n_runs) {
    filename <- paste(res_dir_base, rep - 1, "trace.csv", sep = "/")
    print(filename)
    trace_data[[rep]] <- data.table::fread(filename)
    trace_data[[rep]][, rep := rep - 1]
    trace_data[[rep]][`t_{inf}` >= 0, output_inf := .I]
    trace_data[[rep]] <- cbind(trace_data[[rep]], all_params[rep])
  }

  return(rbindlist(trace_data))
}

# END of functions.R copy

set.seed(12345L)

theme_set(theme_cowplot() + background_grid())

res_dir_base <- "res_special_trace_only"

params <- fread(paste(res_dir_base, "d_m.csv", sep = "/"), sep = " ")

add_constants_and_gammas_to_param_table(params)

trace_data <- read_trace_data(params, res_dir_base)[`t_{inf}` >= 0]

# UI

ui <- fluidPage(
  flowLayout(
    varSelectInput(
      "x_var",
      "x-axis",
      trace_data[, .(
        t,
        `C_u^{tot}`,
        `C_b^{tot}`,
        `C_s^{tot}`,
        `phi_i^{tot}`,
        `phi_m^{tot}`,
        `phi_{C_u}^{tot}`,
        `phi_{C_b}^{tot}`,
        `-F_{phi_i}(x=1)`,
        `-F_{phi_{C_b}}(x=0)`
        )],
      selected = "C_b^{tot}"
    ),
    varSelectInput(
      "y_var",
      "y-axis",
      trace_data[, .(
        t,
        `C_u^{tot}`,
        `C_b^{tot}`,
        `C_s^{tot}`,
        `phi_i^{tot}`,
        `phi_m^{tot}`,
        `phi_{C_u}^{tot}`,
        `phi_{C_b}^{tot}`,
        `-F_{phi_i}(x=1)`,
        `-F_{phi_{C_b}}(x=0)`
        )],
      selected = "phi_{C_b}^{tot}"
    ),
    selectInput(
      "colour_by",
      "Colour by",
      c("t_{inf}", "rep", "j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag", "gamma", "pe"),
      selected = "pe"
    ),
    selectInput(
      "linetype_by",
      "Linetype by",
      c("none", "t_{inf}", "rep", "j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag", "gamma", "pe"),
    ),
    column(
      12,
      selectInput(
        "alpha_by",
        "Alpha by",
        c("none", "t_{inf}", "rep", "j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag", "gamma", "pe")
      ),
      checkboxInput(
        "reverse_alpha",
        "Reverse alpha direction",
        value = FALSE
      )
    ),
    checkboxGroupInput(
      "j_phi_i_i_factor",
      "j_phi_i_i_factor",
      c("2" = "2", "1000" = "1000"),
      selected = c("2")
    ),
    checkboxGroupInput(
      "m_i_factor",
      "m_i_factor",
      c("2" = "2", "1000" = "1000"),
      selected = c("2")
    ),
    checkboxGroupInput(
      "t_j_phi_i_lag",
      "t_j_phi_i_lag",
      c("0" = "0", "25" = "25"),
      selected = c("0")
    ),
    checkboxGroupInput(
      "gamma",
      "gamma",
      c("0" = "0", "1" = "1"),
      selected = c("0", "1")
    ),
    checkboxGroupInput(
      "pe",
      "pe",
      c(
        "-5" = "-5",
        "-4" = "-4",
        "-3" = "-3",
        "-2" = "-2",
        "-1" = "-1",
        "0" = "0",
        "1" = "1",
        "2" = "2",
        "3" = "3",
        "4" = "4",
        "5" = "5"
      ),
      selected = c("-5", "-3", "-1", "1", "3", "5"),
      inline = TRUE
    )
  ),
  plotOutput(outputId = "distPlot", height = "1200"),
  theme = bslib::bs_theme(version = 5)
)

# Server

server <- function(input, output) {
  output$distPlot <- renderPlot({
    trace_data_subset <- trace_data[
      j_phi_i_i_factor %in% as.numeric(input$j_phi_i_i_factor) &
        m_i_factor %in% as.numeric(input$m_i_factor) &
        t_j_phi_i_lag %in% as.numeric(input$t_j_phi_i_lag) &
        gamma %in% as.numeric(input$gamma) &
        pe %in% as.numeric(input$pe)
      ]

    p <- ggplot(trace_data_subset,
      aes(
        x = !!input$x_var,
        y = !!input$y_var,
        group = rep
      ),
    )

    # colour
    if (input$colour_by == "t_{inf}") {
      p <- p + aes(colour = `t_{inf}`) +
        scale_colour_distiller(palette = "RdBu")
    } else {
      p <- p + aes(colour = factor(!!sym(input$colour_by))) +
        guides(colour = guide_legend(ncol = 1))
    }

    # linetype
    if (input$linetype_by != "none") {
      p <- p + aes(linetype = factor(!!sym(input$linetype_by)))
    }

    # alpha
    if (input$alpha_by != "none") {
      p <- p + aes(alpha = !!sym(input$alpha_by))

      if (input$reverse_alpha) {
        p <- p + scale_alpha_continuous(trans = "reverse")
      }
    }

    p <- p +
      geom_path(
        data = function(d) d[state == 1],
        size = 1.5
      ) +
      geom_path(
        data = function(d) d[state == 0],
        size = 1.5
      )

    p
  })
}

shinyApp(ui = ui, server = server)

library(ggplot2)
library(cowplot)
library(plotly)
library(shiny)

source("scripts/sensitivity/functions.R")

theme_set(theme_cowplot())

res_dir_base <- "res_special"

params <- expand.grid(
  c(2, 1000),
  c(2, 1000),
  c(0, 25),
  c(0, 1),
  -5:5 %>% as.numeric
) %>%
  as.data.table %>%
  setnames(c("j_phi_i_i_factor", "m_i_factor", "t_j_phi_i_lag", "gamma", "pe"))

# we can luckily use this function without it overwriting our specified pe
# values, as within the function it does params[, pe := pe], which just sets
# them to themselves, if they have already been assigned (if they haven't, the
# pe on the RHS references the global variable)
add_constants_and_gammas_to_param_table(params)

#write_json_params_files(params)
#fwrite(params, paste(res_dir_base, "d_m.csv", sep = "/"), sep = " ")

#params_test <- fread(paste(res_dir_base, "d_m.csv", sep = "/"), sep = " ")

#simulate_all(params)

trace_data_full <- read_trace_data(params, res_dir_base)
trace_data <- trace_data_full[`t_{inf}` >= 0]

trace_data_long <- melt(
  trace_data,
  measure.vars = c(
    "C_u^{tot}",
    "C_b^{tot}",
    "C_s^{tot}",
    "phi_i^{tot}",
    "phi_m^{tot}",
    "phi_{C_u}^{tot}",
    "phi_{C_b}^{tot}",
    "-F_{phi_i}(x=1)",
    "-F_{phi_{C_b}}(x=0)"
  )
)

ui <- fluidPage(
  flowLayout(
    varSelectInput(
      "colour_by",
      "Colour by",
      params[, .(j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma, pe)],
      selected = "pe"
    ),
    varSelectInput(
      "linetype_by",
      "Linetype by",
      params[, .(j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma, pe)]
    ),
    varSelectInput(
      "alpha_by",
      "Alpha by",
      params[, .(j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma, pe)]
    ),
    checkboxGroupInput(
      "j_phi_i_i_factor",
      "j_phi_i_i_factor",
      c("2" = "2", "1000" = "1000"),
      selected = c("2", "1000")
    ),
    checkboxGroupInput(
      "m_i_factor",
      "m_i_factor",
      c("2" = "2", "1000" = "1000"),
      selected = c("2", "1000")
    ),
    checkboxGroupInput(
      "t_j_phi_i_lag",
      "t_j_phi_i_lag",
      c("0" = "0", "25" = "25"),
      selected = c("0", "25")
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
      selected = c("-5", "-2", "2", "5"),
      inline = TRUE
    )
  ),
  plotOutput(outputId = "distPlot", height = "1200"),
  theme = bslib::bs_theme(version = 5)
)

server <- function(input, output) {
  output$distPlot <- renderPlot({
    ggplot(trace_data_long[
      j_phi_i_i_factor %in% as.numeric(input$j_phi_i_i_factor) &
        m_i_factor %in% as.numeric(input$m_i_factor) &
        t_j_phi_i_lag %in% as.numeric(input$t_j_phi_i_lag) &
        gamma %in% as.numeric(input$gamma) &
        pe %in% as.numeric(input$pe)
      ],
      aes(
        x = `t_{inf}`,
        y = value,
        group = rep,
        colour = factor(!!input$colour_by),
        linetype = factor(!!input$linetype_by),
        alpha = factor(!!input$alpha_by)
      )
    ) +
    geom_line(size = 1.5) +
    scale_alpha_manual(values = c(0.75, 0.3)) +
    facet_wrap(vars(variable), scales = "free")
  })
}

shinyApp(ui = ui, server = server)

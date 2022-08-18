library(ggplot2)
library(cowplot)

source("../../sensitivity/functions.R")

theme_set(theme_cowplot() + background_grid())

res_dir_base <- "../../../res_special"

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

ui <- fluidPage(
  flowLayout(
    sliderInput(
      inputId = "time",
      label = "Time since inflammation:",
      min = 0,
      max = 1500,
      value = 500),
    varSelectInput(
      "colour_by",
      "Colour by",
      params[, .(j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma, pe)],
      selected = "pe"
    ),
    varSelectInput(
      "linetype_by",
      "Linetype by",
      params[, .(j_phi_i_i_factor, m_i_factor, t_j_phi_i_lag, gamma, pe)],
      selected = "gamma"
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
    ),
    checkboxGroupInput(
      "variables",
      "Variables",
      c("C_u", "C_b", "C_s", "phi_i", "phi_m", "phi_C_u", "phi_C_b"),
      selected = c("C_b", "phi_C_b"),
      inline = TRUE
    )
  ),
  plotOutput(outputId = "homeostasisPlot", height = "600"),
  plotOutput(outputId = "laterPlot", height = "600"),
  theme = bslib::bs_theme(version = 5)
)

server <- function(input, output) {
  output$homeostasisPlot <- renderPlot({
    data <- read_spatial_data(params, 0, res_dir_base)
    data_long <- melt(
      data,
      measure.vars = c("C_u", "C_b", "C_s", "phi_i", "phi_m", "phi_C_u", "phi_C_b")
    )
    ggplot(data_long[
      variable %in% input$variables &
        j_phi_i_i_factor == 2 &
        m_i_factor == 2 &
        t_j_phi_i_lag == 0 &
        gamma %in% as.numeric(input$gamma) &
        pe %in% as.numeric(input$pe)
      ],
      aes(
        x = x,
        y = value,
        group = rep,
        colour = factor(!!input$colour_by),
        linetype = factor(!!input$linetype_by)
      )
    ) +
    geom_line(size = 1.5) +
    facet_wrap(
      vars(variable), scales = "free",
      labeller = as_labeller(function(v) all_labels[v], label_parsed)
    ) +
    labs(y = NULL)
  })
  output$laterPlot <- renderPlot({
    data <- read_spatial_data(params, input$time, res_dir_base)
    data_long <- melt(
      data,
      measure.vars = c("C_u", "C_b", "C_s", "phi_i", "phi_m", "phi_C_u", "phi_C_b")
    )
    ggplot(data_long[
      variable %in% input$variables &
        j_phi_i_i_factor %in% as.numeric(input$j_phi_i_i_factor) &
        m_i_factor %in% as.numeric(input$m_i_factor) &
        t_j_phi_i_lag %in% as.numeric(input$t_j_phi_i_lag) &
        gamma %in% as.numeric(input$gamma) &
        pe %in% as.numeric(input$pe)
      ],
      aes(
        x = x,
        y = value,
        group = rep,
        colour = factor(!!input$colour_by),
        linetype = factor(!!input$linetype_by)
      )
    ) +
    geom_line(size = 1.5) +
    facet_wrap(
      vars(variable), scales = "free",
      labeller = as_labeller(function(v) all_labels[v], label_parsed)
    ) +
    labs(y = NULL)
  })
}

shinyApp(ui = ui, server = server)

plot_fastStructure <-
function(df) {
  df |>
    ggplot() +
    geom_line(
      aes(
        x              = run,
        y              = k,
        color          = model
      ),
      linewidth = 1
    ) +
    scale_colour_manual(
      "model",
      values = c(
        structure      = "#9d60ff",
        likelihood     = "#ffc08c"
      ),
      labels = c(
        "Maximizes \n Likelihood \n", "Explain \n Structure"
      )
    ) +
    labs(
      x                = "Run",
      y                = "K",
      title            = "fastStructure",
      caption          = "algorithm runs for choices of K ranging from 1 to 30"
    ) +
    hrbrthemes::theme_ipsum(
      base_family = "",
      axis_text_size = 12,
      axis_title_size = 14,
      plot_margin = margin(
        10, 10, 10, 10
      ),
      grid = TRUE,
      grid_col = "#fabbe2"
    ) +
    theme(
      panel.grid.major = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      panel.grid.minor = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      legend.text = element_text(
        size = 12
      ),
      legend.title = element_text(
        size           = 14,
        face           = "bold"
      ),
      legend.position = "right"
    )
}
plot_fastStructure <-
function(df) {
  df |>
    ggplot() +
    geom_line(
      aes(
        x              = run,
        y              = k,
        color          = model
      ),
      linewidth = 1
    ) +
    scale_colour_manual(
      "model",
      values = c(
        structure      = "#9d60ff",
        likelihood     = "#ffc08c"
      ),
      labels = c(
        "Maximizes \n Likelihood \n", "Explain \n Structure"
      )
    ) +
    labs(
      x                = "Run",
      y                = "K",
      title            = "fastStructure",
      caption          = "algorithm runs for choices of K ranging from 1 to 30"
    ) +
    hrbrthemes::theme_ipsum(
      base_family = "",
      axis_text_size = 12,
      axis_title_size = 14,
      plot_margin = margin(
        10, 10, 10, 10
      ),
      grid = TRUE,
      grid_col = "#fabbe2"
    ) +
    theme(
      panel.grid.major = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      panel.grid.minor = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      legend.text = element_text(
        size = 12
      ),
      legend.title = element_text(
        size           = 14,
        face           = "bold"
      ),
      legend.position = "right"
    )
}
plot_fastStructure <-
function(df) {
  df |>
    ggplot() +
    geom_line(
      aes(
        x              = run,
        y              = k,
        color          = model
      ),
      linewidth = 1
    ) +
    scale_colour_manual(
      "model",
      values = c(
        structure      = "#9d60ff",
        likelihood     = "#ffc08c"
      ),
      labels = c(
        "Maximizes \n Likelihood \n", "Explain \n Structure"
      )
    ) +
    labs(
      x                = "Run",
      y                = "K",
      title            = "fastStructure",
      caption          = "algorithm runs for choices of K ranging from 1 to 30"
    ) +
    hrbrthemes::theme_ipsum(
      base_family = "",
      axis_text_size = 12,
      axis_title_size = 14,
      plot_margin = margin(
        10, 10, 10, 10
      ),
      grid = TRUE,
      grid_col = "#fabbe2"
    ) +
    theme(
      panel.grid.major = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      panel.grid.minor = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      legend.text = element_text(
        size = 12
      ),
      legend.title = element_text(
        size           = 14,
        face           = "bold"
      ),
      legend.position = "right"
    )
}
plot_fastStructure <-
function(df) {
  df |>
    ggplot() +
    geom_line(
      aes(
        x              = run,
        y              = k,
        color          = model
      ),
      linewidth = 1
    ) +
    scale_colour_manual(
      "model",
      values = c(
        structure      = "#9d60ff",
        likelihood     = "#ffc08c"
      ),
      labels = c(
        "Maximizes \n Likelihood \n", "Explain \n Structure"
      )
    ) +
    labs(
      x                = "Run",
      y                = "K",
      title            = "fastStructure",
      caption          = "algorithm runs for choices of K ranging from 1 to 30"
    ) +
    hrbrthemes::theme_ipsum(
      base_family = "",
      axis_text_size = 12,
      axis_title_size = 14,
      plot_margin = margin(
        10, 10, 10, 10
      ),
      grid = TRUE,
      grid_col = "#fabbe2"
    ) +
    theme(
      panel.grid.major = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      panel.grid.minor = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      legend.text = element_text(
        size = 12
      ),
      legend.title = element_text(
        size           = 14,
        face           = "bold"
      ),
      legend.position = "right"
    )
}
plot_fastStructure <-
function(df) {
  df |>
    ggplot() +
    geom_line(
      aes(
        x              = run,
        y              = k,
        color          = model
      ),
      linewidth = 1
    ) +
    scale_colour_manual(
      "model",
      values = c(
        structure      = "#9d60ff",
        likelihood     = "#ffc08c"
      ),
      labels = c(
        "Maximizes \n Likelihood \n", "Explain \n Structure"
      )
    ) +
    labs(
      x                = "Run",
      y                = "K",
      title            = "fastStructure",
      caption          = "algorithm runs for choices of K ranging from 1 to 30"
    ) +
    hrbrthemes::theme_ipsum(
      base_family = "",
      axis_text_size = 12,
      axis_title_size = 14,
      plot_margin = margin(
        10, 10, 10, 10
      ),
      grid = TRUE,
      grid_col = "#fabbe2"
    ) +
    theme(
      panel.grid.major = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      panel.grid.minor = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      legend.text = element_text(
        size = 12
      ),
      legend.title = element_text(
        size           = 14,
        face           = "bold"
      ),
      legend.position = "right"
    )
}
plot_fastStructure <-
function(df) {
  df |>
    ggplot() +
    geom_line(
      aes(
        x              = run,
        y              = k,
        color          = model
      ),
      linewidth = 1
    ) +
    scale_colour_manual(
      "model",
      values = c(
        structure      = "#9d60ff",
        likelihood     = "#ffc08c"
      ),
      labels = c(
        "Maximizes \n Likelihood \n", "Explain \n Structure"
      )
    ) +
    labs(
      x                = "Run",
      y                = "K",
      title            = "fastStructure",
      caption          = "algorithm runs for choices of K ranging from 1 to 30"
    ) +
    hrbrthemes::theme_ipsum(
      base_family = "",
      axis_text_size = 12,
      axis_title_size = 14,
      plot_margin = margin(
        10, 10, 10, 10
      ),
      grid = TRUE,
      grid_col = "#fabbe2"
    ) +
    theme(
      panel.grid.major = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      panel.grid.minor = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      legend.text = element_text(
        size = 12
      ),
      legend.title = element_text(
        size           = 14,
        face           = "bold"
      ),
      legend.position = "right"
    )
}
plot_fastStructure <-
function(df) {
  df |>
    ggplot() +
    geom_line(
      aes(
        x              = run,
        y              = k,
        color          = model
      ),
      linewidth = 1
    ) +
    scale_colour_manual(
      "model",
      values = c(
        structure      = "#9d60ff",
        likelihood     = "#ffc08c"
      ),
      labels = c(
        "Maximizes \n Likelihood \n", "Explain \n Structure"
      )
    ) +
    labs(
      x                = "Run",
      y                = "K",
      title            = "fastStructure",
      caption          = "algorithm runs for choices of K ranging from 1 to 30"
    ) +
    hrbrthemes::theme_ipsum(
      base_family = "",
      axis_text_size = 12,
      axis_title_size = 14,
      plot_margin = margin(
        10, 10, 10, 10
      ),
      grid = TRUE,
      grid_col = "#fabbe2"
    ) +
    theme(
      panel.grid.major = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      panel.grid.minor = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      legend.text = element_text(
        size = 12
      ),
      legend.title = element_text(
        size           = 14,
        face           = "bold"
      ),
      legend.position = "right"
    )
}
plot_fastStructure <-
function(df) {
  df |>
    ggplot() +
    geom_line(
      aes(
        x              = run,
        y              = k,
        color          = model
      ),
      linewidth = 1
    ) +
    scale_colour_manual(
      "model",
      values = c(
        structure      = "#9d60ff",
        likelihood     = "#ffc08c"
      ),
      labels = c(
        "Maximizes \n Likelihood \n", "Explain \n Structure"
      )
    ) +
    labs(
      x                = "Run",
      y                = "K",
      title            = "fastStructure",
      caption          = "algorithm runs for choices of K ranging from 1 to 30"
    ) +
    hrbrthemes::theme_ipsum(
      base_family = "",
      axis_text_size = 12,
      axis_title_size = 14,
      plot_margin = margin(
        10, 10, 10, 10
      ),
      grid = TRUE,
      grid_col = "#fabbe2"
    ) +
    theme(
      panel.grid.major = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      panel.grid.minor = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      legend.text = element_text(
        size = 12
      ),
      legend.title = element_text(
        size           = 14,
        face           = "bold"
      ),
      legend.position = "right"
    )
}
plot_fastStructure <-
function(df) {
  df |>
    ggplot() +
    geom_line(
      aes(
        x              = run,
        y              = k,
        color          = model
      ),
      linewidth = 1
    ) +
    scale_colour_manual(
      "model",
      values = c(
        structure      = "#9d60ff",
        likelihood     = "#ffc08c"
      ),
      labels = c(
        "Maximizes \n Likelihood \n", "Explain \n Structure"
      )
    ) +
    labs(
      x                = "Run",
      y                = "K",
      title            = "fastStructure",
      caption          = "algorithm runs for choices of K ranging from 1 to 30"
    ) +
    hrbrthemes::theme_ipsum(
      base_family = "",
      axis_text_size = 12,
      axis_title_size = 14,
      plot_margin = margin(
        10, 10, 10, 10
      ),
      grid = TRUE,
      grid_col = "#fabbe2"
    ) +
    theme(
      panel.grid.major = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      panel.grid.minor = element_line(
        linetype       = "dashed",
        linewidth      = 0.2
      ),
      legend.text = element_text(
        size = 12
      ),
      legend.title = element_text(
        size           = 14,
        face           = "bold"
      ),
      legend.position = "right"
    )
}

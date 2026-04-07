my_theme <-
function() {
  theme_minimal(base_size = 12, base_family = "") +
    theme(
      panel.grid.major = element_line(
        linetype = "dashed",
        linewidth = 0.2,
        color = "pink"
      ),
      panel.grid.minor = element_line(
        linetype = "dashed",
        linewidth = 0.2,
        color = "pink"
      ),
      # Customize the x-axis label
      axis.title.x = element_text(
        angle          = 0,
        hjust          = 1,
        face           = "bold"
      ),
      # Customize the y-axis label
      axis.title.y = element_text(
        angle          = 90,
        hjust          = 1,
        face           = "bold"
      )
    )
}

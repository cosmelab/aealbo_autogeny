my_theme <- function() {
  theme_minimal(base_size = 12, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial"),  # This line sets the font for all text elements to Arial
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(
        angle = 0,
        hjust = 1,
        face = "bold",
        family = "Arial"  # This line sets the font for the x-axis title to Arial
      ),
      axis.title.y = element_text(
        angle = 90,
        hjust = 1,
        face = "bold",
        family = "Arial"  # This line sets the font for the y-axis title to Arial
      )
    )
}

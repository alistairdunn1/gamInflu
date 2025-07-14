plot_influence_single <- function(influence_data, x, ...) {
  if (nrow(influence_data) == 0) {
    warning("No data to plot")
    return(ggplot())
  }

  # Get clean level names
  clean_levels <- extract_clean_levels(names(x$focus_effects$relative_effects), x$focus)

  p <- ggplot(influence_data, aes(x = level_number, y = influence, color = term)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
    geom_line(aes(group = term)) +
    geom_point() +
    scale_x_continuous(
      breaks = seq_along(clean_levels),
      labels = clean_levels
    ) + # Now uses clean level names
    scale_color_discrete(name = "Terms") +
    labs(x = paste("Focus term:", x$focus), y = "Influence")

  return(p)
}

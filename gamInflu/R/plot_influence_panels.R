plot_influence_panels <- function(influence_data, x, ...) {
  if (nrow(influence_data) == 0) {
    warning("No data to plot")
    return(ggplot())
  }

  # Get clean level names
  clean_levels <- extract_clean_levels(names(x$focus_effects$relative_effects), x$focus)

  p <- ggplot(influence_data, aes(x = level_number, y = influence)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
    geom_line(color = "black") +
    geom_point(color = "black") +
    facet_wrap(~term, ncol = 1, scales = "free_y") +
    scale_x_continuous(breaks = seq_along(clean_levels), labels = clean_levels) +
    labs(x = paste("Focus term:", x$focus), y = "Influence")

  return(p)
}

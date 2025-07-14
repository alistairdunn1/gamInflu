#' Step plot for GAM Terms
#'
#' Plot the step change of terms on the parametric focus term
#'
#' @param step_effects A data frame containing the step effects
#' @param x A gam_influence object with calculations performed
#' @param show_previous Logical indicating whether to show previous steps
#' @param ncol (integer) Number of columns for faceting
#' @return A ggplot object
#' @keywords internal
#'
plot_step_panels <- function(step_effects, x, show_previous = TRUE, ncol = 1) {
  plot_data <- prepare_step_data(step_effects, x, show_previous)
  labs <- unique(plot_data$term_added)[-1]

  if (nrow(plot_data) == 0) {
    warning("No data to plot")
    return(ggplot())
  }

  # Get clean level names
  clean_levels <- extract_clean_levels(names(x$focus_effects$relative_effects), x$focus)

  # Create the plot
  p <- ggplot(plot_data, aes(x = level, y = effect)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray50", alpha = 0.7) +

    # Add previous steps as colour lines
    {
      if (show_previous) {
        list(
          geom_line(data = subset(plot_data, is_previous), aes(group = interaction(step_number, line_type), colour = line_type), alpha = 0.5),
          geom_point(data = subset(plot_data, is_previous), aes(colour = line_type), alpha = 0.5)
        )
      }
    } +

    # Add current step
    geom_line(data = subset(plot_data, is_current), aes(group = step_number), colour = "black") +
    geom_point(data = subset(plot_data, is_current), colour = "black") +

    # Facet by step
    facet_wrap(~step_name, ncol = ncol, scales = "free_y", strip.position = "left") +

    # Modified x-axis to use clean level names
    scale_x_discrete(labels = clean_levels) +
    scale_colour_manual(values = as.factor(unique(plot_data$line_type)[-1]), labels = paste("+", unique(plot_data$term_added)[-1], sep = "")) +
    theme(legend.title = element_blank()) +
    labs(x = x$focus, y = "Relative Effect")

  return(p)
}

#' Step plot for GAM Terms
#'
#' Plot the step change of terms on the parametric focus term
#'
#' @param step_effects A data frame containing the step effects
#' @param x A gam_influence object with calculations performed
#' @param show_previous Logical indicating whether to show previous steps
#' @return A ggplot object
#' @keywords internal
#'
plot_step_single <- function(step_effects, x, show_previous = TRUE) {
  if (length(step_effects) == 0) {
    warning("No data to plot")
    return(ggplot())
  }

  # Get clean level names
  clean_levels <- extract_clean_levels(names(x$focus_effects$relative_effects), x$focus)

  # Prepare data for single panel
  plot_data <- data.frame()

  for (i in seq_along(step_effects)) {
    step <- step_effects[[i]]
    if (!is.null(step$effects) && length(step$effects) > 0) {
      step_data <- data.frame(
        step_number = i,
        step_name = paste(i, ". ", step$term_added, sep = ""),
        step_labels = paste("+", step$term_added, sep = ""),
        term_added = step$term_added,
        level = clean_levels,
        level_number = seq_along(clean_levels),
        effect = step$effects,
        is_final = i == length(step_effects)
      )
      plot_data <- rbind(plot_data, step_data)
    }
  }

  if (nrow(plot_data) == 0) {
    warning("No valid step effects to plot")
    return(ggplot())
  }

  # Create the plot
  p <- ggplot(plot_data, aes(x = level, y = effect)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "gray50", alpha = 0.7)

  if (!show_previous) {
    p <- p +
      geom_line(data = subset(plot_data, is_final), aes(group = step_labels), colour = "black") +
      geom_point(data = subset(plot_data, is_final), aes(group = step_labels), colour = "black")
  } else {
    p <- p +
      geom_line(aes(group = step_labels, colour = factor(step_labels))) +
      geom_point(aes(group = step_labels, colour = factor(step_labels))) +
      guides(colour = guide_legend(override.aes = list())) +
      theme(legend.title = element_blank())
  }

  p <- p +
    scale_x_discrete() + # Let ggplot handle the discrete scale automatically
    labs(x = x$focus, y = "Relative Effect")

  return(p)
}

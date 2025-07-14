#' Step Plot for GAM Influence
#'
#' Create a step plot showing how the parametric focus term changes as
#' smooth terms are added to the model. Each panel shows the cumulative
#' effect with all previous steps displayed as light gray reference lines.
#'
#' @param x A gam_influence object with calculations performed
#' @param panels Whether to use separate panels (TRUE) or single panel (FALSE)
#' @param show_previous Whether to show all previous steps as colour (TRUE).
#'        This helps visualize the cumulative effect of adding each term.
#' @param ncol Number of columns to use for the panels in the plot (default = 1)
#'
#' @details
#' When \code{show_previous = TRUE} (default), each panel shows:
#' \itemize{
#'   \item All previous model steps as colour lines
#'   \item Current step highlighted in black with larger points
#'   \item Reference line at y = 1 (no effect)
#' }
#' This visualization makes it easy to see how each additional term
#' affects the focus term estimates relative to simpler models.
#'
#' @return A ggplot object
#' @export
plot_step <- function(x, panels = TRUE, show_previous = TRUE, ncol = 1) {
  if (!x$calculated) {
    stop("Must run calculate() first")
  }

  # Get focus term coefficients for each step
  step_effects <- calculate_step_effects(x)

  if (panels) {
    return(plot_step_panels(step_effects, x, show_previous = show_previous, ncol = ncol))
  } else {
    return(plot_step_single(step_effects, x, show_previous = show_previous))
  }
}

#' @title Combined Step and Influence Plot
#' @description Creates a combined plot showing both the stepwise index and term influence for a `gam_influence` object.
#' The plots are aligned horizontally so that each stepwise panel on the left corresponds to the same row
#' as the influence panels on the right.
#' @param obj A `gam_influence` object containing calculated indices.
#' @return A patchwork plot combining the stepwise index and term influence plots with proper alignment.
#' @importFrom patchwork plot_layout
#' @importFrom ggplot2 ggplot labs theme_void
#' @export
plot_step_and_influence <- function(obj) {
  p_step <- plot_stepwise_index(obj)
  p_influ <- plot_term_influence(obj)

  # Create a blank plot to align with the "Intercept" panel from stepwise plot
  blank_plot <- ggplot2::ggplot() +
    ggplot2::labs(title = "") +
    ggplot2::theme_void()

  # Use patchwork to combine plots with blank plot at top of influence column
  # This ensures horizontal alignment between stepwise panels and influence panels
  p_step + (blank_plot / p_influ) + patchwork::plot_layout(guides = "collect")
}

#' @title Combined Step and Influence Plot
#' @description Creates a combined plot showing both the stepwise index and term influence for a `gam_influence` object.
#' @param obj A `gam_influence` object containing calculated indices.
#' @return A patchwork plot combining the stepwise index and term influence plots.
#' @importFrom patchwork plot_layout
#' @export
plot_step_and_influence <- function(obj) {
  p_step <- plot_stepwise_index(obj)
  p_influ <- plot_term_influence(obj)

  # Use patchwork to combine plots side-by-side
  p_step + p_influ + patchwork::plot_layout(guides = "collect")
}

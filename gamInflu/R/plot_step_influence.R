#' Plot Step and Influence Together (ggplot)
#'
#' Side-by-side step and influence plots
#'
#' @param x gam_influence object
#' @param panels Logical, whether to use panel layout
#' @param ... Additional arguments
#' @return Combined ggplot object
#' @export
#' 
plot_step_influence <- function(x, panels = FALSE, ...) {
  
  if (!x$calculated) {
    stop("Must run calculate() first")
  }
  
  # Create individual plots
  step_plot <- plot_step(x, panels = panels, ...)
  influence_plot <- plot_influence(x, panels = panels, ...)
  
  # Combine plots
  combined_plot <- grid.arrange(
    step_plot, influence_plot,
    ncol = 2
  )
  
  return(combined_plot)
}

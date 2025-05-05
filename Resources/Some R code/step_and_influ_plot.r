#' Step and influence plots side by side
#'
#' @param obj The influence object
#' @param ... Additional arguments passed to plotting functions
#' @return A combined plot
#' @export
step_and_influ_plot <- function(obj, ...) {
  # Create both plots
  step_p <- step_plot(obj, panels = TRUE, ...)
  influ_p <- influ_plot(obj, panels = TRUE, ...)
  
  # Combine them side by side
  combined_plot <- step_p + influ_p + plot_layout(ncol = 2)
  return(combined_plot)
}


#' Plot CDI Distribution Panel (ggplot)
#'
#' @param x gam_influence object
#' @param variable Variable name for distribution
#' @param ... Additional plotting arguments
#' @return ggplot object
#' @keywords internal
#'
plot_cdi_distribution <- function(x, variable, ...) {
  # Get data for distribution plot
  focus_var_data <- x$model$model[[x$focus]]
  smooth_var_data <- x$data[[variable]]

  # Handle factor vs numeric variables differently
  if (is.factor(smooth_var_data) || is.character(smooth_var_data)) {
    return(plot_distribution_factor(focus_var_data, smooth_var_data, x$focus, variable))
  } else {
    return(plot_distribution_numeric(focus_var_data, smooth_var_data, x$focus, variable))
  }
}

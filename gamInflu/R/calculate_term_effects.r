#' Calculate Effects for a Specific Term
#'
#' @param model GAM model
#' @param term Term name
#' @param data Data used for fitting
#' @return List with term effects
#' @keywords internal
#' 
calculate_term_effects <- function(model, term, data) {
  # Check if it's a smooth term
  smooth_labels <- sapply(model$smooth, function(x) x$label)

  if (term %in% smooth_labels) {
    # Handle smooth term
    return(calculate_smooth_effects(model, term, data))
  } else {
    # Handle parametric term
    return(calculate_parametric_effects(model, term, data))
  }
}

#' Calculate Influence Statistics
#'
#' Perform influence calculations for the GAM model
#'
#' @param x A gam_influence object
#' @param ... Additional arguments
#' @return The gam_influence object with calculated statistics
#' @keywords internal
#' 
calculate_gam_influence <- function(x, ...) {
  if (x$calculated) {
    message("Calculations already performed.")
    return(x)
  }

  # Validate that focus term is parametric
  if (!is_parametric_term(x$model, x$focus)) {
    stop("Focus term must be a parametric term, not a smooth term")
  }

  # Expand by-variables into separate terms
  x$expanded_terms <- expand_by_variables(x$model)

  # Calculate effects for parametric focus term
  x$focus_effects <- calculate_parametric_focus_effects(x$model, x$focus, x$data)

  # Calculate model progression (step-wise addition)
  x$model_progression <- calculate_model_progression(x)

  # Calculate influence of smooth terms on parametric focus term
  x$influences <- calculate_smooth_influences(x)

  # Calculate goodness of fit statistics
  x$fit_stats <- calculate_fit_statistics(x)

  x$calculated <- TRUE
  return(x)
}

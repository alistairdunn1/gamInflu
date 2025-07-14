#' Calculate Smooth Term Influences on Parametric Focus Term
#'
#' @param x gam_influence object
#' @return List with influence calculations per focus level
#' @keywords internal
#'
calculate_smooth_influences <- function(x) {
  influences <- list()

  # Get predictions for all terms
  pred_terms <- predict(x$model, type = "terms")

  # Get focus term levels
  focus_levels <- names(x$focus_effects$relative_effects)

  for (term_name in names(x$expanded_terms)) {
    term_info <- x$expanded_terms[[term_name]]

    # Skip if this is related to the focus term itself
    if (is_focus_related(term_name, x$focus)) {
      next
    }

    # Calculate per-level influences (matching original Influ.r methodology)
    influences[[term_name]] <- calculate_per_level_influence(x, term_info, pred_terms, focus_levels)
  }

  return(influences)
}

#' Get Model Summary
#'
#' @param x gam_influence object
#' @return List with model summary information
#' @export
#'
get_model_summary <- function(x) {
  if (!x$calculated) {
    stop("Must run calculate() first")
  }

  return(list(
    focus_term = x$focus,
    n_smooth_terms = length(x$expanded_terms),
    model_class = class(x$model)[1],
    fit_statistics = x$fit_stats,
    r2 = x$model_progression
  ))
}

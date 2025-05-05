#' Extract effects for a term
#'
#' @param obj The influence object
#' @param model The model
#' @param term The term to extract effects for
#' @return Effects for the term
#' @export
get_effects <- function(obj, model = obj$model, term = obj$focus) {
  coeffs <- get_coeffs(obj, model, term)
  exp(coeffs - mean(coeffs))
}


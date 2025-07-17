#' @title Get Model Terms
#' @description Returns the terms used in a fitted model as a character vector.
#' @param obj A `gam_influence` object or a fitted model object.
#' @param full Logical; if TRUE, returns the full term expressions (e.g., s(term, ...)), otherwise returns variable names.
#' @return A character vector of model terms.
#' @export
get_terms <- function(obj, full = FALSE) {
  if ("gam_influence" %in% class(obj)) {
    if (full && !is.null(obj$model$formula)) {
      return(attr(terms(formula(obj$model)), "term.labels"))
    }
    return(obj$terms)
  } else {
    stop("Input must be a gam_influence object.")
  }
}

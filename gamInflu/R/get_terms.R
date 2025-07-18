#' @title Get Model Terms
#' @description Returns the terms used in a fitted model as a character vector.
#' @param obj A `gam_influence` object or a fitted model object.
#' @param full Logical; if TRUE, returns the full term expressions (e.g., s(term, ...)), otherwise returns variable names.
#' @return A character vector of model terms.
#' @export
get_terms <- function(obj, full = FALSE) {
  if ("gam_influence" %in% class(obj)) {
    res <- attr(terms(formula(obj$model)), "term.labels")
    if (full) {
      return(res)
    } else {
      # Extract variable names from the terms, excluding by variables and second terms in interactions
      return(unique(unlist(lapply(res, function(x) {
        vars <- all.vars(rlang::parse_expr(x))
        # Remove variables used in by= argument
        if (grepl("by\\s*=", x)) {
          by_var <- sub(".*by\\s*=\\s*([^,\\)]+).*", "\\1", x)
          vars <- setdiff(vars, by_var)
        }
        # Only keep the first variable (for interactions or smooths with multiple vars)
        vars[1]
      }))))
    }
    return(obj$terms)
  } else {
    stop("Input must be a gam_influence object.")
  }
}

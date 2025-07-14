#' Get terms from a gam_influence object
#' @param x A gam_influence object
#' @return Character vector of term names
#' @export
#' 
get_terms <- function(x) {
  if (!inherits(x, "gam_influence")) {
    stop("Object must be of class 'gam_influence'")
  }
  
  if (!x$calculated) {
    warning("Calculations not yet performed. Run calculate() first.")
    return(character(0))
  }
  
  if (is.null(x$expanded_terms)) {
    warning("No expanded_terms found. Available components: ", paste(names(x), collapse = ", "))
    return(character(0))
  }
  
  return(names(x$expanded_terms))
}

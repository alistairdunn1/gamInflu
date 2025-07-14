#' Get Per-Level Influences for a Specific Term
#'
#' Extract the per-level influence values for detailed analysis
#'
#' @param x gam_influence object
#' @param term Term name
#' @return Named vector of influences per focus level
#' @export
#' 
get_per_level_influences <- function(x, term) {
  
  if (!x$calculated) {
    stop("Must run calculate() first")
  }
  
  if (!term %in% names(x$influences)) {
    stop("Term '", term, "' not found in influences. Available terms: ",
         paste(names(x$influences), collapse = ", "))
  }
  
  return(x$influences[[term]]$per_level)
}

#' Reset Calculations
#'
#' Remove calculated results to allow recalculation with different parameters
#'
#' @param x gam_influence object
#' @return gam_influence object with calculations reset
#' @export
reset_calculations <- function(x) {
  
  if (!validate_gam_influence(x)) {
    stop("Invalid gam_influence object")
  }
  
  # Remove calculated components
  calculated_components <- c("focus_effects", "model_progression", 
                           "influences", "fit_stats", "expanded_terms")
  
  for (comp in calculated_components) {
    if (comp %in% names(x)) {
      x[[comp]] <- NULL
    }
  }
  
  x$calculated <- FALSE
  
  return(x)
}

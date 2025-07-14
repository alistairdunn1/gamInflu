#' Summary Method for gam_influence Objects
#'
#' @param object A gam_influence object
#' @param ... Additional arguments
#' @export
summary.gam_influence <- function(object, ...) {
  
  if (!object$calculated) {
    stop("Must run calculate() first")
  }
  
  cat("GAM Influence Analysis Summary\n")
  cat("=============================\n\n")
  
  # Print model info
  print(object)
  
  # Add more detailed summary information here
  cat("\nSmooth terms:\n")
  for (term in object$terms_info$smooth_terms) {
    info <- object$terms_info$smooth_info[[term]]
    cat("  ", term, " (edf:", round(info$edf, 2), ")\n")
  }
}

#' Print Method for gam_influence Objects
#'
#' @param x A gam_influence object
#' @param ... Additional arguments
#' @export
print.gam_influence <- function(x, ...) {
  cat("GAM Influence Analysis Object\n")
  cat("============================\n\n")
  cat("Focus term:", x$focus, "\n")
  cat("Response variable:", x$response, "\n")
  cat("Number of smooth terms:", x$terms_info$n_smooth, "\n")
  cat("Number of parametric terms:", x$terms_info$n_parametric, "\n")
  cat("Calculations performed:", ifelse(x$calculated, "Yes", "No"), "\n\n")

  if (x$calculated) {
    cat("Model fit statistics:\n")
    cat("  R-squared:", round(x$fit_stats$r_squared, 3), "\n")
    cat("  Deviance explained:", round(x$fit_stats$deviance_explained, 3), "\n")
    cat("  AIC:", round(x$fit_stats$aic, 2), "\n")
  }
}

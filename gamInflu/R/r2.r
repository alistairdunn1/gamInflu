#' @title Get Model Progression Statistics
#' @description Extracts the model progression summary table from a calculated `gam_influence` object.
#' This table includes statistics like AIC, R-squared, and deviance explained for each step of the model build.
#' @param x A `gam_influence` object that has been processed by `calculate_influence()`.
#' @return A data frame containing the model progression statistics.
#' @export
r2 <- function(x) {
  UseMethod("r2")
}

#' @title Model Progression Statistics for gam_influence
#' @description Method for extracting the model progression summary table from a calculated `gam_influence` object.
#' @param x A `gam_influence` object that has been processed by `calculate_influence()`.
#' @return A data frame containing the model progression statistics.
#' @export
r2.gam_influence <- function(x) {
  if (is.null(x$calculated$summary)) {
    stop("Model progression has not been calculated. Please run `calculate_influence()` first.", call. = FALSE)
  }
  return(x$calculated$summary)
}

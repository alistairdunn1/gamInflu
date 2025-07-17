#' Get Model Progression Statistics
#'
#' Extracts the model progression summary table from a calculated `gam_influence` object.
#' This table includes statistics like AIC, R-squared, and deviance explained
#' for each step of the model build.
#'
#' @param obj A `gam_influence` object that has been processed by `calculate_influence()`.
#' @param ... Unused.
#' @return A data frame containing the model progression statistics.
#' @export
#'
r2 <- function(obj, ...) {
  UseMethod("r2")
}

#' @export
#' @describeIn r2 Method for `gam_influence` class.
r2.gam_influence <- function(obj, ...) {
  if (is.null(obj$calculated$summary)) {
    stop("Model progression has not been calculated. Please run `calculate_influence()` first.", call. = FALSE)
  }
  return(obj$calculated$summary)
}

#' Summary Method for 'gam_influence' Objects
#'
#' Prints the summary table of term contributions, including changes in model
#' fit (R-squared) and the calculated influence metrics.
#'
#' @param object A `gam_influence` object processed by `calculate_influence()`.
#' @param ... Unused, for S3 consistency.
#' @return Prints a formatted summary table to the console.
#' @export
#'
summary.gam_influence <- function(object, ...) {
  if (is.null(object$calculated)) {
    stop("Calculations not performed. Please run `calculate_influence()` first.", call. = FALSE)
  }
  cat("Influence Analysis Summary\n")
  cat("==========================\n")
  cat("Model Formula:", deparse(formula(object$model)), "\n")
  cat("Focus Term:", object$focus, "\n\n")

  # Format the summary data frame for printing
  summary_to_print <- object$calculated$summary
  numeric_cols <- sapply(summary_to_print, is.numeric)
  summary_to_print[numeric_cols] <- lapply(summary_to_print[numeric_cols], round, 4)

  print(summary_to_print, row.names = FALSE)
}

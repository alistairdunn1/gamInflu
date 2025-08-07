#' @title Summary Method for gam_influence Objects
#' @description Prints the summary table of term contributions, including changes in model fit (R-squared) and the calculated influence metrics.
#' @param object A `gam_influence` object processed by `calculate_influence()`.
#' @param ... Additional arguments (not used).
#' @return Prints a formatted summary table to the console (invisible indices).
#' @importFrom stats formula
#' @export
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

  invisible(object$calculated$indices)
}

#' @title Print Method for Indices from gam_influence Object
#' @description Extracts the focus indices from a `gam_influence` object, including
#' unstandardised and standardised indices, confidence intervals, and coefficients of variation (CV).
#' @param x A `gam_influence` object that has been processed by `calculate_influence()`.
#' @param ... Additional arguments (not used).
#' @return A data frame containing the focus indices with descriptive column names.
#' #' @export
print.gam_influence <- function(x, ...) {
  if (is.null(x$calculated$indices)) {
    stop("No indices found. Please run calculate_influence() first.", call. = FALSE)
  }
  cat("Focus term indices from gam_influence Object\n")
  cat("=======================================\n")
  indices_to_print <- x$calculated$indices
  numeric_cols <- sapply(indices_to_print, is.numeric)
  indices_to_print[numeric_cols] <- lapply(indices_to_print[numeric_cols], round, 4)
  print(indices_to_print, row.names = FALSE)
  invisible(indices_to_print)
}

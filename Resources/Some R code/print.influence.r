#' Print method for influence objects
#'
#' @param x The influence object
#' @param ... Additional arguments
#' @export
print.influence <- function(x, ...) {
  cat("Influence object for model:\n")
  cat("  Response:", x$response, "\n")
  cat("  Focus term:", x$focus, "\n")
  cat("  Terms:", paste(x$terms, collapse=", "), "\n")
  
  # Print interaction information if available
  if (length(x$interactions) > 0) {
    cat("\nInteractions:\n")
    for (term in names(x$interactions)) {
      vars <- paste(x$interactions[[term]]$vars, collapse=", ")
      cat("  ", term, ":", vars, "\n")
    }
  }
  
  if (!is.null(x$summary)) {
    cat("\nSummary statistics:\n")
    print(x$summary)
  } else {
    cat("\nCall calc() to generate summary statistics and enable plotting.\n")
  }
}


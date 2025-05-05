# -----------------------------------------------------------------------------
# Print and Summary Methods
# -----------------------------------------------------------------------------

#' Print influence_gam Object
#'
#' @param x An object of class `influence_gam`.
#' @param ... Additional arguments (unused).
#' @method print influence_gam
#' @export
print.influence_gam <- function(x, ...) {
  cat("--- Influence Analysis for GAM ---\n")
  cat("Model Call: ", deparse(x$model$call), "\n")
  cat("Focus Term: ", x$focus, " (Variable: ", x$focus_var, ")\n")
  cat("Response:   ", x$response, "\n")
  cat("Family:     ", x$model$family$family, " (Link: ", x$model$family$link, ")\n")
  cat("Status:     ", if(x$calculated) "Calculated" else "Not Calculated (Run calculate_influence())", "\n")
  if (x$calculated && !is.null(x$summary)) {
      cat("\nSummary Table (First 6 Terms):\n")
      print(utils::head(x$summary)) # Use utils::head
  } else if (x$calculated) {
      cat("\nSummary table is missing.\n")
  }
  cat("----------------------------------\n")
  invisible(x)
}

#' Summarize influence_gam Object
#'
#' Provides a summary, including the model comparison table and snippets
#' of indices and influences.
#'
#' @param object An object of class `influence_gam`.
#' @param ... Additional arguments (unused).
#' @method summary influence_gam
#' @export
summary.influence_gam <- function(object, ...) {
  if (!object$calculated) {
    stop("Calculations not performed. Run 'calculate_influence()' before summarizing.")
  }
  cat("--- Summary of Influence Analysis ---\n")
  cat("Focus Term: ", object$focus, " (Variable: ", object$focus_var, ")\n")
  cat("Response:   ", object$response, "\n")
  cat("Family:     ", object$model$family$family, " (Link: ", object$model$family$link, ")\n")


  if(!is.null(object$summary)){
      cat("\nModel Comparison (Terms Added Sequentially):\n")
      # Round only numeric columns before printing
      summary_to_print <- object$summary
      numeric_cols <- sapply(summary_to_print, is.numeric)
      summary_to_print[numeric_cols] <- lapply(summary_to_print[numeric_cols], round, 3)
      print(summary_to_print)
  } else {
      cat("\nModel comparison summary table is missing.\n")
  }

  if(!is.null(object$indices)){
      cat("\nStandardised vs Unstandardised Index (First 6 Levels):\n")
      index_summary <- object$indices[, c("level", "unstan", "stan", "stanLower", "stanUpper")]
      numeric_idx_cols <- sapply(index_summary, is.numeric)
      index_summary[numeric_idx_cols] <- lapply(index_summary[numeric_idx_cols], round, 3)
      print(utils::head(index_summary))
  } else {
      cat("\nIndices table is missing.\n")
  }

  if(!is.null(object$influences) && ncol(object$influences) > 1){
      cat("\nInfluence of Terms (First 6 Levels, Link Scale):\n")
      influ_summary <- object$influences
      numeric_influ_cols <- sapply(influ_summary, is.numeric)
      influ_summary[numeric_influ_cols] <- lapply(influ_summary[numeric_influ_cols], round, 3)
      print(utils::head(influ_summary))
  } else {
       cat("\nInfluences table is missing or contains no terms.\n")
  }

  cat("--------------------------------------\n")
  invisible(object)
}
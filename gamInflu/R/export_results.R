#' Updated export_results function to include per-level influences
#'
#' @param x gam_influence object
#' @return Data frame with influence results including per-level data
#' @export
#' 
export_results <- function(x) {
  
  if (!x$calculated) {
    stop("Must run calculate() first")
  }
  
  # Prepare results data frame
  terms <- names(x$influences)
  
  if (length(terms) == 0) {
    warning("No influence terms found")
    return(data.frame())
  }
  
  results_df <- data.frame(
    term = terms,
    type = sapply(terms, function(t) x$expanded_terms[[t]]$type),
    stringsAsFactors = FALSE
  )
  
  # Add overall influence scores (following original Influ.r)
  results_df$overall_influence <- sapply(terms, function(t) {
    inf <- x$influences[[t]]$overall
    if (is.null(inf)) NA else inf
  })
  
  # Add trend influence scores (following original Influ.r)
  results_df$trend_influence <- sapply(terms, function(t) {
    inf <- x$influences[[t]]$trend
    if (is.null(inf)) NA else inf
  })
  
  # Add model progression info if available
  if (!is.null(x$model_progression)) {
    results_df$model_r_squared <- summary(x$model)$r.sq
  }
  
  return(results_df)
}

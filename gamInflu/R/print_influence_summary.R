#' Print detailed influence summary
#'
#' @param x gam_influence object
#' @export
print_influence_summary <- function(x) {
  
  if (!x$calculated) {
    stop("Must run calculate() first")
  }
  
  cat("GAM Influence Analysis - Detailed Summary\n")
  cat("========================================\n\n")
  
  cat("Focus term:", x$focus, "\n")
  cat("Focus levels:", paste(names(x$focus_effects$relative_effects), collapse = ", "), "\n\n")
  
  for (term_name in names(x$influences)) {
    cat("Term:", term_name, "\n")
    cat("  Type:", x$expanded_terms[[term_name]]$type, "\n")
    cat("  Overall influence:", round(x$influences[[term_name]]$overall, 4), "\n")
    cat("  Trend influence:", round(x$influences[[term_name]]$trend, 4), "\n")
    cat("  Per-level influences:\n")
    
    per_level <- x$influences[[term_name]]$per_level
    for (level in names(per_level)) {
      cat("    ", level, ": ", round(per_level[level], 4), "\n")
    }
    cat("\n")
  }
}
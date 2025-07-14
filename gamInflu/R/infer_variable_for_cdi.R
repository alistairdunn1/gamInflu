#' Infer Variable Name for CDI Plot
#'
#' @param x gam_influence object
#' @param term_info Term information
#' @return Variable name
#' @keywords internal
#' 
infer_variable_for_cdi <- function(x, term_info) {
  
  # Get the main variable from the smooth term
  if (length(term_info$variables) > 0) {
    # For tensor products, use the first variable
    return(term_info$variables[1])
  }
  
  # If we can't infer, try to extract from original label
  label <- term_info$original_label
  
  # Extract variable name from labels like "s(year)", "te(x,y)", etc.
  var_match <- regmatches(label, gregexpr("(?<=\\()([^,)]+)", label, perl = TRUE))[[1]]
  
  if (length(var_match) > 0) {
    return(var_match[1])
  }
  
  stop("Cannot infer variable name for term. Please specify 'variable' argument.")
}

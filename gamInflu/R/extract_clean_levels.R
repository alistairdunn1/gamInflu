#' Create clean labels for terms
#'
#' @param coef_names Coefficient names from the GAM model
#' @param focus_term Focus term name
#' @return Vector of cleaned up names
#' @keywords internal
#' 
extract_clean_levels <- function(coef_names, focus_term) {
  if (focus_term == "(Intercept)") {
    return("Intercept")
  }
  
  # Remove the focus term prefix from coefficient names
  clean_names <- gsub(paste0("^", focus_term), "", coef_names)
  
  # If no prefix was removed (i.e., the name didn't start with focus_term),
  # keep the original name
  clean_names[clean_names == coef_names] <- coef_names[clean_names == coef_names]
  
  # Handle cases where clean_names is empty string (base level)
  clean_names[clean_names == ""] <- paste0(focus_term, "_baseline")
  
  return(clean_names)
}

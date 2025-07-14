#' Validate GAM Influence Object
#'
#' @param x Object to validate
#' @return Logical indicating if valid
#' @keywords internal
validate_gam_influence <- function(x) {
  
  if (!inherits(x, "gam_influence")) {
    return(FALSE)
  }
  
  required_components <- c("model", "data", "response", "focus", "terms_info")
  
  if (!all(required_components %in% names(x))) {
    return(FALSE)
  }
  
  return(TRUE)
}

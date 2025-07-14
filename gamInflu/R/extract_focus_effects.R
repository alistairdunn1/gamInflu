#' Extract Focus Term Effects from Model
#'
#' @param model GAM model
#' @param focus_term Focus term name
#' @return Vector of relative effects
#' @keywords internal
#' 
extract_focus_effects <- function(model, focus_term) {
  
  coefs <- model$coefficients
  
  if (focus_term == "(Intercept)") {
    focus_coefs <- coefs["(Intercept)"]
  } else {
    focus_pattern <- paste0("^", focus_term)
    focus_indices <- grep(focus_pattern, names(coefs))
    focus_coefs <- coefs[focus_indices]
  }
  
  if (length(focus_coefs) == 0) {
    return(NULL)
  }
  
  # Calculate relative effects
  if (length(focus_coefs) > 1) {
    base_effect <- mean(focus_coefs)
    relative_effects <- exp(focus_coefs - base_effect)
  } else {
    relative_effects <- exp(focus_coefs)
  }
  
  return(relative_effects)
}

#' calculate_parametric_focus_effects
#'
#' Calculate influence of a term for each level of the focus term,
#'
#' @param x gam_influence object
#' @param focus_term Focus term
#' @param data Data used for the model
#' @return List with influence metrics per focus level
#' @keywords internal
#' 
calculate_parametric_focus_effects <- function(model, focus_term, data) {
  
  # Get coefficients for focus term
  coefs <- model$coefficients
  
  # Handle different types of parametric terms
  if (focus_term == "(Intercept)") {
    focus_coefs <- coefs["(Intercept)"]
    focus_indices <- 1
  } else {
    focus_pattern <- paste0("^", focus_term)
    focus_indices <- grep(focus_pattern, names(coefs))
    focus_coefs <- coefs[focus_indices]
  }
  
  if (length(focus_coefs) == 0) {
    stop("Focus term '", focus_term, "' not found in model coefficients")
  }
  
  # Get standard errors
  vcov_mat <- vcov(model)
  focus_ses <- sqrt(diag(vcov_mat)[focus_indices])
  
  # Calculate relative effects (similar to original package)
  if (length(focus_coefs) > 1) {
    # Factor variable - calculate relative to mean
    base_effect <- mean(focus_coefs)
    relative_effects <- exp(focus_coefs - base_effect)
    relative_lower <- exp(focus_coefs - base_effect - 2 * focus_ses)
    relative_upper <- exp(focus_coefs - base_effect + 2 * focus_ses)
  } else {
    # Single coefficient
    relative_effects <- exp(focus_coefs)
    relative_lower <- exp(focus_coefs - 2 * focus_ses)
    relative_upper <- exp(focus_coefs + 2 * focus_ses)
  }
  
  # Extract clean level names
  clean_levels <- extract_clean_levels(names(focus_coefs), focus_term)
  names(relative_effects) <- clean_levels
  names(relative_lower) <- clean_levels
  names(relative_upper) <- clean_levels
  
  return(list(
    coefficients = focus_coefs,
    standard_errors = focus_ses,
    relative_effects = relative_effects,
    confidence_lower = relative_lower,
    confidence_upper = relative_upper,
    term_levels = clean_levels  # Now contains clean level names
  ))
}

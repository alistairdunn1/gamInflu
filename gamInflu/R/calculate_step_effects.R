#' Calculate Focus Term Effects at Each Step
#'
#' @param x gam_influence object
#' @return List with effects at each step
#' @keywords internal
#' 
calculate_step_effects <- function(x) {
  
  original_formula <- formula(x$model)
  all_terms <- attr(terms(x$model), "term.labels")
  
  # Separate parametric and smooth terms
  param_terms <- all_terms[!grepl("s\\(|te\\(|ti\\(", all_terms)]
  smooth_terms <- all_terms[grepl("s\\(|te\\(|ti\\(", all_terms)]
  
  step_effects <- list()
  current_terms <- character(0)
  
  # Add parametric terms first
  for (i in seq_along(param_terms)) {
    current_terms <- c(current_terms, param_terms[i])
    
    if (param_terms[i] == x$focus || any(grepl(x$focus, current_terms))) {
      # Focus term is now in the model
      new_formula_str <- paste("~", paste(current_terms, collapse = " + "))
      new_formula <- update(original_formula, formula(new_formula_str))
      
      tryCatch({
        model_step <- gam(new_formula, data = x$data, family = x$model$family)
        focus_effects <- extract_focus_effects(model_step, x$focus)
        step_effects[[length(step_effects) + 1]] <- list(
          step = length(step_effects) + 1,
          term_added = param_terms[i],
          effects = focus_effects
        )
      }, error = function(e) {
        warning("Could not extract effects for step: ", param_terms[i])
      })
    }
  }
  
  # Add smooth terms
  for (i in seq_along(smooth_terms)) {
    current_terms <- c(current_terms, smooth_terms[i])
    new_formula_str <- paste("~", paste(current_terms, collapse = " + "))
    new_formula <- update(original_formula, formula(new_formula_str))
    
    tryCatch({
      model_step <- gam(new_formula, data = x$data, family = x$model$family)
      
      if (any(grepl(x$focus, current_terms))) {
        focus_effects <- extract_focus_effects(model_step, x$focus)
        step_effects[[length(step_effects) + 1]] <- list(
          step = length(step_effects) + 1,
          term_added = smooth_terms[i],
          effects = focus_effects
        )
      }
    }, error = function(e) {
      warning("Could not extract effects for step: ", smooth_terms[i])
    })
  }
  
  return(step_effects)
}

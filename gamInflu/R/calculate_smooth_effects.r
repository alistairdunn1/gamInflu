#' Calculate Smooth Term Effects
#'
#' @param model GAM model
#' @param term Smooth term label
#' @param data Data
#' @return List with smooth effects
#' @keywords internal
#' 
calculate_smooth_effects <- function(model, term, data) {
  # Find the smooth term
  smooth_idx <- which(sapply(model$smooth, function(x) x$label) == term)
  smooth_obj <- model$smooth[[smooth_idx]]

  # Get predictions for this smooth term
  pred_terms <- predict(model, type = "terms", se.fit = TRUE)

  # Extract this specific term's predictions
  term_pred <- pred_terms$fit[, smooth_idx]
  term_se <- pred_terms$se.fit[, smooth_idx]

  # Get the main variable for this smooth
  main_var <- smooth_obj$term[1] # First variable for multidimensional smooths

  # Calculate effects by levels of main variable
  if (smooth_obj$by != "NA") {
    # Handle by-variables
    by_var <- smooth_obj$by
    effects_data <- data.frame(
      main_var = data[[main_var]],
      by_var = data[[by_var]],
      effect = term_pred,
      se = term_se
    )
  } else {
    effects_data <- data.frame(
      main_var = data[[main_var]],
      effect = term_pred,
      se = term_se
    )
  }

  return(list(
    type = "smooth",
    smooth_info = smooth_obj,
    effects_data = effects_data,
    edf = model$edf[smooth_idx]
  ))
}

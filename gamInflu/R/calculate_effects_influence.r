#' Calculate Effects-based Influence
#'
#' @param x gam_influence object
#' @param term_info Information about the term
#' @return Numeric influence score
#' @keywords internal
#' 
calculate_effects_influence <- function(x, term_info) {
  # Get predictions for this smooth term
  pred_terms <- predict(x$model, type = "terms")
  smooth_idx <- term_info$smooth_index

  if (term_info$type == "by_smooth") {
    # For by-variables, subset to relevant observations
    by_var <- term_info$by_variable
    by_level <- term_info$by_level
    relevant_obs <- x$model$model[[by_var]] == by_level
    term_effects <- pred_terms[relevant_obs, smooth_idx]
  } else {
    term_effects <- pred_terms[, smooth_idx]
  }

  # For tensor products, combine the influence across dimensions
  if (is_tensor_product(term_info$original_label)) {
    influence_score <- calculate_tensor_influence(term_effects, term_info)
  } else {
    # Regular smooth: calculate influence as root mean square of effects
    influence_score <- sqrt(mean(term_effects^2))
  }

  return(influence_score)
}

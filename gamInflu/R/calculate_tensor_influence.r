#' Calculate Tensor Product Influence
#'
#' For tensor products, calculate combined influence measure
#'
#' @param term_effects Vector of term effects
#' @param term_info Term information
#' @return Combined influence score
#' @keywords internal
#' 
calculate_tensor_influence <- function(term_effects, term_info) {
  # For tensor products, we want to capture the combined effect
  # across all dimensions. Use RMS as a summary measure.

  n_vars <- length(term_info$variables)

  if (n_vars == 1) {
    # Not actually a tensor product
    return(sqrt(mean(term_effects^2)))
  }

  # For multi-dimensional tensor products:
  # Weight by the number of dimensions and calculate combined influence
  combined_influence <- sqrt(mean(term_effects^2)) * sqrt(n_vars)

  return(combined_influence)
}

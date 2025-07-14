#' Prepare Influence Data for ggplot (CORRECTED)
#'
#' Now properly extracts per-level influences instead of repeating single values
#'
#' @param x gam_influence object
#' @return Data frame with influence scores per level
#' @keywords internal
prepare_influence_data <- function(x) {
  
  terms <- names(x$influences)
  
  if (length(terms) == 0) {
    warning("No influence terms found")
    return(data.frame())
  }
  
  # Prepare data frame
  influence_data <- data.frame()
  focus_levels <- names(x$focus_effects$relative_effects)
  
  for (i in seq_along(terms)) {
    term_name <- terms[i]
    term_inf <- x$influences[[term_name]]
    
    if (!is.null(term_inf$per_level)) {
      # Use the per-level influences (this is the correction!)
      per_level_influences <- term_inf$per_level
      
      # Ensure we have influences for all focus levels
      influence_scores <- rep(NA, length(focus_levels))
      names(influence_scores) <- focus_levels
      
      # Fill in available influences
      for (level in names(per_level_influences)) {
        if (level %in% focus_levels) {
          influence_scores[level] <- per_level_influences[level]
        }
      }
      
      # Replace any remaining NAs with 1 (no influence)
      influence_scores[is.na(influence_scores)] <- 1
      
    } else {
      # Fallback: no influence
      influence_scores <- rep(1, length(focus_levels))
      names(influence_scores) <- focus_levels
    }
    
    # Create data for each focus level
    term_data <- data.frame(
      term = term_name,
      level = focus_levels,
      level_number = seq_along(focus_levels),
      influence = as.numeric(influence_scores)
    )
    
    influence_data <- rbind(influence_data, term_data)
  }
  
  return(influence_data)
}

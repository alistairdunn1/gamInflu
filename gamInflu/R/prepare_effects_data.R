#' Prepare Effects Data for ggplot
#'
#' @param x gam_influence object
#' @param terms Terms to include
#' @return Data frame for ggplot
#' @keywords internal
prepare_effects_data <- function(x, terms) {
  
  effects_data <- data.frame()
  
  # Get predictions for all terms
  pred_terms <- predict(x$model, type = "terms", se.fit = TRUE)
  
  for (term in terms) {
    term_info <- x$expanded_terms[[term]]
    smooth_idx <- term_info$smooth_index
    
    # Get the main variable
    main_var <- term_info$variables[1]
    var_data <- x$data[[main_var]]
    
    if (term_info$type == "by_smooth") {
      # Handle by-variables
      by_var <- term_info$by_variable
      by_level <- term_info$by_level
      
      # Subset to relevant observations
      relevant_obs <- x$model$model[[by_var]] == by_level
      term_effects <- pred_terms$fit[relevant_obs, smooth_idx]
      term_se <- pred_terms$se.fit[relevant_obs, smooth_idx]
      var_subset <- var_data[relevant_obs]
      
      term_label <- paste0(term_info$original_label, " (", by_level, ")")
    } else {
      term_effects <- pred_terms$fit[, smooth_idx]
      term_se <- pred_terms$se.fit[, smooth_idx]
      var_subset <- var_data
      term_label <- term_info$original_label
    }
    
    # Create data frame for this term
    term_data <- data.frame(
      term = term,
      term_label = term_label,
      var_name = main_var,
      var_value = if (is.numeric(var_subset)) var_subset else as.numeric(as.factor(var_subset)),
      var_type = if (is.numeric(var_subset)) "numeric" else "factor",
      effect = term_effects,
      se = term_se
    )
    
    # For factors, add factor level names
    if (!is.numeric(var_subset)) {
      term_data$var_label <- as.character(var_subset)
    } else {
      term_data$var_label <- as.character(var_subset)
    }
    
    effects_data <- rbind(effects_data, term_data)
  }
  
  return(effects_data)
}

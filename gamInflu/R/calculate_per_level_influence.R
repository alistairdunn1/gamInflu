#' Calculate Per-Level Influence (CORRECTED)
#'
#' Calculate influence of a term for each level of the focus term,
#' following the original Influ.r methodology
#'
#' @param x gam_influence object
#' @param term_info Information about the term
#' @param pred_terms Predictions from predict(model, type = "terms")
#' @param focus_levels Levels of the focus term
#' @return List with influence metrics per focus level
#' @keywords internal
#' 
calculate_per_level_influence <- function(x, term_info, pred_terms, focus_levels) {
  
  smooth_idx <- term_info$smooth_index
  
  # Get the focus term data from the model
  focus_data <- x$model$model[[x$focus]]
  
  if (term_info$type == "by_smooth") {
    # For by-variables, subset to relevant observations
    by_var <- term_info$by_variable
    by_level <- term_info$by_level
    relevant_obs <- x$model$model[[by_var]] == by_level
    
    term_effects <- pred_terms[relevant_obs, smooth_idx]
    focus_subset <- focus_data[relevant_obs]
  } else {
    term_effects <- pred_terms[, smooth_idx]
    focus_subset <- focus_data
  }
  
  # Aggregate term effects by focus term levels
  # This matches the original: aggregate(list(value = preds[, paste("fit", term, sep = ".")]), 
  #                                      list(level = preds[, focus]), mean)
  influence_by_level <- aggregate(
    list(value = term_effects),
    list(level = focus_subset),
    mean
  )
  
  # Ensure we have values for all focus levels
  all_levels_df <- data.frame(level = levels(focus_data))
  influence_by_level <- merge(all_levels_df, influence_by_level, by = "level", all.x = TRUE)
  
  # Replace any NAs with 0 (for levels where this term has no effect)
  influence_by_level$value[is.na(influence_by_level$value)] <- 0
  
  # Calculate influence metrics following original Influ.r formulas
  values <- influence_by_level$value
  
  # Original overall influence: exp(mean(abs(value))) - 1
  overall_influence <- exp(mean(abs(values))) - 1
  
  # Original trend influence: exp(cov(1:length(value), value) / var(1:length(value))) - 1
  if (length(values) > 1 && var(1:length(values)) > 0) {
    trend_influence <- exp(cov(1:length(values), values) / var(1:length(values))) - 1
  } else {
    trend_influence <- 0
  }
  
  # Per-level relative influences (standardized to show relative effect)
  # Convert to exponential scale for interpretation (matching original plots)
  per_level_influences <- exp(values)
  names(per_level_influences) <- influence_by_level$level
  
  return(list(
    per_level = per_level_influences,
    overall = overall_influence,
    trend = trend_influence,
    raw_values = values
  ))
}

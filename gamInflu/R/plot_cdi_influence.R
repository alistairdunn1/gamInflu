#' Plot CDI Influence Panel (ggplot)
#'
#' @param x gam_influence object
#' @param term Term name
#' @param ... Additional plotting arguments
#' @return ggplot object
#' @keywords internal
#' 
plot_cdi_influence <- function(x, term, ...) {
  
  # Get influence scores for this term
  term_influences <- x$influences[[term]]
  
  # Use appropriate influence measure
  if (!is.null(term_influences$effects_influence)) {
    influence_scores <- term_influences$effects_influence
  } else {
    influence_scores <- NA
  }
  
  # Create data frame
  focus_levels <- names(x$focus_effects$relative_effects)
  n_focus_levels <- length(focus_levels)
  
  # Create influence scores for each focus level (placeholder approach)
  if (length(influence_scores) == 1) {
    influence_scores <- rep(influence_scores, n_focus_levels)
  }
  
  # Standardize influence to have mean 1
  mean_influence <- mean(influence_scores, na.rm = TRUE)
  if (mean_influence > 0) {
    influence_scores <- influence_scores / mean_influence
  }
  
  infl_data <- data.frame(
    level = focus_levels,
    level_number = seq_len(n_focus_levels),
    influence = influence_scores
  )
  
  p <- ggplot(infl_data, aes(x = influence, y = level_number)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
    geom_line(color = "darkred") +
    geom_point(color = "darkred") +
    scale_y_continuous(breaks = infl_data$level_number, 
                       labels = infl_data$level,
                       position = "right") +
    labs(
      x = "Standardized Influence",
      y = "",
      #subtitle = "Mean = 1"
    )
  
  return(p)
}

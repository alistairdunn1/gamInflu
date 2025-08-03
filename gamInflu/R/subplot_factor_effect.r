#' @title Subplot for Factor Effects
#' @description Internal function to plot predicted effects for factor terms.
#' @param obj A `gam_influence` object.
#' @param t The term name.
#' @param preds_df Data frame of predictions.
#' @param se_df Data frame of standard errors.
#' @param term_vars Variables in the term.
#' @param se_col Standard error column(s).
#' @param islog Logical; if TRUE, exponentiate effects.
#' @noRd
subplot_factor_effect <- function(obj, t, term_vars, cdi = FALSE) {
  message("Plotting factor effects for term: ", t)
  islog <- isTRUE(obj$islog)
  preds_df <- obj$calculated$predictions
  se_df <- obj$calculated$prediction_se

  effect <- preds_df[[t]]
  se <- se_df[[t]]
  if (islog) {
    lower <- exp(effect - 1.96 * se)
    upper <- exp(effect + 1.96 * se)
    effect <- exp(effect)
    ylim <- c(0, NA)
  } else {
    lower <- effect - 1.96 * se
    upper <- effect + 1.96 * se
    ylim <- c(NA, NA)
  }

  # Handle potential dimension mismatch between data and predictions (subset analysis)
  # Use only the levels that correspond to the predictions
  if (length(effect) != nrow(obj$data)) {
    # This indicates subset analysis was performed
    # The predictions should correspond to the factor levels present in the subset
    # Get the unique factor levels from the subset data that was used for predictions
    message("Detected subset analysis: adjusting factor levels to match predictions")

    # Get all unique levels from the current data (which is the full data after restore)
    all_levels <- unique(obj$data[[term_vars[1]]][!is.na(obj$data[[term_vars[1]]])])

    # If we have the right number of levels, use them
    if (length(all_levels) == length(effect)) {
      df <- data.frame(level = all_levels, effect = effect, lower = lower, upper = upper)
    } else {
      # Create factor levels based on the length of predictions
      # This assumes the predictions are in the same order as sorted factor levels
      sorted_levels <- sort(all_levels)
      if (length(sorted_levels) >= length(effect)) {
        used_levels <- sorted_levels[seq_len(length(effect))]
        df <- data.frame(level = used_levels, effect = effect, lower = lower, upper = upper)
      } else {
        # Fallback: use numeric indices with warning
        df <- data.frame(level = paste("Level", seq_along(effect)), effect = effect, lower = lower, upper = upper)
        warning("Could not match factor levels to predictions. Using generic labels.")
      }
    }
  } else {
    # Normal case - dimensions match
    df <- data.frame(level = obj$data[[term_vars[1]]], effect = effect, lower = lower, upper = upper)
  }

  df <- df[!is.na(df$level), ]

  p_coef <- ggplot2::ggplot(df, ggplot2::aes(x = level, y = effect)) +
    ggplot2::geom_point(size = 3, colour = "royalblue") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), colour = "royalblue", alpha = 0.5, width = 0.2, na.rm = TRUE) +
    ggplot2::ylim(ylim) +
    ggplot2::labs(y = "Partial effect")
  if (cdi) {
    p_coef <- p_coef + ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank()
    )
  } else {
    p_coef <- p_coef +
      ggplot2::xlab(term_vars[1])
  }
  return(p_coef)
}

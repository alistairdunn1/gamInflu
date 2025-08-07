#' @title Subplot for Continuous Effects
#' @description Internal function to plot predicted effects for continuous terms.
#' @param obj A `gam_influence` object.
#' @param t The term name.
#' @param term_vars Variables in the term.
#' @param cdi Logical indicating if the plot is for CDI (Cumulative Distribution Influence).
#' @return A ggplot object showing the continuous effects.
#' @noRd
subplot_continuous_effect <- function(obj, t, term_vars, cdi = FALSE) {
  message("Plotting continuous effects for term: ", t)
  islog <- isTRUE(obj$islog)
  preds_df <- obj$calculated$predictions
  se_df <- obj$calculated$prediction_se

  matching_cols <- names(preds_df)[
    vapply(names(preds_df), function(nm) {
      grepl(term_vars[1], nm)
    }, logical(1))
  ]
  preds_df$effect <- preds_df[[matching_cols[1]]]
  preds_df$lower <- preds_df$effect - 1.96 * se_df[[matching_cols[1]]]
  preds_df$upper <- preds_df$effect + 1.96 * se_df[[matching_cols[1]]]

  if (islog) {
    preds_df$effect <- exp(preds_df$effect)
    preds_df$lower <- exp(preds_df$lower)
    preds_df$upper <- exp(preds_df$upper)
    ylim <- c(0, NA_real_)
  } else {
    ylim <- c(NA_real_, NA_real_)
  }

  # Handle potential dimension mismatch between data and predictions (subset analysis)
  if (length(preds_df$effect) != nrow(obj$data)) {
    # This indicates subset analysis was performed
    message("Detected subset analysis: using prediction data length for continuous effects")

    # For continuous variables, we need to ensure the x values match the predictions
    # Use the first n values where n = length of predictions
    x_values <- obj$data[[term_vars[1]]][seq_len(length(preds_df$effect))]
    df <- data.frame(
      x = x_values,
      effect = preds_df$effect,
      lower = preds_df$lower,
      upper = preds_df$upper
    )
  } else {
    # Normal case - dimensions match
    df <- data.frame(
      x = obj$data[[term_vars[1]]],
      effect = preds_df$effect,
      lower = preds_df$lower,
      upper = preds_df$upper
    )
  }
  p_coef <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = effect)) +
    ggplot2::geom_line(colour = "royalblue") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "royalblue") +
    ggplot2::labs(y = "Partial effect") +
    ggplot2::ylim(ylim)
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

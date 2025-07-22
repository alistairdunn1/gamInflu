#' @title Subplot for Continuous Effects
#' @description Internal function to plot predicted effects for continuous terms.
#' @param obj A `gam_influence` object.
#' @param t The term name.
#' @param term_vars Variables in the term.
#' @param cdi Logical indicating if the plot is for a continuous derivative influence (CDI).
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
    ylim <- c(0, NA)
  } else {
    ylim <- c(NA, NA)
  }
  df <- data.frame(
    x = obj$data[[term_vars[1]]], effect = preds_df$effect, lower = preds_df$lower, upper = preds_df$upper
  )
  p_coef <- ggplot(df, aes(x = x, y = effect)) +
    geom_line(colour = "royalblue") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "royalblue") +
    labs(y = "Partial effect")
  if (cdi) {
    p_coef <- p_coef + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), legend.title = element_blank())
  } else {
    p_coef <- p_coef +
      xlab(term_vars[1])
  }
  return(p_coef)
}

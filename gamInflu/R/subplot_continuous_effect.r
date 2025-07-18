#' @title Subplot for Continuous Effects
#' @description Internal function to plot predicted effects for continuous terms.
#' @param obj A `gam_influence` object.
#' @param t The term name.
#' @param preds_df Data frame of predictions.
#' @param se_df Data frame of standard errors.
#' @param term_vars Variables in the term.
#' @param se_col Standard error column(s).
#' @param islog Logical; if TRUE, exponentiate effects.
#' @noRd
subplot_continuous_effect <- function(obj, t, term_vars, cdi = FALSE) {
  islog <- isTRUE(obj$islog)
  obj_data <- obj$data
  preds_df <- obj$calculated$predictions
  se_df <- obj$calculated$prediction_se

  effect <- preds_df[[t]]
  se <- if (length(se_col) > 0) se_df[[se_col[1]]] else NA
  ymin <- effect - 1.96 * se
  ymax <- effect + 1.96 * se
  if (islog) {
    effect <- exp(effect)
    ymin <- exp(ymin)
    ymax <- exp(ymax)
  }
  df <- data.frame(
    x = preds_df[[term_vars[1]]], effect = effect, ymin = ymin, ymax = ymax
  )
  ggplot(df, aes(x = x, y = effect)) +
    geom_line(colour = "royalblue") +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, fill = "royalblue") +
    labs(x = t, y = "Effect")
}

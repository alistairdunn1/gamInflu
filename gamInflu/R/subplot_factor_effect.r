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
  df <- data.frame(level = obj$data[[term_vars[1]]], effect = effect, lower = lower, upper = upper)

  p_coef <- ggplot(df, aes(x = level, y = effect)) +
    geom_point(size = 3, colour = "royalblue") +
    geom_errorbar(aes(ymin = lower, ymax = upper), colour = "royalblue", alpha = 0.5, width = 0.2, na.rm = TRUE) +
    ylim(ylim) +
    labs(y = "Partial effect")
  if (cdi) {
    p_coef <- p_coef + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), legend.title = element_blank())
  } else {
    p_coef <- p_coef +
      xlab(term_vars[1])
  }
  return(p_coef)
}

#' @title Subplot for Focus Effects
#' @description Internal function to plot predicted effects for focus terms.
#' @param obj A `gam_influence` object.
#' @param t The term name.
#' @param preds_df Data frame of predictions.
#' @param se_df Data frame of standard errors.
#' @param term_vars Variables in the term.
#' @param se_col Standard error column(s).
#' @param islog Logical; if TRUE, exponentiate effects.
#' @noRd
subplot_focus_effect <- function(obj, t, term_vars, cdi = FALSE) {
  if (t != obj$focus) {
    message(paste(t, " is not a focus term", sep = ""))
    return(NULL)
  }

  islog <- isTRUE(obj$islog)
  obj_data <- obj$data
  by_var <- sub(".*by\\s*=\\s*([^,\\)]+).*", "\\1", t)
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

  # Convert df$level to numeric if it is numeric-like (e.g., factor or character with numeric values)
  if (!is.numeric(df$level) && all(!is.na(as.numeric(as.character(df$level))))) {
    df$level <- as.numeric(as.character(df$level))
  }

  p_coef <- ggplot(df, aes(x = level, y = effect)) +
    geom_point(colour = "royalblue") +
    geom_line(colour = "royalblue") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "royalblue", alpha = 0.2) +
    ylim(ylim) +
    labs(y = "Index")
  if (cdi) {
    p_coef <- p_coef + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), legend.title = element_blank())
  } else {
    p_coef <- p_coef +
      xlab(term_vars[1])
  }
  return(p_coef)
}

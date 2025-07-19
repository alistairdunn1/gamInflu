#' @title Subplot for Random Effects
#' @description Internal function to plot predicted effects for random effect terms (e.g., s(term, bs="re")).
#' @param obj A `gam_influence` object.
#' @param t The term name.
#' @param type Plot type: "point", "bar", or "violin".
#' @param preds_df Data frame of predictions.
#' @param se_df Data frame of standard errors.
#' @param term_vars Variables in the term.
#' @param se_col Standard error column(s).
#' @param islog Logical; if TRUE, exponentiate effects.
#' @noRd
subplot_random_effect <- function(obj, t, term_vars, type = "point", cdi = FALSE) {
  islog <- isTRUE(obj$islog)
  obj_data <- obj$data
  by_var <- sub(".*by\\s*=\\s*([^,\\)]+).*", "\\1", t)
  terms_vec <- get_terms(obj, full = TRUE)
  term_list <- if (is.null(term)) terms_vec else term

  preds_df <- obj$calculated$predictions
  se_df <- obj$calculated$prediction_se

  matching_cols <- names(preds_df)[
    vapply(names(preds_df), function(nm) {
      grepl(term_vars[1], nm)
    }, logical(1))
  ]

  levels <- obj$data[[term_vars[1]]]
  effect <- preds_df[[matching_cols[1]]]
  se <- se_df[[matching_cols[1]]]

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

  if (type == "violin") {
    p_coef <- ggplot(df, aes(x = level, y = effect)) +
      geom_violin(fill = "royalblue") +
      labs(y = "Random Effect")
  } else if (type == "bar") {
    p_coef <- ggplot(df, aes(x = level, y = effect)) +
      geom_bar(stat = "identity", fill = "royalblue") +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, na.rm = TRUE) +
      labs(y = "Random Effect")
  } else if (type == "point") {
    p_coef <- ggplot(df, aes(x = level, y = effect)) +
      geom_point(size = 3, colour = "royalblue") +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, na.rm = TRUE) +
      labs(y = "Random Effect")
  } else {
    message("Unsupported type for random effect plot. Use 'point', 'bar', or 'violin'.")
    p_coef <- plot_spacer()
  }
  if (cdi) {
    p_coef <- p_coef + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), legend.title = element_blank())
  } else {
    p_coef <- p_coef +
      xlab(term_vars[1])
  }
  return(p_coef)
}

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
subplot_random_effect <- function(obj, t, type, term_vars, cdi = FALSE) {
  islog <- isTRUE(obj$islog)
  obj_data <- obj$data
  by_var <- sub(".*by\\s*=\\s*([^,\\)]+).*", "\\1", t)
  main_var <- term_vars[1]
  preds_df <- obj$calculated$predictions
  se_df <- obj$calculated$prediction_se

  if (!(t %in% colnames(preds_df)) || !(term_vars[1] %in% colnames(preds_df))) {
    warning(paste("Random effect term", t, "not found in predictions."))
    return(ggplot())
  }
  re_levels <- preds_df[[term_vars[1]]]
  re_effect <- preds_df[[t]]
  re_se <- if (length(se_col) > 0 && se_col[1] %in% colnames(se_df)) se_df[[se_col[1]]] else NA
  if (islog) {
    re_effect <- exp(re_effect)
    re_se <- ifelse(!is.na(re_se), exp(re_se), re_se)
  }
  df <- data.frame(level = re_levels, effect = re_effect, se = re_se)
  if (type == "violin") {
    ggplot(df, aes(x = level, y = effect)) +
      geom_violin(fill = "royalblue") +
      labs(title = t, x = term_vars[1], y = "Random Effect")
  } else if (type == "bar") {
    ggplot(df, aes(x = level, y = effect)) +
      geom_bar(stat = "identity", fill = "royalblue") +
      geom_errorbar(aes(ymin = effect - 1.96 * se, ymax = effect + 1.96 * se), width = 0.2, na.rm = TRUE) +
      labs(title = t, x = term_vars[1], y = "Random Effect")
  } else {
    ggplot(df, aes(x = level, y = effect)) +
      geom_point(size = 3, colour = "royalblue") +
      geom_errorbar(aes(ymin = effect - 1.96 * se, ymax = effect + 1.96 * se), width = 0.2, na.rm = TRUE) +
      labs(x = t, y = "Random Effect")
  }
}

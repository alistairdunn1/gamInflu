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
  preds_df <- obj$calculated$predictions
  se_df <- obj$calculated$prediction_se

  # Find columns in preds_df that match the random effect term
  matching_cols <- names(preds_df)[
    vapply(names(preds_df), function(nm) {
      grepl(term_vars[1], nm) && grepl(t, nm)
    }, logical(1))
  ]

  # Defensive: If no matching columns, return blank plot
  if (length(matching_cols) == 0) {
    warning("No matching columns found for random effect term.")
    return(plot_spacer())
  }
  if (length(matching_cols) > 1) {
    warning("Multiple matching columns found for random effect; using the first one.")
  }

  # Extract levels, effect, and se
  levels <- preds_df[[term_vars[1]]]
  effect <- preds_df[[matching_cols[1]]]
  se <- if (matching_cols[1] %in% names(se_df)) se_df[[matching_cols[1]]] else NA

  # Defensive: If effect is empty, return blank plot
  if (length(effect) == 0) {
    return(plot_spacer())
  }

  # Calculate confidence intervals and transform if log
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

  df <- data.frame(level = levels, effect = effect, lower = lower, upper = upper)

  # Ensure level is a factor for plotting
  if (!is.factor(df$level)) df$level <- factor(df$level)

  # Plot according to type
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

  # Remove axis labels for CDI plots
  if (cdi) {
    p_coef <- p_coef +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank()
      )
  } else {
    p_coef <- p_coef + xlab(term_vars[1])
  }
  return(p_coef)
}

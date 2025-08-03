#' @title Subplot for By-Variable Effects
#' @description Internal function to plot predicted effects for by-variable terms (e.g., s(var1, by=var2)).
#' @param obj A `gam_influence` object.
#' @param t The term name.
#' @param preds_df Data frame of predictions.
#' @param se_df Data frame of standard errors.
#' @param term_vars Variables in the term.
#' @param se_col Standard error column(s).
#' @param islog Logical; if TRUE, exponentiate effects.
#' @noRd
subplot_by_variable <- function(obj, t, term_vars, cdi = FALSE) {
  message("Plotting by-variable effects for term: ", t)
  islog <- isTRUE(obj$islog)
  by_var <- sub(".*by\\s*=\\s*([^,\\)]+).*", "\\1", t)
  preds_df <- obj$calculated$predictions
  se_df <- obj$calculated$prediction_se

  matching_cols <- names(preds_df)[
    vapply(names(preds_df), function(nm) {
      grepl(term_vars[1], nm) && grepl(by_var, nm)
    }, logical(1))
  ]

  # Handle potential dimension mismatch between data and predictions (subset analysis)
  if (nrow(preds_df) != nrow(obj$data)) {
    message("Detected subset analysis: adjusting by-variable data to match predictions")

    # Use only the first n rows of data where n = number of prediction rows
    n_pred <- nrow(preds_df)
    by_var_data <- obj$data[[by_var]][seq_len(n_pred)]
    this_term_data <- obj$data[[term_vars[1]]][seq_len(n_pred)]

    pred_long <- cbind(preds_df, by_var = by_var_data, this_term = this_term_data) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(matching_cols), names_to = "term", values_to = "effect")
  } else {
    # Normal case - dimensions match
    pred_long <- cbind(preds_df, by_var = obj$data[[by_var]], this_term = obj$data[[term_vars[1]]]) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(matching_cols), names_to = "term", values_to = "effect")
  }

  pred_long <- pred_long[mapply(grepl, pred_long$by_var, pred_long$term), ]
  pred_long <- aggregate(x = pred_long$effect, by = list(var1 = pred_long$by_var, var2 = pred_long$this_term), FUN = mean)

  # Apply the same logic for standard errors
  if (nrow(se_df) != nrow(obj$data)) {
    n_se <- nrow(se_df)
    by_var_data <- obj$data[[by_var]][seq_len(n_se)]
    this_term_data <- obj$data[[term_vars[1]]][seq_len(n_se)]

    se_long <- cbind(se_df, by_var = by_var_data, this_term = this_term_data) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(matching_cols), names_to = "term", values_to = "effect")
  } else {
    se_long <- cbind(se_df, by_var = obj$data[[by_var]], this_term = obj$data[[term_vars[1]]]) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(matching_cols), names_to = "term", values_to = "effect")
  }
  se_long <- se_long[mapply(grepl, se_long$by_var, se_long$term), ]
  se_long <- aggregate(x = se_long$effect, by = list(var1 = se_long$by_var, var2 = se_long$this_term), FUN = mean)
  pred_long$lower <- pred_long$x - 1.96 * se_long$x
  pred_long$upper <- pred_long$x + 1.96 * se_long$x

  if (islog) {
    pred_long$x <- exp(pred_long$x)
    pred_long$lower <- exp(pred_long$lower)
    pred_long$upper <- exp(pred_long$upper)
    ylim <- c(0, NA)
  } else {
    ylim <- c(NA, NA)
  }
  p_coef <- ggplot2::ggplot(pred_long, ggplot2::aes(x = var2, y = x, group = var1)) +
    ggplot2::geom_line(ggplot2::aes(colour = var1)) +
    ggplot2::geom_ribbon(ggplot2::aes(fill = var1, ymin = lower, ymax = upper), alpha = 0.2, show.legend = FALSE) +
    ggplot2::labs(y = "Partial effect", colour = by_var) +
    ggplot2::ylim(ylim)
  if (cdi) {
    p_coef <- p_coef + ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.position = "top",
      legend.key.size = ggplot2::unit(0.8, "lines"),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.8))
    )
  } else {
    p_coef <- p_coef +
      ggplot2::xlab(term_vars[1]) +
      ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "bottom")
  }
  return(p_coef)
}

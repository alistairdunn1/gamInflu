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
  pred_long <- cbind(preds_df, by_var = obj$data[[by_var]], this_term = obj$data[[term_vars[1]]]) %>% pivot_longer(cols = all_of(matching_cols), names_to = "term", values_to = "effect")
  pred_long <- pred_long[mapply(grepl, pred_long$by_var, pred_long$term), ]
  pred_long <- aggregate(x = pred_long$effect, by = list(var1 = pred_long$by_var, var2 = pred_long$this_term), FUN = mean)
  se_long <- cbind(se_df, by_var = obj$data[[by_var]], this_term = obj$data[[term_vars[1]]]) %>% pivot_longer(cols = all_of(matching_cols), names_to = "term", values_to = "effect")
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
  p_coef <- ggplot(pred_long, aes(x = var2, y = x, group = var1)) +
    geom_line(aes(colour = var1)) +
    geom_ribbon(aes(fill = var1, ymin = lower, ymax = upper), alpha = 0.2, show.legend = FALSE) +
    labs(y = "Partial effect", colour = by_var) +
    ylim(ylim)
  if (cdi) {
    p_coef <- p_coef + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), legend.title = element_blank(), legend.position = "top")
  } else {
    p_coef <- p_coef +
      xlab(term_vars[1]) +
      theme(legend.title = element_blank(), legend.position = "bottom")
  }
  return(p_coef)
}

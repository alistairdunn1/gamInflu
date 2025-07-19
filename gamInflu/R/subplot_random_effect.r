#' @title Subplot for Random Effects
#' @description Internal function to plot predicted effects for random effect terms (e.g., s(term, bs="re")).
#' @param obj A `gam_influence` object.
#' @param t The term name.
#' @param re_type Plot type: "points", "qq", "hist", or "caterpillar".
#' @param preds_df Data frame of predictions.
#' @param se_df Data frame of standard errors.
#' @param term_vars Variables in the term.
#' @param se_col Standard error column(s).
#' @param islog Logical; if TRUE, exponentiate effects.
#' @noRd
subplot_random_effect <- function(obj, t, term_vars, re_type = "qq", cdi = FALSE) {
  # Plot according to type
  re_type <- pmatch(re_type, c("points", "qq", "hist", "caterpillar"), nomatch = NA)
  if (is.na(re_type)) {
    message("Unsupported re_type for random effect plot. Use 'points', 'qq', 'hist', or 'caterpillar'")
    p_coef <- plot_spacer()
  } else if (re_type == 1) {
    p_coef <- subplot_random_effect_points(obj, term_vars)
  } else if (re_type == 2) {
    p_coef <- subplot_random_effect_qq(obj, term_vars)
  } else if (re_type == 3) {
    p_coef <- subplot_random_effect_histogram(obj, term_vars)
  } else if (re_type == 4) {
    p_coef <- subplot_random_effect_caterpillar(obj, term_vars)
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

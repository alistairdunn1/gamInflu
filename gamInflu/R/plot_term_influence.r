#' @title Term Influence Plot
#' @description Creates an influence plot showing the influence of each non-focus term on the focus term's index.
#' @param obj A `gam_influence` object containing calculated indices.
#' @return A ggplot object showing the influence of each non-focus term.
#' @export
#' @describeIn plot.gam_influence Creates an influence plot.
#' Shows the influence of each non-focus term on the focus term's index.
plot_term_influence <- function(obj) {
  df <- obj$calculated$influences
  if (is.null(df) || nrow(df) == 0) {
    return(ggplot() +
      labs(title = "Term Influence Plot", subtitle = "No non-focus terms to plot."))
  }

  # Convert level to numeric if possible
  if (is.factor(df$level) && all(!is.na(as.numeric(as.character(levels(df$level)))))) {
    df$level <- as.numeric(as.character(df$level))
  }

  # Get model term order from formula
  if (!is.null(obj$model$formula)) {
    term_order <- attr(terms(obj$model), "term.labels")
    # Only keep terms present in df
    term_order <- term_order[term_order %in% unique(df$term)]
    df$term <- factor(df$term, levels = term_order)
  } else {
    df$term <- factor(df$term, levels = unique(df$term))
  }

  ggplot(df, aes(x = level, y = influence, group = 1)) +
    geom_line(colour = "royalblue", alpha = 0.5) +
    geom_point(colour = "royalblue") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    facet_wrap(~term, ncol = 1, scales = "fixed") +
    ylim(0, NA) +
    labs(x = obj$focus, y = "Influence")
}

#' @title Term Influence Plot
#' @description Creates an influence plot showing the influence of each non-focus term on the focus term's index.
#' @param obj A `gam_influence` object containing calculated indices.
#' @return A ggplot object showing the influence of each non-focus term.
#' @importFrom ggplot2 ggplot aes geom_hline geom_line geom_point facet_wrap labs theme element_text ylim
#' @importFrom rlang .data
#' @export
#' @describeIn plot.gam_influence Creates an influence plot.
#' Shows the influence of each non-focus term on the focus term's index.
plot_term_influence <- function(obj) {
  df <- obj$calculated$influences
  if (is.null(df) || nrow(df) == 0) {
    return(ggplot() +
      labs(subtitle = "No non-focus terms to plot."))
  }

  # Convert level to numeric if possible
  if (is.factor(df$level) && all(!is.na(as.numeric(as.character(levels(df$level)))))) {
    df$level <- as.numeric(as.character(df$level))
  }

  # Get model term order from formula
  if (!is.null(obj$model$formula)) {
    term_order <- attr(terms(obj$model), "term.labels")
    # Robust matching to handle whitespace differences
    matched_terms <- sapply(term_order, function(formula_term) {
      available_terms <- unique(df$term)
      # First try exact match
      exact_match <- available_terms[available_terms == formula_term]
      if (length(exact_match) > 0) {
        return(exact_match[1])
      }
      # Try whitespace-stripped match
      formula_stripped <- strip_whitespace(formula_term)
      available_stripped <- strip_whitespace(available_terms)
      match_idx <- which(available_stripped == formula_stripped)
      if (length(match_idx) > 0) {
        return(available_terms[match_idx[1]])
      }
      return(NA)
    })
    # Only keep terms that actually matched
    matched_terms <- matched_terms[!is.na(matched_terms)]
    df$term <- factor(df$term, levels = matched_terms)
  } else {
    df$term <- factor(df$term, levels = unique(df$term))
  }

  ggplot(df, aes(x = .data$level, y = .data$influence, group = 1)) +
    geom_line(colour = "royalblue", alpha = 0.5) +
    geom_point(colour = "royalblue") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    facet_wrap(~term, ncol = 1, scales = "fixed") +
    ylim(0, NA) +
    labs(x = obj$focus, y = "Influence")
}

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

  # Get model term order from formula using the same method as calculate_influence
  if (!is.null(obj$model)) {
    # Use the same parse_formula_terms function that's used in calculate_influence
    parse_formula_terms <- function(frm) {
      rhs <- deparse(formula(frm)[[3]])
      rhs <- paste(rhs, collapse = " ")
      chars <- strsplit(rhs, "")[[1]]
      depth_paren <- 0
      buf <- ""
      tokens <- c()
      for (ch in chars) {
        if (ch == "(") depth_paren <- depth_paren + 1
        if (ch == ")") depth_paren <- max(0, depth_paren - 1)
        if (ch == "+" && depth_paren == 0) {
          token <- trimws(buf)
          if (nzchar(token)) tokens <- c(tokens, token)
          buf <- ""
        } else {
          buf <- paste0(buf, ch)
        }
      }
      last_token <- trimws(buf)
      if (nzchar(last_token)) tokens <- c(tokens, last_token)
      tokens <- tokens[nzchar(tokens)]
      tokens
    }

    term_order <- parse_formula_terms(obj$model)
    # Remove focus term from the order
    term_order <- setdiff(term_order, obj$focus)

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

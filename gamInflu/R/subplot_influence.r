#' @title Subplot for Influence Plot
#' @description Internal function to plot the influence for a term in CDI plots.
#' @param obj A `gam_influence` object.
#' @param term The term name.
#' @param focus_var The focus variable name.
#' @param cdi Logical indicating if the plot is for CDI (Cumulative Distribution Influence).
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline labs theme element_text facet_wrap ylim coord_flip
#' @importFrom rlang .data
#' @noRd

# Helper function for robust term matching
match_term_robust <- function(search_term, available_terms) {
  # First try exact match
  exact_match <- available_terms[available_terms == search_term]
  if (length(exact_match) > 0) {
    return(exact_match[1])
  }

  # If no exact match, try matching after stripping ALL whitespace
  search_stripped <- strip_whitespace(search_term)
  available_stripped <- strip_whitespace(available_terms)

  match_idx <- which(available_stripped == search_stripped)
  if (length(match_idx) > 0) {
    return(available_terms[match_idx[1]])
  }

  # No match found
  return(NULL)
}

subplot_influence <- function(obj, term, focus_var, cdi = FALSE) {
  message("Plotting influence for term: ", term)

  # --- Data Preparation ---
  focus_var <- obj$focus
  model_terms <- get_terms(obj, full = FALSE)

  # For CDI plots, we expect a single term
  if (cdi) {
    # Identify variables associated with the term (can be 1 or 2 for smooths)
    term_vars <- all.vars(rlang::parse_expr(term))
    term_vars_in_model <- term_vars[term_vars %in% model_terms]
    if (length(term_vars_in_model) == 0) {
      stop(paste0(
        "None of the variables in the supplied term ('", term, "') match the model term.\n",
        "Valid terms are: ", paste(setdiff(model_terms, obj$focus), collapse = ", ")
      ))
    }
    if (any(term_vars_in_model == focus_var)) {
      stop("CDI plot cannot be generated for the focus term itself.")
    }

    # Fix: Robust matching against the full term name, handling whitespace differences
    available_terms <- unique(obj$calculated$influences$term)
    matched_term <- match_term_robust(term, available_terms)

    if (is.null(matched_term)) {
      stop(
        "No influence data found for term: ", term, ". Available terms: ",
        paste(available_terms, collapse = ", ")
      )
    }

    influ_data <- obj$calculated$influences[obj$calculated$influences$term == matched_term, ]
  } else {
    # For multi-term influence plots, use all available influence data
    influ_data <- obj$calculated$influences
    if (is.null(influ_data) || nrow(influ_data) == 0) {
      return(ggplot2::ggplot() +
        ggplot2::labs(subtitle = "No influence data available."))
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
        available_terms <- unique(influ_data$term)
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
      influ_data$term <- factor(influ_data$term, levels = matched_terms)
    } else {
      influ_data$term <- factor(influ_data$term, levels = unique(influ_data$term))
    }
  }

  # Apply transformations
  influ_data$influence <- exp(influ_data$influence)
  influ_data$influence <- influ_data$influence / mean(influ_data$influence)

  # Safety check: ensure we have valid influence data after transformations
  if (all(is.na(influ_data$influence))) {
    warning("All influence values are NA for term: ", term)
    influ_data$influence <- rep(1, nrow(influ_data)) # Default to no influence
  }

  # Ensure proper ordering of levels for plotting
  # This prevents issues where "1992" appears before "2010" etc.
  if (is.factor(influ_data$level)) {
    numeric_levels <- suppressWarnings(as.numeric(as.character(levels(influ_data$level))))
    if (!any(is.na(numeric_levels))) {
      # If all levels can be converted to numeric, reorder the factor
      ordered_levels <- levels(influ_data$level)[order(numeric_levels)]
      influ_data$level <- factor(influ_data$level, levels = ordered_levels)
    }
  } else if (is.character(influ_data$level)) {
    # Try to convert character to numeric for ordering
    numeric_values <- suppressWarnings(as.numeric(influ_data$level))
    if (!any(is.na(numeric_values))) {
      influ_data$level <- factor(influ_data$level, levels = unique(influ_data$level)[order(numeric_values)])
    }
  }

  ylim <- c(
    pmin(0.75, min(influ_data$influence, na.rm = TRUE)),
    pmax(1.25, max(influ_data$influence, na.rm = TRUE))
  )

  if (cdi || length(unique(influ_data$term)) == 1) {
    # Single term plot (CDI or single term)
    ggplot2::ggplot(influ_data, ggplot2::aes(x = .data$level, y = .data$influence, group = .data$term)) +
      ggplot2::geom_line(colour = "royalblue", alpha = 0.5) +
      ggplot2::geom_point(colour = "royalblue") +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
      ggplot2::labs(y = "Influence") +
      ggplot2::theme(
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank()
      ) +
      ggplot2::ylim(ylim) +
      ggplot2::coord_flip()
  } else {
    # Multi-term plot with faceting (like plot_term_influence)
    ggplot2::ggplot(influ_data, ggplot2::aes(x = .data$level, y = .data$influence, group = 1)) +
      ggplot2::geom_line(colour = "royalblue", alpha = 0.5) +
      ggplot2::geom_point(colour = "royalblue") +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
      ggplot2::facet_wrap(~ .data$term, ncol = 1, scales = "fixed") +
      ggplot2::ylim(0, NA) +
      ggplot2::labs(x = obj$focus, y = "Influence") +
      ggplot2::theme(strip.text = ggplot2::element_text(hjust = 0))
  }
}

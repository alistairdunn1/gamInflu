#' @title Subplot for Influence Plot
#' @description Internal function to plot the influence for a term in CDI plots.
#' @param obj A `gam_influence` object.
#' @param term The term name.
#' @param focus_var The focus variable name.
#' @param cdi Logical indicating if the plot is for CDI (Cumulative Distribution Influence).
#' @return A ggplot object.
#' @noRd

# Helper function for robust term matching
match_term_robust <- function(search_term, available_terms) {
  # First try exact match
  exact_match <- available_terms[available_terms == search_term]
  if (length(exact_match) > 0) {
    return(exact_match[1])
  }

  # If no exact match, try matching after stripping ALL whitespace
  strip_whitespace <- function(x) gsub("\\s+", "", x)
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

  influ_data <- obj$calculated$influences[obj$calculated$influences$term == matched_term, ] # Apply transformations
  influ_data$influence <- exp(influ_data$influence)
  influ_data$influence <- influ_data$influence / mean(influ_data$influence)

  # Safety check: ensure we have valid influence data after transformations
  if (all(is.na(influ_data$influence))) {
    warning("All influence values are NA for term: ", term)
    influ_data$influence <- rep(1, nrow(influ_data)) # Default to no influence
  }

  ylim <- c(
    pmin(0.75, min(influ_data$influence, na.rm = TRUE)),
    pmax(1.25, max(influ_data$influence, na.rm = TRUE))
  )

  if (length(unique(influ_data$term)) == 1) {
    ggplot(influ_data, aes(x = level, y = influence, group = term)) +
      geom_line(colour = "royalblue", alpha = 0.5) +
      geom_point(colour = "royalblue") +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(y = "Influence") +
      ggplot2::theme(
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank()
      ) +
      ylim(ylim) +
      ggplot2::coord_flip()
  } else {
    ggplot(influ_data, aes(x = level, y = influence, group = term, colour = term)) +
      geom_line(alpha = 0.5) +
      geom_point() +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(y = "Influence", colour = "Term") +
      ggplot2::theme(
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank()
      ) +
      ylim(ylim) +
      ggplot2::coord_flip() +
      ggplot2::theme(
        legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
        legend.key.size = grid::unit(0.8, "lines"),
        legend.text = ggplot2::element_text(size = ggplot2::rel(0.8))
      )
  }
}

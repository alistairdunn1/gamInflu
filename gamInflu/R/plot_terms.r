#' @title Plot Predicted Effects for Model Terms
#' @description Plots the predicted effects for each model term, with options for single or all terms,
#' including by-variable panels and random effects.
#' @param obj A `gam_influence` object containing calculated predictions and standard errors.
#' @param term The character name of the model term to plot (e.g., `"s(temp)"`). Alternatively, it can
#' be a numeric index of the term in the model. If a numeric index is provided, it will be converted
#' to the corresponding term name. If `NULL`, all terms will be plotted.
#' @param re_type Character; for random effects, one of "points", "qq", "hist", or "caterpillar".
#' @param cdi Logical indicating if the plot is for CDI (Cumulative Distribution Influence).
#' @return A ggplot object (or patchwork if multiple terms).
#' @export
plot_terms <- function(obj, term = NULL, re_type = "points", cdi = FALSE) {
  # --- Setup ---
  if (is.numeric(term) && length(term) == 1 && term == as.integer(term)) {
    all_terms <- get_terms(obj, full = TRUE)
    if (term > length(all_terms) || term < 1) {
      stop("Term index out of bounds. There are ", length(all_terms), " valid terms:\n  ", paste(all_terms, collapse = "\n  "))
    }
    term <- all_terms[term]
  }

  preds_df <- obj$calculated$predictions
  se_df <- obj$calculated$prediction_se
  terms_vec <- get_terms(obj, full = TRUE)
  term_list <- if (is.null(term)) terms_vec else term
  islog <- isTRUE(obj$islog)

  plots <- lapply(term_list, function(t) {
    term_vars <- all.vars(rlang::parse_expr(t))
    is_random <- grepl('bs\\s*=\\s*"re"', t)
    is_by <- grepl("by\\s*=", t)
    is_factor <- is.factor(obj$data[[term_vars[1]]])
    is_tensor2d <- length(term_vars) == 2 && !is_by && !is_random && !is_factor

    if (t == obj$focus) {
      subplot_focus_effect(obj, t, term_vars, cdi = cdi)
    } else if (is_tensor2d) {
      subplot_tensor2d_effect(obj, t, term_vars, cdi = cdi)
    } else if (length(term_vars) == 1 && !is_by && !is_random && !is_factor) {
      subplot_continuous_effect(obj, t, term_vars, cdi = cdi)
    } else if (length(term_vars) == 2 && is_by && !is_random && !is_factor) {
      subplot_by_variable(obj, t, term_vars, cdi = cdi)
    } else if (is_random && !is_by) {
      subplot_random_effect(obj, t, term_vars, re_type = re_type, cdi = cdi)
    } else if (length(term_vars) == 1 && !is_by && !is_random && is_factor) {
      subplot_factor_effect(obj, t, term_vars, cdi = cdi)
    } else {
      message(
        "Plotting for terms of type '", t, "' with variables ",
        paste(term_vars, collapse = ", "),
        " is not yet implemented yet.\n   [is_random = ", is_random, ", is_by = ", is_by, ", is_factor = ", is_factor, ", is_tensor2d = ", is_tensor2d, "]"
      )
      patchwork::plot_spacer()
    }
  })

  if (length(plots) == 1) {
    return(plots[[1]])
  } else {
    patchwork::wrap_plots(plots)
  }
}

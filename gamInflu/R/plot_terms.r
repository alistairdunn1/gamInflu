#' @title Plot Predicted Effects for Model Terms
#' @description Plots the predicted effects for each model term, with options for single or all terms,
#' including by-variable panels and random effects.
#' @param obj A `gam_influence` object containing calculated predictions and standard errors.
#' @param term The character name of the model term to plot (e.g., `"s(temp)"`). Alternatively, it can
#' be a numeric index of the term in the model. If a numeric index is provided, it will be converted
#' to the corresponding term name. If `NULL`, all terms will be plotted.
#' @param type Character; for random effects, one of "point", "bar", or "violin".
#' @return A ggplot object (or patchwork if multiple terms).
#' @export
plot_terms <- function(obj, term = NULL, type = "point") {
  # --- Setup ---
  if (is.numeric(term) && length(term) == 1 && term == as.integer(term)) {
    all_terms <- get_terms(obj, full = TRUE)
    if (term > length(all_terms) || term < 1) {
      stop("Term index out of bounds. There are ", length(all_terms), " valid terms:\n  ", paste(all_terms, collapse = "\n  "), call. = FALSE)
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

    if (t == obj$focus) {
      message("Plotting focus term: ", t, ".")
      subplot_focus_effect(obj, t, term_vars, cdi = FALSE)
    } else if (length(term_vars) == 1 && !is_by && !is_random && !is_factor) {
      message("Plotting term: ", t, " as a continuous effect.")
      subplot_continuous_effect(obj, t, term_vars, cdi = FALSE)
    } else if (length(term_vars) == 2 && is_by && !is_random && !is_factor) {
      message("Plotting term: ", t, " as a continuous by-variable effect.")
      subplot_by_variable(obj, t, term_vars, cdi = FALSE)
    } else if (is_random && !is_by) {
      message("Plotting term: ", t, " as a random effect.")
      subplot_random_effect(obj, t, term_vars, type = type, cdi = FALSE)
    } else if (length(term_vars) == 1 && !is_by && !is_random && is_factor) {
      message("Plotting term: ", t, " as a factor effect.")
      subplot_factor_effect(obj, t, term_vars, cdi = FALSE)
    } else {
      message("Plotting term:", t, " as an unsupported effect.")
    }
  })

  if (length(plots) == 1) {
    return(plots[[1]])
  } else {
    patchwork::wrap_plots(plots)
  }
}

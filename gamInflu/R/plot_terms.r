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
  islog <- obj$islog

  plots <- lapply(term_list, function(t) {
    term_vars <- all.vars(rlang::parse_expr(t))
    is_random <- grepl('bs\\s*=\\s*"re"', t)
    is_by <- grepl("by\\s*=", t)
    se_col <- if (!is.null(se_df)) grep(paste0("(^|\\(|\\:)", t, "($|\\)|\\:|\\d+)"), colnames(se_df), value = TRUE) else character(0)

    if (is_random) {
      # Random effect plot
      if (!(t %in% colnames(preds_df)) || !(term_vars[1] %in% colnames(preds_df))) {
        warning(paste("Random effect term", t, "not found in predictions."))
        return(ggplot())
      }
      re_levels <- preds_df[[term_vars[1]]]
      re_effect <- preds_df[[t]]
      re_se <- if (length(se_col) > 0 && se_col[1] %in% colnames(se_df)) se_df[[se_col[1]]] else NA
      if (islog) {
        re_effect <- exp(re_effect)
        re_se <- if (!is.na(re_se)) exp(re_se) else re_se
      }
      df <- data.frame(level = re_levels, effect = re_effect, se = re_se)
      if (type == "violin") {
        p <- ggplot(df, aes(x = level, y = effect)) +
          geom_violin(fill = "grey80") +
          labs(title = t, x = term_vars[1], y = "Random Effect")
      } else if (type == "bar") {
        p <- ggplot(df, aes(x = level, y = effect)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          geom_errorbar(aes(ymin = effect - 1.96 * se, ymax = effect + 1.96 * se), width = 0.2, na.rm = TRUE) +
          labs(title = t, x = term_vars[1], y = "Random Effect")
      } else {
        p <- ggplot(df, aes(x = level, y = effect)) +
          geom_point(size = 3, colour = "steelblue") +
          geom_errorbar(aes(ymin = effect - 1.96 * se, ymax = effect + 1.96 * se), width = 0.2, na.rm = TRUE) +
          labs(title = t, x = term_vars[1], y = "Random Effect")
      }
    } else if (is_by) {
      # By-variable panel
      by_var <- sub(".*by\\s*=\\s*([^,\\)]+).*", "\\1", t)
      main_var <- term_vars[1]
      df <- preds_df
      effect <- df[[t]]
      ymin <- effect - 1.96 * if (length(se_col) > 0) se_df[[se_col[1]]] else 0
      ymax <- effect + 1.96 * if (length(se_col) > 0) se_df[[se_col[1]]] else 0
      if (islog) {
        effect <- exp(effect)
        ymin <- exp(ymin)
        ymax <- exp(ymax)
      }
      df$effect <- effect
      df$ymin <- ymin
      df$ymax <- ymax
      p <- ggplot(df, aes(x = !!rlang::sym(main_var), y = effect, colour = !!rlang::sym(by_var))) +
        geom_line(aes(group = !!rlang::sym(by_var))) +
        geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = !!rlang::sym(by_var)), alpha = 0.2, colour = NA) +
        labs(title = t, x = main_var, y = "Effect", colour = by_var, fill = by_var) +
        facet_wrap(vars(!!rlang::sym(by_var)))
    } else if (is.factor(obj$data[[term_vars[1]]])) {
      # Factor variable
      effect <- preds_df[[t]]
      se <- if (length(se_col) > 0) se_df[[se_col[1]]] else NA
      if (islog) {
        effect <- exp(effect)
        se <- ifelse(!is.na(se), exp(se), se)
      }
      df <- data.frame(
        level = obj$data[[term_vars[1]]], effect = effect, se = se
      )
      p <- ggplot(df, aes(x = level, y = effect)) +
        geom_point(size = 3, colour = "royalblue") +
        geom_errorbar(aes(ymin = effect - 1.96 * se, ymax = effect + 1.96 * se), colour = "royalblue", alpha = 0.5, width = 0.2, na.rm = TRUE) +
        labs(title = t, x = term_vars[1], y = "Effect")
    } else {
      # Continuous variable
      effect <- preds_df[[t]]
      se <- if (length(se_col) > 0) se_df[[se_col[1]]] else NA
      ymin <- effect - 1.96 * se
      ymax <- effect + 1.96 * se
      if (islog) {
        effect <- exp(effect)
        ymin <- exp(ymin)
        ymax <- exp(ymax)
      }
      df <- data.frame(
        x = preds_df[[term_vars[1]]], effect = effect, ymin = ymin, ymax = ymax
      )
      p <- ggplot(df, aes(x = x, y = effect)) +
        geom_line(colour = "royalblue") +
        geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, fill = "royalblue") +
        labs(title = t, x = term_vars[1], y = "Effect")
    }
    p
  })

  if (length(plots) == 1) {
    return(plots[[1]])
  } else {
    patchwork::wrap_plots(plots)
  }
}

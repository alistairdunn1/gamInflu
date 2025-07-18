#' @describeIn plot.gam_influence Creates a Coefficient-Distribution-Influence (CDI) plot.
#' This is a multi-panel plot that visualizes a term's effect, its data
#' distribution, and its influence on the focus term.
#'
#' @param obj A `gam_influence` object.
#' @param term The character name of the model term to plot (e.g., `"s(temp)"`). Alternatively, it can be a numeric index of the term in the model. If a numeric index is provided, it will be converted to the corresponding term name.
#' @return A patchwork ggplot object.
#' @export
#'
plot_cdi <- function(obj, term) {
  # --- Setup ---
  if (is.numeric(term) && length(term) == 1 && term == as.integer(term)) {
    all_terms <- get_terms(obj, full = TRUE)
    if (term > length(all_terms) || term < 1) {
      stop("Term index out of bounds. There are ", length(all_terms), " valid terms:\n  ", paste(all_terms, collapse = "\n  "), call. = FALSE)
    }
    term <- all_terms[term]
  }

  # --- Data Preparation ---
  preds_df <- obj$calculated$predictions
  se_df <- obj$calculated$prediction_se
  obj_data <- obj$data
  focus_var <- obj$focus
  model_terms <- get_terms(obj, full = FALSE)

  # Identify variables associated with the term (can be 1 or 2 for smooths)
  # This handles terms like `s(temp)`, `te(lon,lat)`, and `s(day, by=vessel)`
  term_vars <- all.vars(rlang::parse_expr(term))

  # Extract the variable name(s) from the supplied term (e.g., s(temp) -> temp)
  # and check against obj$terms (which are variable names)
  # Allow for multiple variables in smooths (e.g., te(lon,lat))
  term_vars_in_model <- term_vars[term_vars %in% model_terms]
  if (length(term_vars_in_model) == 0) {
    stop(
      paste0(
        "None of the variables in the supplied term ('", term, "') match the model term.\n",
        "Valid terms are: ", paste(setdiff(model_terms, obj$focus), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (any(term_vars_in_model == focus_var)) {
    stop("CDI plot cannot be generated for the focus term itself.", call. = FALSE)
  }

  # --- Plot 1: Coefficient Plot (The term's effect) ---
  call_obj <- rlang::parse_expr(term)
  call_list <- as.list(call_obj)
  is_by_factor <- "by" %in% names(call_list) && !is.null(call_list$by)
  is_random_effect <- grepl('bs\\s*=\\s*"re"', term)

  if (length(term_vars) == 1 && !is_by_factor && !is_random_effect) {
    var1 <- term_vars[1]
    # Try to get standard errors for this term
    se_df <- obj$calculated$prediction_se
    se_col <- if (!is.null(se_df)) grep(paste0("(^|\\(|\\:)", term, "($|\\)|\\:|\\d+)"), colnames(se_df), value = TRUE) else character(0)
    if (is.factor(obj_data[[var1]])) {
      # If the variable is a factor, aggregate by its levels
      p_coef_data <- aggregate(
        cbind(effect = preds_df[[term]], se = if (length(se_col) > 0) se_df[[se_col[1]]] else NA),
        by = list(var_level = as.character(obj_data[[var1]])),
        FUN = mean
      )
      p_coef_data$lower <- p_coef_data$effect - 1.96 * p_coef_data$se
      p_coef_data$upper <- p_coef_data$effect + 1.96 * p_coef_data$se
      if (obj$islog) {
        p_coef_data$lower <- exp(p_coef_data$lower)
        p_coef_data$upper <- exp(p_coef_data$upper)
        p_coef_data$effect <- exp(p_coef_data$effect)
      }
      p_coef <- ggplot(p_coef_data, aes(x = var_level, y = effect)) +
        geom_point() +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, na.rm = TRUE) +
        labs(y = "Coefficient") +
        theme(
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank()
        )
    } else {
      # Otherwise, treat as numeric
      p_coef_data <- aggregate(
        cbind(effect = preds_df[[term]], se = if (length(se_col) > 0) se_df[[se_col[1]]] else NA),
        by = list(var_level = obj_data[[var1]]),
        FUN = mean
      )
      p_coef_data$lower <- p_coef_data$effect - 1.96 * p_coef_data$se
      p_coef_data$upper <- p_coef_data$effect + 1.96 * p_coef_data$se
      p_coef <- ggplot(p_coef_data, aes(x = var_level, y = effect)) +
        geom_line() +
        geom_point() +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, na.rm = TRUE) +
        labs(y = "Coefficient") +
        theme(
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank()
        )
    }
  } else if (length(term_vars) >= 2 && !is_by_factor && !is_random_effect) {
    # Case 2: 2D smooth term, e.g., te(lon, lat)
    var1 <- term_vars[1]
    var2 <- term_vars[2]
    p_coef <- ggplot(preds_df, aes(x = !!rlang::sym(var1), y = !!rlang::sym(var2), fill = !!rlang::sym(term))) +
      geom_raster() +
      scale_fill_viridis_c() +
      labs(x = NULL, y = var2, fill = "Effect")
  } else if (length(term_vars) == 2 && is_by_factor && !is_random_effect) {
    # Case 3: By-variable panel, e.g., s(temp, by=vessel)
    factor_var <- call_list$by
    matching_cols <- names(preds_df)[
      vapply(names(preds_df), function(nm) {
        grepl(term_vars[1], nm) && grepl(factor_var, nm)
      }, logical(1))
    ]
    pred_long <- cbind(preds_df, factor_var = obj_data[[factor_var]], this_term = obj_data[[term_vars[1]]]) %>% pivot_longer(cols = all_of(matching_cols), names_to = "term", values_to = "effect")
    pred_long <- pred_long[mapply(grepl, pred_long$factor_var, pred_long$term), ]
    pred_long <- aggregate(x = pred_long$effect, by = list(var1 = pred_long$factor_var, var2 = pred_long$this_term), FUN = mean)
    se_long <- cbind(se_df, factor_var = obj_data[[factor_var]], this_term = obj_data[[term_vars[1]]]) %>% pivot_longer(cols = all_of(matching_cols), names_to = "term", values_to = "effect")
    se_long <- se_long[mapply(grepl, se_long$factor_var, se_long$term), ]
    se_long <- aggregate(x = se_long$effect, by = list(var1 = se_long$factor_var, var2 = se_long$this_term), FUN = mean)
    pred_long$lower <- pred_long$x - 1.96 * se_long$x
    pred_long$upper <- pred_long$x + 1.96 * se_long$x

    p_coef <- ggplot(pred_long, aes(x = var2, y = x, group = var1)) +
      geom_line(aes(colour = var1)) +
      geom_ribbon(aes(fill = var1, ymin = lower, ymax = upper), alpha = 0.2, show.legend = FALSE) +
      labs(y = "Effect", colour = factor_var) +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top"
      )
  } else if (is_random_effect) {
    # Case 4: Random effect
    message("CDI plot not fully implemented for random effects.")
    p_coef <- plot_spacer()
  } else {
    # Case 5: Other term types (e.g., by variable, interaction, etc.)
    message("CDI plot not fully implemented for this term type. Please check the term format.")
    p_coef <- plot_spacer()
  }
  # --- Plot 2: Distribution Plot ---
  p_dist_data <- obj_data %>%
    dplyr::group_by(.data[[focus_var]], .data[[term_vars[1]]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      n = ((n * 8) / max(n, na.rm = TRUE)),
      n = if (length(unique(n)) == 1) 1 else n
    )

  p_dist <- ggplot(p_dist_data, aes(x = !!rlang::sym(focus_var), y = !!rlang::sym(term_vars[1]))) +
    geom_point(colour = "royalblue", size = p_dist_data$n, alpha = 0.5) +
    coord_flip() +
    labs(y = term_vars[1], x = focus_var)

  # --- Plot 3: Influence Plot ---
  # Find the influence row(s) for any of the variables in the term
  influ_data <- subset(obj$calculated$influences, term %in% term_vars_in_model) %>%
    dplyr::mutate(influence = exp(influence)) %>%
    dplyr::mutate(influence = influence / mean(influence))
  p_influ <- ggplot(influ_data, aes(x = level, y = influence, group = 1)) +
    geom_line(colour = "royalblue", alpha = 0.5) +
    geom_point(colour = "royalblue") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(y = "Influence") +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank()
    ) +
    coord_flip()

  # --- Combine with Patchwork: p_coef on top left, p_dist bottom left, p_influ right ---
  layout <- "
  AB
  CD
  "

  patchwork::wrap_plots(
    A = p_coef,
    B = plot_spacer(),
    C = p_dist,
    D = p_influ,
    design = layout,
    widths = c(1, 0.5),
    heights = c(1, 2)
  )
}

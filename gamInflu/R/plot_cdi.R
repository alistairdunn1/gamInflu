#' @describeIn plot.gam_influence Creates a Coefficient-Distribution-Influence (CDI) plot.
#' This is a multi-panel plot that visualizes a term's effect, its data
#' distribution, and its influence on the focus term.
#'
#' @param obj A `gam_influence` object.
#' @param term The character name of the model term to plot (e.g., `"s(temp)"`).
#' @return A patchwork ggplot object.
#' @export
#'
plot_cdi <- function(obj, term) {
  # --- Data Preparation ---
  preds_df <- obj$calculated$predictions
  obj_data <- obj$data
  focus_var <- obj$focus

  # Identify variables associated with the term (can be 1 or 2 for smooths)
  # This handles terms like `s(temp)`, `te(lon,lat)`, and `s(day, by=vessel)`
  term_vars <- all.vars(rlang::parse_expr(term))

  # Extract the variable name(s) from the supplied term (e.g., s(temp) -> temp)
  # and check against obj$terms (which are variable names)
  # Allow for multiple variables in smooths (e.g., te(lon,lat))
  term_vars_in_model <- term_vars[term_vars %in% obj$terms]
  if (length(term_vars_in_model) == 0) {
    stop(
      paste0(
        "None of the variables in the supplied term ('", term, "') match the model terms: ",
        paste(obj$terms, collapse = ", ")
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

  if (length(term_vars) == 1 && !is_by_factor) {
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
        by = list(var_level = preds_df[[var1]]),
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
  } else if (length(term_vars) >= 2 && !is_by_factor) {
    # Case 2: 2D smooth term, e.g., te(lon, lat)
    var1 <- term_vars[1]
    var2 <- term_vars[2]
    p_coef <- ggplot(preds_df, aes_string(x = var1, y = var2, fill = term)) +
      geom_raster() +
      scale_fill_viridis_c() +
      labs(x = NULL, y = var2, fill = "Effect")
  } else {
    # Case 3: Nested factor (by-variable), e.g., s(day, by = vessel)
    by_value <- as.list(rlang::parse_expr(term))$by
    factor_var <- as.character(by_value)
    numeric_var <- term_vars[!term_vars %in% all.vars(by_value)]

    p_coef <- ggplot(preds_df, aes_string(x = numeric_var, y = paste0("exp(", term, ")"), colour = factor_var)) +
      geom_line(aes_string(group = factor_var)) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(y = "Effect", colour = factor_var) +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
      ) +
      facet_wrap(vars(!!rlang::sym(factor_var)))
  }
  # --- Plot 2: Distribution Plot ---
  p_dist_data <- obj_data %>%
    dplyr::group_by(.data[[focus_var]], .data[[term_vars[1]]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(n = n / max(n, na.rm = TRUE))

  p_dist <- ggplot(p_dist_data, aes_string(x = focus_var, y = term_vars[1])) +
    geom_point(fill = "blue", size = p_dist_data$n * 10, alpha = 0.5) +
    coord_flip() +
    labs(y = term_vars[1], x = focus_var)

  # --- Plot 3: Influence Plot ---
  # Find the influence row(s) for any of the variables in the term
  influ_data <- subset(obj$calculated$influences, term %in% term_vars_in_model) %>%
    mutate(influence = exp(influence)) %>%
    mutate(influence = influence / mean(influence))
  p_influ <- ggplot(influ_data, aes(x = level, y = influence, group = 1)) +
    geom_line() +
    geom_point() +
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

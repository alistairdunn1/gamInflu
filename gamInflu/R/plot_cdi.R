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
  term_vars <- all.vars(rlang::parse_expr(term))
  term_vars_in_model <- term_vars[term_vars %in% model_terms]
  if (length(term_vars_in_model) == 0) {
    stop(paste0(
      "None of the variables in the supplied term ('", term, "') match the model term.\n",
      "Valid terms are: ", paste(setdiff(model_terms, obj$focus), collapse = ", ")
    ), call. = FALSE)
  }
  if (any(term_vars_in_model == focus_var)) {
    stop("CDI plot cannot be generated for the focus term itself.", call. = FALSE)
  }

  # --- Plot 1: Coefficient Plot (The term's effect) ---
  call_obj <- rlang::parse_expr(term)
  call_list <- as.list(call_obj)
  is_by_factor <- "by" %in% names(call_list) && !is.null(call_list$by)
  is_random_effect <- grepl('bs\\s*=\\s*"re"', term)

  if (is_random_effect) {
    p_coef <- subplot_random_effect(obj, term, "point", term_vars, cdi = TRUE)
  } else if (is_by_factor) {
    p_coef <- subplot_by_variable(obj, term, term_vars, cdi = TRUE)
  } else if (is.factor(obj_data[[term_vars[1]]])) {
    p_coef <- subplot_factor_effect(obj, term, term_vars, cdi = TRUE)
  } else if (length(term_vars) == 1 && is.numeric(obj_data[[term_vars[1]]])) {
    p_coef <- subplot_continuous_effect(obj, term, term_vars, cdi = TRUE)
  } else {
    message("Unsupported term type for CDI plot. Only continuous, factor, or random effects are supported.", call. = FALSE)
    p_coef <- plot_spacer()
  }

  # --- Plot 2: Distribution Plot ---
  p_dist <- subplot_distribution(obj_data, term, focus_var)

  # --- Plot 3: Influence Plot ---
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

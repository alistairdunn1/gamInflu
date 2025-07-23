#' @describeIn plot.gam_influence Creates a Coefficient-Distribution-Influence (CDI) plot.
#' This is a multi-panel plot that visualizes a term's effect, its data
#' distribution, and its influence on the focus term.
#'
#' @param obj A `gam_influence` object.
#' @param term The character name of the model term to plot (e.g., `"s(temp)"`). Alternatively, it can be a numeric index of the term in the model. If a numeric index is provided, it will be converted to the corresponding term name.
#' @param re_type Character; for random effects, one of "points", "qq", "hist", or "caterpillar".
#' @return A patchwork ggplot object.
#' @export
#'
plot_cdi <- function(obj, term, re_type = "points") {
  # --- Setup ---
  if (is.numeric(term) && length(term) == 1 && term == as.integer(term)) {
    all_terms <- get_terms(obj, full = TRUE)
    if (term > length(all_terms) || term < 1) {
      stop("Term index out of bounds. There are ", length(all_terms), " valid terms:\n  ", paste(all_terms, collapse = "\n  "))
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
    ))
  }
  if (any(term_vars_in_model == focus_var)) {
    stop("CDI plot cannot be generated for the focus term itself.")
  }

  # --- Plot 1: Coefficient Plot (The term's effect) ---
  p_coef <- plot_terms(obj, term = term, re_type = re_type, cdi = TRUE)

  # --- Plot 2: Distribution Plot ---
  p_dist <- subplot_distribution(obj = obj, term = term, focus_var = focus_var)

  # --- Plot 3: Influence Plot ---
  p_influ <- subplot_influence(obj = obj, term = term, focus_var = focus_var, cdi = TRUE)

  # --- Combine with Patchwork: p_coef on top left, p_dist bottom left, p_influ right ---
  layout <- "
  AB
  CD
  "

  patchwork::wrap_plots(
    A = p_coef,
    B = patchwork::plot_spacer(),
    C = p_dist,
    D = p_influ,
    design = layout,
    widths = c(1, 0.5),
    heights = c(1, 2)
  )
}

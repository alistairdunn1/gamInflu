#' @title Coefficient-Distribution-Influence (CDI) Plot
#' @description Creates a Coefficient-Distribution-Influence (CDI) plot.
#' This is a multi-panel plot that visualises a term's effect, its data
#' distribution, and its influence on the focus term.
#' @param obj A `gam_influence` object.
#' @param term The character name of the model term to plot (e.g., `"s(temp)"`). Alternatively, it can be a numeric index of the term in the model. If a numeric index is provided, it will be converted to the corresponding term name.
#' @param re_type Character; for random effects, one of "points", "qq", "hist", or "caterpillar".
#' @param show_summary Logical; whether to show the model summary table in the top right corner (default: FALSE).
#' @return A patchwork ggplot object.
#' @importFrom rlang parse_expr
#' @importFrom patchwork wrap_plots plot_spacer
#' @export
plot_cdi <- function(obj, term, re_type = "points", show_summary = FALSE) {
  # --- Setup ---
  if (is.numeric(term) && length(term) == 1 && term == as.integer(term)) {
    all_terms <- get_terms(obj, full = TRUE)
    if (term > length(all_terms) || term < 1) {
      stop("Term index out of bounds. There are ", length(all_terms), " valid terms:\n  ", paste(all_terms, collapse = "\n  "), call. = FALSE)
    }
    term <- all_terms[term]
  }

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
    ), call. = FALSE)
  }
  if (any(term_vars_in_model == focus_var)) {
    stop("CDI plot cannot be generated for the focus term itself.", call. = FALSE)
  }

  # --- Plot 1: Coefficient Plot (The term's effect) ---
  p_coef <- plot_terms(obj, term = term, re_type = re_type, cdi = TRUE)

  # --- Plot 2: Distribution Plot ---
  p_dist <- plot_term_distribution(obj = obj, term = term, by = FALSE)

  # --- Plot 3: Influence Plot ---
  p_influ <- subplot_influence(obj = obj, term = term, focus_var = focus_var, cdi = TRUE)

  # --- Combine plots based on show_summary parameter ---
  if (show_summary) {
    # --- Plot 4: Model Summary Table ---
    p_summary <- create_model_summary_table(obj)
  } else {
    p_summary <- patchwork::plot_spacer()
  }
  # --- Combine with Patchwork: p_coef on top left, p_dist bottom left, p_influ right ---
  layout <- "
    AB
    CD
    "

  patchwork::wrap_plots(
    A = p_coef,
    B = p_summary,
    C = p_dist,
    D = p_influ,
    design = layout,
    widths = c(1, 0.5),
    heights = c(1, 2)
  )
}


#' @title Create Model Summary Table for CDI Plot
#' @description Creates a table showing model deviance and R² information
#' @param obj A gam_influence object
#' @return A ggplot object containing the summary table
#' @importFrom ggplot2 ggplot aes geom_text theme_void theme element_blank element_rect
create_model_summary_table <- function(obj) {
  # Extract model information
  model <- obj$model

  # Get deviance
  model_deviance <- round(deviance(model), 2)
  null_deviance <- round(model$null.deviance, 2)

  # Get R² (adjusted if available)
  model_summary <- summary(model)
  r_squared <- if (!is.null(model_summary$r.sq)) {
    round(model_summary$r.sq, 3)
  } else {
    NA
  }

  adj_r_squared <- if (!is.null(model_summary$dev.expl)) {
    round(model_summary$dev.expl, 3)
  } else {
    NA
  }

  # Create data frame for the table
  table_data <- data.frame(
    Label = c("Deviance:", "Null Dev:", "Dev Expl:", "R2:"),
    Value = c(
      as.character(model_deviance),
      as.character(null_deviance),
      if (!is.na(adj_r_squared)) paste0(as.character(adj_r_squared * 100), "%") else "N/A",
      if (!is.na(r_squared)) as.character(r_squared) else "N/A"
    ),
    y = c(5, 4, 3, 2) # Position from top to bottom
  )

  # Create ggplot table
  ggplot2::ggplot(table_data, ggplot2::aes(x = 0.1, y = .data$y, label = paste(.data$Label, .data$Value))) +
    ggplot2::geom_text(hjust = 0, size = 2.5) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0.5, 6.5)
}

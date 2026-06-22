#' @title Residual-Implied Coefficients Plot
#' @description Creates a multi-panel plot comparing, for each level of a modelled
#'   factor term, the term's *actual* fitted coefficient against the *implied*
#'   coefficient derived from the residuals at each level of the focus term. The
#'   implied coefficient is the actual term coefficient plus the mean standardised
#'   residual for each `term-level x focus-level` cell. Where the implied
#'   coefficient for a factor level drifts away from its (flat) actual coefficient
#'   across the focus term, that level is not consistent with the focus trend.
#'
#'   This is the companion of [plot_implied_residuals()]: instead of a non-modelled
#'   variable compared against the focus trend, it examines a term already in the
#'   model and compares the residual-implied coefficient against the model's own
#'   fitted coefficient for that term.
#' @param obj A `gam_influence` object with calculated indices (see [calculate_influence()]).
#' @param term The modelled factor term to plot. Either the character name of the
#'   term (e.g. `"area"`) or its numeric index among the model terms.
#' @param ylim (Optional) The y-axis limits for the plot.
#' @param n.exclude (Optional) Minimum number of records required to include a
#'   `term-level x focus-level` cell in the plot.
#' @return A ggplot object (multi-panel plot, faceted by term level).
#' @details
#' The actual coefficient is taken from the model's fitted partial effects
#' (`obj$calculated$predictions`) and is constant across the focus term for each
#' level. The implied coefficient adds the mean standardised residual within each
#' `term-level x focus-level` cell, so it reveals how the residuals would shift
#' that level's coefficient at each focus level. A level that is consistent with
#' the focus trend has implied points scattered around its flat actual line; a
#' level that systematically departs from the line in some focus levels is
#' inconsistent with the focus trend.
#' @seealso [plot_implied_residuals()]
#' @importFrom dplyr group_by summarise n %>%
#' @importFrom ggplot2 ggplot aes geom_errorbar geom_point geom_line geom_col scale_colour_manual scale_fill_manual facet_wrap ylab xlab coord_cartesian
#' @importFrom rlang .data parse_expr
#' @importFrom stats sd formula
#' @export
plot_residual_implied_coefficients <- function(obj, term, ylim = NULL, n.exclude = 0) {
  # --- Checks ---
  if (is.null(obj$data)) stop("No data found in gam_influence object.", call. = FALSE)
  if (is.null(obj$calculated) || is.null(obj$calculated$predictions)) {
    stop("No calculated predictions found. Run calculate_influence() first.", call. = FALSE)
  }
  data <- obj$data
  focus_var <- obj$focus
  islog <- isTRUE(obj$islog)

  # --- Resolve term (allow numeric index, as in plot_cdi) ---
  if (is.numeric(term) && length(term) == 1 && term == as.integer(term)) {
    all_terms <- get_terms(obj, full = TRUE)
    if (term > length(all_terms) || term < 1) {
      stop("Term index out of bounds. There are ", length(all_terms), " valid terms:\n  ",
        paste(all_terms, collapse = "\n  "),
        call. = FALSE
      )
    }
    term <- all_terms[term]
  }

  # Identify the variable(s) behind the term and validate it is a modelled factor
  model_vars <- get_terms(obj, full = FALSE)
  term_vars <- all.vars(rlang::parse_expr(term))
  term_vars_in_model <- term_vars[term_vars %in% model_vars]
  if (length(term_vars_in_model) == 0) {
    stop(paste0(
      "None of the variables in the supplied term ('", term, "') match a model term.\n",
      "Valid terms are: ", paste(setdiff(model_vars, focus_var), collapse = ", ")
    ), call. = FALSE)
  }
  term_var <- term_vars_in_model[1]
  if (term_var == focus_var) {
    stop("Cannot plot residual-implied coefficients for the focus term itself.", call. = FALSE)
  }
  if (!term_var %in% names(data)) {
    stop(paste0("Term variable '", term_var, "' not found in the data."), call. = FALSE)
  }
  if (!is.factor(data[[term_var]])) {
    stop(paste0("Term '", term_var, "' must be a modelled factor term."), call. = FALSE)
  }

  # --- Actual fitted coefficient per record (constant within each level) ---
  preds_df <- obj$calculated$predictions
  term_cols <- find_term_columns(term, colnames(preds_df), group_by_bysmooth = TRUE)
  if (length(term_cols) == 0) {
    stop(paste0("Could not find fitted partial-effect columns for term '", term, "'."), call. = FALSE)
  }
  # On the link scale; summed across columns to mirror calculate_influence().
  actual_coef <- rowSums(preds_df[, term_cols, drop = FALSE])

  # --- Residuals on the model's modelled-response (link) scale ---
  # Evaluate the model's LHS so observed values are on the same scale as the
  # fitted values and the term coefficients (e.g. log(cpue) for islog models).
  # Using the raw response here would mismatch the scale and produce large,
  # one-sided residuals.
  predicted <- predict_response(obj$model, data)$fit
  observed <- as.numeric(eval(stats::formula(obj$model)[[2]], envir = data))
  residual <- observed - predicted

  records <- data.frame(
    level = data[[term_var]],
    focus_var = data[[focus_var]],
    actual_coef = actual_coef,
    residual = residual
  )

  # --- Aggregate per term-level x focus-level cell ---
  df <- records %>%
    dplyr::group_by(.data$level, .data$focus_var) %>%
    dplyr::summarise(
      mean_resid = mean(.data$residual, na.rm = TRUE),
      se_resid = sd(.data$residual, na.rm = TRUE) / sqrt(dplyr::n()),
      actual_coef = mean(.data$actual_coef, na.rm = TRUE),
      n_records = dplyr::n(),
      .groups = "drop"
    )

  if (n.exclude > 0) {
    df <- df[df$n_records >= n.exclude, ]
  }

  # --- Implied coefficient and uncertainty (link scale) ---
  df$implied_coef <- df$actual_coef + df$mean_resid
  df$q025 <- df$implied_coef - 1.96 * df$se_resid
  df$q975 <- df$implied_coef + 1.96 * df$se_resid

  # --- Scale handling: exponentiate to the coefficient scale shown elsewhere ---
  if (islog) {
    df$actual_coef <- exp(df$actual_coef)
    df$implied_coef <- exp(df$implied_coef)
    df$q025 <- exp(df$q025)
    df$q975 <- exp(df$q975)
  }

  # --- y-axis limits ---
  if (is.null(ylim)) {
    ylim <- range(c(df$q025, df$q975, df$actual_coef, df$implied_coef), na.rm = TRUE)
    if (islog) ylim[1] <- max(0, ylim[1])
  }

  # --- Count bars scaled into the lower portion of the panel ---
  y_range <- diff(ylim)
  df$N <- (df$n_records / max(df$n_records, na.rm = TRUE)) * (y_range * 0.3) + ylim[1]

  # Convert focus to numeric where possible for a continuous x-axis
  if (is.factor(df$focus_var) && all(!is.na(suppressWarnings(as.numeric(as.character(levels(df$focus_var))))))) {
    df$focus_var <- as.numeric(as.character(df$focus_var))
  }

  # Interpretation guidance
  message(
    "Interpreting the residual-implied coefficients plot:\n",
    "  - Each panel is a level of '", term_var, "'; the x-axis is the focus term '", focus_var, "'.\n",
    "  - The black 'Actual' line is the model's fitted coefficient for that level (flat: assumed\n",
    "    constant across the focus term).\n",
    "  - Blue 'Implied' points/intervals are that coefficient adjusted by the mean residual in\n",
    "    each cell.\n",
    "  - Where Implied points sit away from the flat Actual line (line outside the interval), that\n",
    "    level is inconsistent with the focus trend in those focus levels.\n",
    "  - Grey bars show the number of records per cell (wider support = more reliable points)."
  )

  # --- Plot ---
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$focus_var)) +
    ggplot2::geom_col(ggplot2::aes(y = .data$N, fill = "Records"), alpha = 0.4, width = 1) +
    ggplot2::geom_line(ggplot2::aes(y = .data$actual_coef, colour = "Actual", group = 1)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$q025, ymax = .data$q975, colour = "Implied"),
      alpha = 0.6, group = 1
    ) +
    ggplot2::geom_point(ggplot2::aes(y = .data$implied_coef, colour = "Implied")) +
    ggplot2::scale_colour_manual(name = NULL, values = c("Implied" = "royalblue", "Actual" = "black")) +
    ggplot2::scale_fill_manual(name = NULL, values = c("Records" = "grey60")) +
    ggplot2::facet_wrap(~ .data$level, scales = "fixed") +
    ggplot2::ylab("Implied coefficient residual") +
    ggplot2::xlab(focus_var) +
    ggplot2::coord_cartesian(ylim = ylim)

  return(p)
}

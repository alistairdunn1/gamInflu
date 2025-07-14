#' @title Plot CDI Coefficients
#' @description Plot CDI Coefficients
#' @param x A gam_influence object
#' @param term The term to plot
#' @return A ggplot object for the CDI plot
#' @keywords internal
#'
plot_cdi_coefficients <- function(x, term, ...) {
  # Get term information
  term_info <- x$expanded_terms[[term]]

  if (is.null(term_info)) {
    stop("Term '", term, "' not found in expanded terms")
  }

  # Get predictions for this smooth term
  pred_terms <- predict(x$model, type = "terms", se.fit = TRUE)
  smooth_idx <- term_info$smooth_index

  # Extract term effects and standard errors
  if (term_info$type == "by_smooth") {
    # Handle by-variables
    by_var <- term_info$by_variable
    by_level <- term_info$by_level

    # Subset to relevant observations
    relevant_obs <- x$model$model[[by_var]] == by_level
    term_effects <- pred_terms$fit[relevant_obs, smooth_idx]
    term_se <- pred_terms$se.fit[relevant_obs, smooth_idx]

    # Get the main variable values for this subset
    main_var <- term_info$variables[1]
    var_values <- x$data[[main_var]][relevant_obs]

    title_text <- paste0(term_info$original_label, " (", by_level, ")")
  } else {
    term_effects <- pred_terms$fit[, term_info$original_label]
    term_se <- pred_terms$se.fit[, term_info$original_label]

    # Get the main variable values
    main_var <- term_info$variables[1]
    var_values <- x$data[[main_var]]

    title_text <- term_info$original_label
  }

  # Create data frame for plotting
  coef_data <- data.frame(
    var_value = var_values,
    effect = term_effects,
    se = term_se,
    lower = term_effects - 2 * term_se,
    upper = term_effects + 2 * term_se
  )

  # Sort by variable value for better plotting
  coef_data <- coef_data[order(coef_data$var_value), ]

  # Create the plot
  if (is.numeric(var_values)) {
    # For numeric variables, use lines and ribbons
    p <- ggplot(coef_data, aes(x = var_value, y = effect)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "lightblue") +
      geom_line(color = "black", size = 1) +
      labs(x = main_var, y = "Relative effect")
  } else {
    # For factor variables, use points and error bars
    # Aggregate by factor level (take mean if multiple observations per level)
    coef_summary <- coef_data %>%
      group_by(var_value) %>%
      summarise(
        effect = mean(effect), se = sqrt(mean(se^2)), # Root mean square of standard errors
        lower = mean(lower), upper = mean(upper), .groups = "drop"
      ) %>%
      mutate(level_number = row_number())

    p <- ggplot(coef_summary, aes(x = level_number, y = effect)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.8) +
      geom_point(size = 3, color = "black") +
      scale_x_continuous(breaks = coef_summary$level_number, labels = coef_summary$var_value, position = "top") +
      labs(x = "", y = "Relative effect")
  }

  return(p)
}

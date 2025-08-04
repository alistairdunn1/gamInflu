#' @title Print Method for Combined Influence Objects
#' @description Print method for objects of class 'gam_influence_combined'
#' @param x A gam_influence_combined object
#' @param ... Additional arguments (unused)
#' @return Invisibly returns the input object
#' @export
print.gam_influence_combined <- function(x, ...) {
  cat("Combined GAM Influence Analysis (Delta-GLM)\n")
  cat("===========================================\n\n")

  cat("Focus Term:", x$focus_term, "\n")
  cat("Combination Method:", tools::toTitleCase(x$method), "\n")
  cat("Confidence Method:", tools::toTitleCase(x$confidence_method), "\n")
  cat("Number of Levels:", x$diagnostics$n_levels, "\n\n")

  cat("Component Models:\n")
  cat("  Binomial (Probability):", x$diagnostics$binomial_family, "\n")
  cat("  Positive (Conditional):", x$diagnostics$positive_family, "\n\n")

  cat("Combined Index Summary:\n")
  summary_stats <- summary(x$combined_indices$combined_index)
  print(summary_stats)

  cat(
    "\nMean Coefficient of Variation:",
    sprintf("%.3f", x$diagnostics$mean_combined_cv), "\n"
  )

  if (length(x$diagnostics$warnings) > 0) {
    cat("\nWarnings:\n")
    for (warning in x$diagnostics$warnings) {
      cat("  -", warning, "\n")
    }
  }

  cat("\nUse summary() for detailed results or plot() for visualisation.\n")

  invisible(x)
}

#' @title Summary Method for Combined Influence Objects
#' @description Detailed summary for objects of class 'gam_influence_combined'
#' @param object A gam_influence_combined object
#' @param ... Additional arguments (unused)
#' @return Invisibly returns a summary data frame
#' @export
summary.gam_influence_combined <- function(object, ...) {
  cat("Combined GAM Influence Analysis - Detailed Summary\n")
  cat("================================================\n\n")

  # Basic information
  cat("Analysis Configuration:\n")
  cat("  Focus Term:", object$focus_term, "\n")
  cat("  Combination Method:", tools::toTitleCase(object$method), "\n")
  cat("  Confidence Method:", tools::toTitleCase(object$confidence_method), "\n")
  cat("  Rescaled:", ifelse(object$rescaled, "Yes", "No"), "\n")
  cat("  Number of Levels:", object$diagnostics$n_levels, "\n\n")

  # Model information
  cat("Component Models:\n")
  cat("  Binomial Model Family:", object$diagnostics$binomial_family, "\n")
  cat("  Positive Model Family:", object$diagnostics$positive_family, "\n\n")

  # Component ranges
  cat("Component Index Ranges:\n")
  cat(
    "  Binomial (Probability): [",
    sprintf("%.3f", object$diagnostics$binomial_range[1]), ", ",
    sprintf("%.3f", object$diagnostics$binomial_range[2]), "]\n"
  )
  cat(
    "  Positive (Conditional): [",
    sprintf("%.3f", object$diagnostics$positive_range[1]), ", ",
    sprintf("%.3f", object$diagnostics$positive_range[2]), "]\n"
  )
  cat(
    "  Combined Index: [",
    sprintf("%.3f", object$diagnostics$combined_range[1]), ", ",
    sprintf("%.3f", object$diagnostics$combined_range[2]), "]\n\n"
  )

  # Uncertainty measures
  cat("Uncertainty Summary (CV):\n")
  cat("  Binomial Component:", sprintf("%.3f", object$diagnostics$mean_binom_cv), "\n")
  cat("  Positive Component:", sprintf("%.3f", object$diagnostics$mean_positive_cv), "\n")
  cat("  Combined Index:", sprintf("%.3f", object$diagnostics$mean_combined_cv), "\n\n")

  # Component correlation
  cat("Component Correlation:", sprintf("%.3f", object$diagnostics$component_correlation), "\n\n")

  # Detailed results table
  cat("Detailed Results by", object$focus_term, ":\n")
  results_table <- object$combined_indices[, c(
    "level",
    "standardised_index_binom",
    "standardised_index_pos",
    "combined_index",
    "combined_lower",
    "combined_upper",
    "combined_cv"
  )]

  # Round numeric columns for display
  numeric_cols <- sapply(results_table, is.numeric)
  results_table[numeric_cols] <- lapply(results_table[numeric_cols], function(x) round(x, 4))

  # Rename columns for better display
  names(results_table) <- c("Level", "Prob_Index", "Pos_Index", "Combined", "Lower_CI", "Upper_CI", "CV")

  print(results_table, row.names = FALSE)

  # Warnings
  if (length(object$diagnostics$warnings) > 0) {
    cat("\nWarnings:\n")
    for (warning in object$diagnostics$warnings) {
      cat("  -", warning, "\n")
    }
  }

  invisible(results_table)
}

#' @title Plot Method for Combined Influence Objects
#' @description Plot method for objects of class 'gam_influence_combined'
#' @param x A gam_influence_combined object
#' @param type Character. Type of plot:
#'   - "comparison": Single plot showing the components and combined (default)
#'   - "combined": Combined index with confidence intervals
#'   - "components": All three indices (binomial, positive, combined) on separate panels
#' @param show_points Logical. Should individual points be shown? Default TRUE.
#' @param show_ci Logical. Should confidence intervals be shown? Default TRUE.
#' @param ... Additional arguments passed to ggplot2
#' @return A ggplot object
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon labs scale_colour_manual facet_wrap
#' @importFrom patchwork wrap_plots
#' @importFrom rlang .data
#' @export
plot.gam_influence_combined <- function(x, type = c("comparison", "combined", "components"),
                                        show_points = TRUE, show_ci = TRUE, ...) {
  type <- match.arg(type)

  # Prepare data for plotting
  plot_data <- x$combined_indices

  # Convert level to numeric if possible for better plotting
  if (is.factor(plot_data$level) && all(!is.na(as.numeric(as.character(levels(plot_data$level)))))) {
    plot_data$level_numeric <- as.numeric(as.character(plot_data$level))
  } else {
    plot_data$level_numeric <- as.numeric(plot_data$level)
  }

  if (type == "combined") {
    # Plot only the combined index
    p <- create_combined_index_plot(plot_data, x$focus_term, x$method, show_points, show_ci)
  } else if (type == "components") {
    # Plot all three indices on separate panels
    p <- create_components_plot(plot_data, x$focus_term, x$method, show_points, show_ci)
  } else if (type == "comparison") {
    # Side-by-side comparison
    p <- create_comparison_plot(plot_data, x$focus_term, x$method, show_points, show_ci)
  }

  return(p)
}

#' @title Create Combined Index Plot
#' @description Internal function to create combined index plot
#' @param plot_data Data frame with plotting data
#' @param focus_term Name of focus term
#' @param method Combination method
#' @param show_points Logical for showing points
#' @param show_ci Logical for showing confidence intervals
#' @return ggplot object
#' @noRd
create_combined_index_plot <- function(plot_data, focus_term, method, show_points, show_ci) {
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$level_numeric, group = 1))

  # Add confidence ribbon if requested
  if (show_ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$combined_lower,
        ymax = .data$combined_upper
      ),
      alpha = 0.3, fill = "blue"
    )
  }

  # Add line
  p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$combined_index),
    colour = "royalblue"
  )

  # Add points if requested
  if (show_points) {
    p <- p + ggplot2::geom_point(ggplot2::aes(y = .data$combined_index),
      colour = "royalblue"
    )
  }

  # Add reference line
  p <- p + ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50")

  # Formatting
  p <- p + ggplot2::labs(x = focus_term, y = "Relative Index")

  return(p)
}

#' @title Create Components Plot
#' @description Internal function to create multi-panel components plot
#' @param plot_data Data frame with plotting data
#' @param focus_term Name of focus term
#' @param method Combination method
#' @param show_points Logical for showing points
#' @param show_ci Logical for showing confidence intervals
#' @return ggplot object
#' @noRd
create_components_plot <- function(plot_data, focus_term, method, show_points, show_ci) {
  # Reshape data for faceting
  long_data <- data.frame(
    level_numeric = rep(plot_data$level_numeric, 3),
    index = c(
      plot_data$standardised_index_binom,
      plot_data$standardised_index_pos,
      plot_data$combined_index
    ),
    lower = c(
      plot_data$stan_lower_binom,
      plot_data$stan_lower_pos,
      plot_data$combined_lower
    ),
    upper = c(
      plot_data$stan_upper_binom,
      plot_data$stan_upper_pos,
      plot_data$combined_upper
    ),
    component = factor(
      rep(c("Catch Probability", "Positive Catch Rate", "Combined Index"),
        each = nrow(plot_data)
      ),
      levels = c("Catch Probability", "Positive Catch Rate", "Combined Index")
    )
  )

  p <- ggplot2::ggplot(long_data, ggplot2::aes(x = .data$level_numeric, group = 1))

  # Add confidence ribbon if requested
  if (show_ci) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      alpha = 0.3, fill = "steelblue"
    )
  }

  # Add line
  p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$index), colour = "royalblue")

  # Add points if requested
  if (show_points) {
    p <- p + ggplot2::geom_point(ggplot2::aes(y = .data$index), colour = "royalblue")
  }

  # Add reference line
  p <- p + ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50")

  # Facet by component
  p <- p + ggplot2::facet_wrap(~component, scales = "free_y", ncol = 1)

  # Formatting
  p <- p + ggplot2::labs(x = focus_term, y = "Relative Index")

  return(p)
}

#' @title Create Comparison Plot
#' @description Internal function to create side-by-side comparison plot
#' @param plot_data Data frame with plotting data
#' @param focus_term Name of focus term
#' @param method Combination method
#' @param show_points Logical for showing points
#' @param show_ci Logical for showing confidence intervals
#' @return ggplot object
#' @noRd
create_comparison_plot <- function(plot_data, focus_term, method, show_points, show_ci) {
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$level_numeric, group = 1))

  # Add confidence ribbons if requested
  if (show_ci) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data$stan_lower_binom,
          ymax = .data$stan_upper_binom
        ),
        alpha = 0.2, fill = "#E74C3C"
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data$stan_lower_pos,
          ymax = .data$stan_upper_pos
        ),
        alpha = 0.2, fill = "#3498DB"
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data$combined_lower,
          ymax = .data$combined_upper
        ),
        alpha = 0.2, fill = "#2ECC71"
      )
  }

  # Add lines
  p <- p +
    ggplot2::geom_line(ggplot2::aes(y = .data$standardised_index_binom, colour = "Catch probability")) +
    ggplot2::geom_line(ggplot2::aes(y = .data$standardised_index_pos, colour = "Positive catch")) +
    ggplot2::geom_line(ggplot2::aes(y = .data$combined_index, colour = "Combined index"))

  # Add points if requested
  if (show_points) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(y = .data$standardised_index_binom, colour = "Catch probability")) +
      ggplot2::geom_point(ggplot2::aes(y = .data$standardised_index_pos, colour = "Positive catch")) +
      ggplot2::geom_point(ggplot2::aes(y = .data$combined_index, colour = "Combined index"))
  }

  # Add reference line
  p <- p + ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50")

  # Colour scale
  p <- p + ggplot2::scale_colour_manual(
    name = "Component",
    values = c("Catch probability" = "#E74C3C", "Positive catch" = "#3498DB", "Combined index" = "#2ECC71"),
    breaks = c("Catch probability", "Positive catch", "Combined index")
  )

  # Formatting
  p <- p +
    ggplot2::labs(x = focus_term, y = "Relative Index")

  return(p)
}

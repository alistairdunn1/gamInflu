#' @name gam_influence_comparison_methods
#' @title Print Method for gam_influence_comparison
#' @description Prints a summary of the grouped influence analysis results.
#' @param x A gam_influence_comparison object.
#' @param ... Additional arguments (not used).
#' @export
print.gam_influence_comparison <- function(x, ...) {
  cat("GAM Influence Comparison Analysis\n")
  cat("=================================\n\n")
  cat("Focus variable:", x$focus, "\n")
  cat("Grouping variable:", x$grouping_var, "\n")
  cat("Number of groups:", length(x$groups), "\n")
  cat("Groups:", paste(x$groups, collapse = ", "), "\n")
  cat("Model family:", x$family$family, "\n")
  cat("Formula:", deparse(x$formula), "\n\n")

  # Show sample sizes and focus levels per group
  cat("Group Details:\n")
  for (group in x$groups) {
    result <- x$results[[group]]
    n_obs <- nrow(result$data)
    n_focus_levels <- length(unique(result$data[[x$focus]]))
    cat("  ", group, ": ", n_obs, " observations, ", n_focus_levels, " focus levels\n", sep = "")
  }

  cat("\nUse plot() to visualize results or summary() for detailed statistics.\n")
}

#' @title Summary Method for gam_influence_comparison
#' @description Provides detailed summary statistics for grouped influence analysis.
#' @param object A gam_influence_comparison object.
#' @param ... Additional arguments (not used).
#' @export
summary.gam_influence_comparison <- function(object, ...) {
  cat("GAM Influence Comparison Analysis - Detailed Summary\n")
  cat("===================================================\n\n")

  print(object)

  cat("\nModel Fit Statistics by Group:\n")
  cat("------------------------------\n")

  for (group in x$groups) {
    result <- x$results[[group]]

    cat("\nGroup:", group, "\n")

    # Model fit statistics
    if ("summary" %in% names(result$calculated)) {
      summary_stats <- result$calculated$summary
      if (nrow(summary_stats) > 0) {
        final_stats <- summary_stats[nrow(summary_stats), ]
        cat("  R-squared:", round(final_stats$r_sq, 4), "\n")
        cat("  Deviance explained:", round(final_stats$deviance_explained, 4), "\n")
        cat("  AIC:", round(final_stats$aic, 2), "\n")
      }
    }

    # Focus term index range
    if ("indices" %in% names(result$calculated)) {
      indices <- result$calculated$indices
      if ("standardised_index" %in% names(indices)) {
        index_range <- range(indices$standardised_index, na.rm = TRUE)
        cat("  Index range:", round(index_range[1], 3), "to", round(index_range[2], 3), "\n")
      }
    }
  }

  invisible(object)
}

#' @title Plot Method for gam_influence_comparison
#' @description Creates comparative plots for grouped influence analysis results.
#' @param x A gam_influence_comparison object.
#' @param type Character string specifying plot type: "standardisation", "comparison",
#'   "stepwise", "influence", or "terms".
#' @param groups Character vector of specific groups to plot (default: all groups).
#' @param ncol Number of columns for multi-panel plots.
#' @param ... Additional arguments passed to individual plotting functions.
#' @return A ggplot or patchwork object.
#' @details
#' Plot types:
#' - "standardisation": Standardised vs unstandardised indices for each group
#' - "comparison": Direct comparison of standardised indices across groups
#' - "stepwise": Stepwise index progression for each group
#' - "influence": Term influence plots for each group
#' - "terms": Model term effects for each group
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon labs theme_minimal facet_wrap
#' @export
plot.gam_influence_comparison <- function(x, type = c("standardisation", "comparison", "stepwise", "influence", "terms"),
                                          groups = NULL, ncol = NULL, ...) {
  type <- match.arg(type)

  if (is.null(groups)) {
    groups <- x$groups
  } else {
    # Validate requested groups
    missing_groups <- setdiff(groups, x$groups)
    if (length(missing_groups) > 0) {
      stop("Groups not found: ", paste(missing_groups, collapse = ", "))
    }
  }

  if (is.null(ncol)) {
    ncol <- min(3, length(groups))
  }

  if (type == "comparison") {
    # Create combined comparison plot
    plot_comparison_indices(x, groups = groups, ...)
  } else {
    # Create individual plots for each group
    plots <- list()

    for (group in groups) {
      result <- x$results[[group]]

      p <- switch(type,
        "standardisation" = plot_standardisation(result, ...) +
          ggplot2::ggtitle(paste("Group:", group)),
        "stepwise" = plot_stepwise_index(result, ...) +
          ggplot2::ggtitle(paste("Group:", group)),
        "influence" = plot_term_influence(result, ...) +
          ggplot2::ggtitle(paste("Group:", group)),
        "terms" = plot_terms(result, ...) +
          ggplot2::ggtitle(paste("Group:", group))
      )

      plots[[group]] <- p
    }

    if (length(plots) == 1) {
      return(plots[[1]])
    } else {
      return(patchwork::wrap_plots(plots, ncol = ncol))
    }
  }
}

#' @title Plot Comparison of Standardised Indices
#' @description Creates a combined plot comparing standardised indices across groups.
#' @param object A gam_influence_comparison object.
#' @param groups Character vector of groups to include.
#' @param confidence_intervals Logical indicating whether to show confidence intervals.
#' @param ... Additional arguments (not used).
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon labs scale_colour_viridis_d theme_minimal
#' @importFrom dplyr bind_rows
#' @noRd
plot_comparison_indices <- function(object, groups = NULL, confidence_intervals = TRUE, ...) {
  if (is.null(groups)) {
    groups <- object$groups
  }

  # Combine standardised indices from all groups
  combined_data <- list()

  for (group in groups) {
    result <- object$results[[group]]
    indices <- result$calculated$indices

    if ("standardised_index" %in% names(indices)) {
      df <- data.frame(
        level = indices$level,
        index = indices$standardised_index,
        group = group,
        stringsAsFactors = FALSE
      )

      # Add confidence intervals if available
      if (confidence_intervals && "stan_lower" %in% names(indices)) {
        df$lower <- indices$stan_lower
        df$upper <- indices$stan_upper
      }

      combined_data[[group]] <- df
    }
  }

  if (length(combined_data) == 0) {
    stop("No standardised indices found in results")
  }

  plot_data <- dplyr::bind_rows(combined_data)

  # Convert level to numeric if possible
  if (is.factor(plot_data$level) && all(!is.na(as.numeric(as.character(levels(plot_data$level)))))) {
    plot_data$level <- as.numeric(as.character(plot_data$level))
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = level, y = index, colour = group, fill = group)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_viridis_d(name = object$grouping_var) +
    ggplot2::labs(
      x = object$focus,
      y = "Standardised Index",
    )

  # Add confidence intervals if available
  if (confidence_intervals && "lower" %in% names(plot_data)) {
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2) +
      ggplot2::scale_fill_viridis_d(name = object$grouping_var, guide = "none")
  }

  return(p)
}

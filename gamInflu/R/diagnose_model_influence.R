#' @title Comprehensive Model Diagnostics and Influence Summary
#' @description Provides a comprehensive diagnostic summary of GAM model adequacy and influence analysis results.
#' @param obj A gam_influence object that has been processed with calculate_influence.
#' @param residual_type Character. Type of residuals to use for diagnostics.
#' @param qq_test_threshold Numeric. P-value threshold for normality tests (default 0.05).
#' @param funnel_threshold Numeric. Threshold for detecting funnel-shaped residuals (default 2.0).
#' @param influence_threshold Numeric. Threshold for identifying influential years/levels as deviations from 1 (default 0.1).
#' @param index_difference_threshold Numeric. Threshold for flagging large differences (default 0.2).
#' @param detailed_output Logical. If TRUE, provides detailed diagnostic information (default TRUE).
#' @return A list containing diagnostic results
#' @importFrom stats residuals fitted quantile var shapiro.test
#' @export
diagnose_model_influence <- function(obj,
                                     residual_type = c("deviance", "pearson", "response", "working"),
                                     qq_test_threshold = 0.05,
                                     funnel_threshold = 2.0,
                                     influence_threshold = 0.1,
                                     index_difference_threshold = 0.2,
                                     detailed_output = TRUE) {
  # Validate inputs
  if (!inherits(obj, "gam_influence")) {
    stop("Object must be of class 'gam_influence'", call. = FALSE)
  }

  if (is.null(obj$calculated) || isFALSE(obj$calculated)) {
    stop("Calculations not performed. Please run calculate_influence() first.", call. = FALSE)
  }

  residual_type <- match.arg(residual_type)

  # Extract components
  model <- obj$model
  influences <- obj$calculated$influences
  indices <- obj$calculated$indices

  # Calculate residuals and fitted values
  residuals_vec <- residuals(model, type = residual_type)
  fitted_vals <- fitted(model)

  # === RESIDUAL ANALYSIS ===
  fitted_quantiles <- quantile(fitted_vals, probs = c(0.25, 0.75), na.rm = TRUE)
  low_fitted <- fitted_vals <= fitted_quantiles[1]
  high_fitted <- fitted_vals >= fitted_quantiles[2]

  var_low <- var(residuals_vec[low_fitted], na.rm = TRUE)
  var_high <- var(residuals_vec[high_fitted], na.rm = TRUE)
  variance_ratio <- max(var_high, var_low) / min(var_high, var_low)

  outliers <- abs(residuals_vec) > 2 * sd(residuals_vec, na.rm = TRUE)
  n_outliers <- sum(outliers, na.rm = TRUE)

  residual_diagnostics <- list(
    funnel_shape = list(
      detected = variance_ratio > funnel_threshold,
      variance_ratio = variance_ratio
    ),
    summary_stats = list(
      n_outliers = n_outliers,
      skewness = mean((residuals_vec - mean(residuals_vec, na.rm = TRUE))^3, na.rm = TRUE) /
        (sd(residuals_vec, na.rm = TRUE)^3),
      kurtosis = mean((residuals_vec - mean(residuals_vec, na.rm = TRUE))^4, na.rm = TRUE) /
        (sd(residuals_vec, na.rm = TRUE)^4) - 3
    )
  )

  # === DISTRIBUTIONAL TESTS ===
  normality_p <- if (model$family$family == "gaussian" && length(residuals_vec) >= 3 && length(residuals_vec) <= 5000) {
    tryCatch(shapiro.test(residuals_vec)$p.value, error = function(e) NA)
  } else {
    NA
  }

  distributional_tests <- list(
    normality = list(
      p_value = normality_p,
      assumption_violated = if (!is.na(normality_p)) normality_p < qq_test_threshold else FALSE
    ),
    family_specific = list()
  )

  # === INDEX DIFFERENCES ===
  if ("unstandardised_index" %in% names(indices) && "standardised_index" %in% names(indices)) {
    rel_diff <- abs(indices$standardised_index - indices$unstandardised_index) /
      abs(indices$unstandardised_index)
    flagged <- rel_diff > index_difference_threshold & !is.na(rel_diff)

    index_differences <- list(
      summary = list(
        n_flagged = sum(flagged, na.rm = TRUE),
        max_difference = max(rel_diff, na.rm = TRUE),
        mean_difference = mean(rel_diff, na.rm = TRUE),
        flagged_levels = if ("level" %in% names(indices)) indices$level[flagged] else character(0)
      )
    )
  } else {
    index_differences <- list(summary = list(n_flagged = NA))
  }

  # === INFLUENCE IMPACT ===
  if ("influence" %in% names(influences)) {
    # Calculate deviations from 1 (baseline = no influence)
    influence_deviations <- abs(influences$influence - 1)
    high_influence_idx <- influence_deviations > influence_threshold

    influence_impact <- list(
      high_influence_count = sum(high_influence_idx, na.rm = TRUE),
      total_observations = length(influence_deviations),
      high_influence_prop = mean(high_influence_idx, na.rm = TRUE),
      mean_abs_deviation = mean(influence_deviations, na.rm = TRUE),
      max_abs_deviation = max(influence_deviations, na.rm = TRUE)
    )

    # Focus values with high influence
    if (influence_impact$high_influence_count > 0 && "level" %in% names(influences)) {
      influence_impact$high_influence_focus <- influences$level[high_influence_idx]
    }

    # Term-by-term analysis
    if ("term" %in% names(influences)) {
      unique_terms <- unique(influences$term)
      term_details <- list()

      for (term in unique_terms) {
        term_data <- influences[influences$term == term, ]
        # Calculate deviations from 1 (baseline = no influence)
        term_influence_deviations <- abs(term_data$influence - 1)
        term_high_idx <- term_influence_deviations > influence_threshold

        term_stats <- list(
          term = term,
          mean_abs_deviation = mean(term_influence_deviations, na.rm = TRUE),
          max_abs_deviation = max(term_influence_deviations, na.rm = TRUE),
          n_observations = nrow(term_data),
          n_high_influence = sum(term_high_idx, na.rm = TRUE),
          prop_high_influence = mean(term_high_idx, na.rm = TRUE)
        )

        # Focus values with high influence for this term
        if (any(term_high_idx) && "level" %in% names(term_data)) {
          high_focus_data <- term_data[term_high_idx, ]
          # Order by absolute deviation from 1
          high_focus_data <- high_focus_data[order(abs(high_focus_data$influence - 1), decreasing = TRUE), ]

          term_stats$high_influence_focus_values <- high_focus_data$level
          term_stats$strongest_influence_focus <- high_focus_data$level[1]
          term_stats$strongest_influence_value <- high_focus_data$influence[1]
        } else {
          term_stats$high_influence_focus_values <- character(0)
          term_stats$strongest_influence_focus <- NA
          term_stats$strongest_influence_value <- NA
        }

        # Assessment - base on actual threshold exceedances, not just mean
        if (term_stats$n_high_influence > 0) {
          term_stats$assessment <- if (term_stats$prop_high_influence > 0.5) {
            "High influence (affects many focus levels)"
          } else {
            "Moderate influence (affects some focus levels)"
          }
        } else {
          term_stats$assessment <- "Low influence (no focus values exceed threshold)"
        }

        term_details[[term]] <- term_stats
      }

      # Sort terms by mean absolute deviation from 1
      mean_deviations <- sapply(term_details, function(x) x$mean_abs_deviation)
      term_order <- order(mean_deviations, decreasing = TRUE)
      influence_impact$term_details <- term_details[term_order]
    }
  } else {
    influence_impact <- list(note = "Influence column not found")
  }

  # === MODEL ADEQUACY ===
  issues <- character(0)
  recommendations <- character(0)

  if (residual_diagnostics$funnel_shape$detected) {
    issues <- c(issues, "Heteroscedasticity detected")
    recommendations <- c(recommendations, "Consider variance modeling or different family")
  }

  if (!is.na(distributional_tests$normality$assumption_violated) && distributional_tests$normality$assumption_violated) {
    issues <- c(issues, "Normality assumption violated")
    recommendations <- c(recommendations, "Consider non-Gaussian family or transformation")
  }

  if (n_outliers > 0) {
    issues <- c(issues, paste(n_outliers, "potential outliers"))
    recommendations <- c(recommendations, "Investigate outlying observations")
  }

  model_adequacy <- list(
    overall_assessment = if (length(issues) == 0) "Good" else if (length(issues) <= 2) "Moderate concerns" else "Substantial concerns",
    identified_issues = if (length(issues) > 0) issues else "No major issues detected",
    recommendations = if (length(recommendations) > 0) recommendations else "Model appears adequate",
    n_issues = length(issues)
  )

  # === RESULTS ===
  result <- list(
    residual_diagnostics = residual_diagnostics,
    distributional_tests = distributional_tests,
    index_differences = index_differences,
    influence_impact = influence_impact,
    model_adequacy = model_adequacy,
    diagnostic_plots = if (detailed_output) list(note = "Diagnostic plots available with ggplot2") else NULL
  )

  class(result) <- "gam_influence_diagnostics"
  return(result)
}

#' @title Print method for GAM influence diagnostics
#' @description Print summary of GAM influence diagnostics
#' @param x A gam_influence_diagnostics object
#' @param ... Additional arguments (not used)
#' @return Invisibly returns the input object
#' @export
print.gam_influence_diagnostics <- function(x, ...) {
  cat("GAM Influence Model Diagnostics Summary\n")
  cat("=======================================\n\n")

  # Overall Assessment
  cat("OVERALL ASSESSMENT:", x$model_adequacy$overall_assessment, "\n")
  cat("Issues identified:", x$model_adequacy$n_issues, "\n\n")

  # Residual Diagnostics
  cat("RESIDUAL PATTERN ANALYSIS:\n")
  cat("- Funnel shape (heteroscedasticity):", ifelse(x$residual_diagnostics$funnel_shape$detected, "DETECTED", "Not detected"), "\n")
  if (x$residual_diagnostics$funnel_shape$detected) {
    cat("  Variance ratio:", round(x$residual_diagnostics$funnel_shape$variance_ratio, 2), "\n")
  }
  cat("- Outliers (>2 SD):", x$residual_diagnostics$summary_stats$n_outliers, "\n")
  cat("- Residual skewness:", round(x$residual_diagnostics$summary_stats$skewness, 3), "\n")
  cat("- Residual kurtosis:", round(x$residual_diagnostics$summary_stats$kurtosis, 3), "\n\n")

  # Distributional Tests
  cat("DISTRIBUTIONAL ASSUMPTIONS:\n")
  if (!is.na(x$distributional_tests$normality$p_value)) {
    cat("- Normality test p-value:", round(x$distributional_tests$normality$p_value, 4), "\n")
    cat("- Normal assumption:", ifelse(x$distributional_tests$normality$assumption_violated, "VIOLATED", "Satisfied"), "\n")
  }
  cat("\n")

  # Index Differences
  if (!is.na(x$index_differences$summary$n_flagged)) {
    cat("INDEX DIFFERENCES:\n")
    cat("- Levels with large differences:", x$index_differences$summary$n_flagged, "\n")
    cat("- Maximum relative difference:", round(x$index_differences$summary$max_difference, 3), "\n")
    cat("- Mean relative difference:", round(x$index_differences$summary$mean_difference, 3), "\n")
    if (length(x$index_differences$summary$flagged_levels) > 0) {
      cat("- Flagged levels:", paste(x$index_differences$summary$flagged_levels, collapse = ", "), "\n")
    }
    cat("\n")
  }

  # Influence Impact
  if (!is.null(x$influence_impact) && !"note" %in% names(x$influence_impact)) {
    cat("INFLUENCE IMPACT ANALYSIS:\n")
    cat(
      "- High influence observations:", x$influence_impact$high_influence_count,
      "out of", x$influence_impact$total_observations,
      paste0("(", round(100 * x$influence_impact$high_influence_prop, 1), "%)"), "\n"
    )

    if (x$influence_impact$high_influence_count > 0) {
      cat("- Mean absolute deviation from 1:", round(x$influence_impact$mean_abs_deviation, 3), "\n")
      cat("- Maximum absolute deviation from 1:", round(x$influence_impact$max_abs_deviation, 3), "\n")
      if (!is.null(x$influence_impact$high_influence_focus) && length(x$influence_impact$high_influence_focus) > 0) {
        if (is.factor(x$influence_impact$high_influence_focus)) {
          unique_levels <- unique(x$influence_impact$high_influence_focus)
          focus_range <- paste(levels(unique_levels)[c(1, length(levels(unique_levels)))], collapse = " to ")
          cat("- Focus variable levels with high influence:", paste(unique_levels, collapse = ", "), "\n")
        } else {
          focus_range <- range(x$influence_impact$high_influence_focus, na.rm = TRUE)
          cat("- Focus variable range of high influence:", paste(focus_range, collapse = " to "), "\n")
        }
      }
    }

    # Display term summaries
    if (!is.null(x$influence_impact$term_details)) {
      cat("\nTERM-BY-TERM INFLUENCE SUMMARY:\n")
      n_terms_to_show <- length(x$influence_impact$term_details) # Show all terms

      for (i in 1:n_terms_to_show) {
        term_info <- x$influence_impact$term_details[[i]]
        cat("- Term:", term_info$term, "-", term_info$assessment, "\n")

        # Only show detailed diagnostics for terms with influence above threshold
        if (term_info$n_high_influence > 0) {
          cat(
            "  Mean |deviation from 1|:", round(term_info$mean_abs_deviation, 4),
            "| High influence obs.:", term_info$n_high_influence, "/", term_info$n_observations, "\n"
          )

          if (length(term_info$high_influence_focus_values) > 0) {
            focus_display <- if (length(term_info$high_influence_focus_values) <= 3) {
              paste(term_info$high_influence_focus_values, collapse = ", ")
            } else {
              paste0(
                paste(head(term_info$high_influence_focus_values, 3), collapse = ", "),
                " (+", length(term_info$high_influence_focus_values) - 3, " more)"
              )
            }
            cat("  High influence focus values:", focus_display, "\n")

            if (!is.na(term_info$strongest_influence_focus)) {
              cat(
                "  Strongest:", term_info$strongest_influence_focus,
                "(", round(term_info$strongest_influence_value, 4), ")\n"
              )
            }
          }
        }
      }
    }
    cat("\n")
  }

  # Recommendations
  cat("RECOMMENDATIONS:\n")
  for (rec in x$model_adequacy$recommendations) {
    cat("- ", rec, "\n")
  }
  cat("\n")

  if (!is.null(x$diagnostic_plots)) {
    cat("Diagnostic plots available in $diagnostic_plots\n")
  }

  invisible(x)
}

#' @title Summary method for GAM influence diagnostics
#' @description Detailed summary of GAM influence diagnostics
#' @param object A gam_influence_diagnostics object
#' @param ... Additional arguments (not used)
#' @return Invisibly returns the input object
#' @export
summary.gam_influence_diagnostics <- function(object, ...) {
  # First print the regular summary
  print(object)

  # Then add detailed term information if available
  if (!is.null(object$influence_impact$term_details)) {
    cat("\nDETAILED TERM-BY-TERM INFLUENCE ANALYSIS:\n")
    cat("=========================================\n\n")

    for (i in seq_along(object$influence_impact$term_details)) {
      term_info <- object$influence_impact$term_details[[i]]
      cat("TERM:", term_info$term, "\n")
      cat("Assessment:", term_info$assessment, "\n")

      # Only show detailed diagnostics for terms with influence above threshold
      if (term_info$n_high_influence > 0) {
        cat("Statistics:\n")
        cat("  - Mean |deviation from 1|:", round(term_info$mean_abs_deviation, 4), "\n")
        cat("  - Max |deviation from 1|:", round(term_info$max_abs_deviation, 4), "\n")
        cat(
          "  - High influence observations:", term_info$n_high_influence, "out of", term_info$n_observations,
          paste0("(", round(100 * term_info$prop_high_influence, 1), "%)"), "\n"
        )

        if (length(term_info$high_influence_focus_values) > 0) {
          cat("High influence focus values:", paste(term_info$high_influence_focus_values, collapse = ", "), "\n")
          cat(
            "Strongest influence:", term_info$strongest_influence_focus,
            "with influence =", round(term_info$strongest_influence_value, 4), "\n"
          )
        }
      }

      cat("\n")
    }
  }

  invisible(object)
}

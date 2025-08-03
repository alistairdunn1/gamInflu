#' @title Analyse Residual Patterns for Model Adequacy
#' @description Analyses patterns in GAM residuals by fitting models using other variables
#' in the dataset. This helps identify potential covariates that should be considered for
#' inclusion in the original GAM model. The function fits both linear and smooth relationships
#' between residuals and candidate variables, ranking them by explanatory power.
#' @param obj A `gam_influence` object that has been processed by `calculate_influence()`.
#' @param candidate_vars Character vector of variable names to test as potential predictors
#' of residual patterns. If NULL, all numeric variables not already in the model will be tested.
#' @param exclude_vars Character vector of variable names to exclude from analysis (e.g., IDs, dates).
#' @param residual_type Character; type of residuals to use. One of "deviance" (default),
#' "pearson", "response", or "working".
#' @param smooth_terms Logical; if TRUE (default), test smooth terms s(var) in addition to linear terms.
#' @param significance_level Numeric; significance level for determining important relationships (default 0.05).
#' @param min_r_squared Numeric; minimum R-squared threshold for flagging relationships (default 0.1).
#' @param plot_results Logical; if TRUE (default), generate diagnostic plots for significant relationships.
#' @return A list containing:
#'   - `linear_results`: Data frame of linear model results ranked by significance
#'   - `smooth_results`: Data frame of smooth model results (if smooth_terms = TRUE)
#'   - `significant_vars`: Character vector of variables showing significant patterns
#'   - `recommendations`: Text recommendations for model improvement
#'   - `residual_data`: Data frame with residuals and candidate variables
#'   - `plots`: List of diagnostic plots (if plot_results = TRUE)
#' @importFrom stats lm anova glm
#' @importFrom mgcv gam s
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs facet_wrap
#' @export
analyse_residual_patterns <- function(obj,
                                      candidate_vars = NULL,
                                      exclude_vars = NULL,
                                      residual_type = "deviance",
                                      smooth_terms = TRUE,
                                      significance_level = 0.05,
                                      min_r_squared = 0.1,
                                      plot_results = TRUE) {
  # Check if calculations have been performed
  if (is.null(obj$calculated)) {
    stop("Calculations not performed. Please run `calculate_influence()` first.", call. = FALSE)
  }

  # Extract model and data
  model <- obj$model
  data <- obj$data

  # Calculate residuals
  residuals_vec <- switch(residual_type,
    "deviance" = residuals(model, type = "deviance"),
    "pearson" = residuals(model, type = "pearson"),
    "response" = residuals(model, type = "response"),
    "working" = residuals(model, type = "working"),
    stop("Invalid residual_type. Choose from 'deviance', 'pearson', 'response', 'working'.", call. = FALSE)
  )

  # Get variables already in the model
  model_vars <- all.vars(formula(model))
  response_var <- model_vars[1]
  predictor_vars <- model_vars[-1]

  # Determine candidate variables
  if (is.null(candidate_vars)) {
    # Use all numeric variables not already in the model
    numeric_vars <- sapply(data, is.numeric)
    candidate_vars <- names(data)[numeric_vars & !names(data) %in% c(predictor_vars, response_var)]
  }

  # Remove excluded variables
  if (!is.null(exclude_vars)) {
    candidate_vars <- candidate_vars[!candidate_vars %in% exclude_vars]
  }

  # Remove variables with insufficient variation or too many missing values
  candidate_vars <- candidate_vars[sapply(candidate_vars, function(var) {
    if (!var %in% names(data)) {
      return(FALSE)
    }
    var_data <- data[[var]]
    if (sum(!is.na(var_data)) < 10) {
      return(FALSE)
    } # Need at least 10 observations
    if (is.numeric(var_data) && length(unique(var_data[!is.na(var_data)])) < 3) {
      return(FALSE)
    } # Need variation
    return(TRUE)
  })]

  if (length(candidate_vars) == 0) {
    stop("No suitable candidate variables found for analysis.", call. = FALSE)
  }

  # Create analysis dataset
  residual_data <- data.frame(
    residuals = residuals_vec,
    data[candidate_vars],
    stringsAsFactors = FALSE
  )

  # Remove rows with missing data
  complete_cases <- complete.cases(residual_data)
  residual_data <- residual_data[complete_cases, ]

  if (nrow(residual_data) < 20) {
    stop("Insufficient complete cases for analysis (need at least 20).", call. = FALSE)
  }

  # Initialize results
  linear_results <- data.frame(
    variable = character(),
    r_squared = numeric(),
    f_statistic = numeric(),
    p_value = numeric(),
    coefficient = numeric(),
    std_error = numeric(),
    significant = logical(),
    stringsAsFactors = FALSE
  )

  smooth_results <- data.frame(
    variable = character(),
    r_squared = numeric(),
    edf = numeric(),
    f_statistic = numeric(),
    p_value = numeric(),
    significant = logical(),
    stringsAsFactors = FALSE
  )

  # Test linear relationships
  message("Testing linear relationships...")
  for (var in candidate_vars) {
    tryCatch(
      {
        # Fit linear model
        lm_formula <- as.formula(paste("residuals ~", var))
        lm_fit <- lm(lm_formula, data = residual_data)
        lm_summary <- summary(lm_fit)

        # Extract results
        r_sq <- lm_summary$r.squared
        f_stat <- lm_summary$fstatistic[1]
        p_val <- pf(f_stat, lm_summary$fstatistic[2], lm_summary$fstatistic[3], lower.tail = FALSE)
        coef_val <- coef(lm_fit)[2]
        std_err <- lm_summary$coefficients[2, 2]

        linear_results <- rbind(linear_results, data.frame(
          variable = var,
          r_squared = r_sq,
          f_statistic = f_stat,
          p_value = p_val,
          coefficient = coef_val,
          std_error = std_err,
          significant = p_val < significance_level && r_sq > min_r_squared,
          stringsAsFactors = FALSE
        ))
      },
      error = function(e) {
        message("Warning: Could not fit linear model for variable '", var, "': ", e$message)
      }
    )
  }

  # Test smooth relationships (if requested)
  if (smooth_terms) {
    message("Testing smooth relationships...")
    for (var in candidate_vars) {
      tryCatch(
        {
          # Fit GAM with smooth term
          gam_formula <- as.formula(paste("residuals ~ s(", var, ")"))
          gam_fit <- mgcv::gam(gam_formula, data = residual_data)
          gam_summary <- summary(gam_fit)

          # Extract results
          r_sq <- gam_summary$r.sq
          smooth_terms_summary <- gam_summary$s.table
          if (nrow(smooth_terms_summary) > 0) {
            edf_val <- smooth_terms_summary[1, "edf"]
            f_stat <- smooth_terms_summary[1, "F"]
            p_val <- smooth_terms_summary[1, "p-value"]
          } else {
            edf_val <- NA
            f_stat <- NA
            p_val <- 1
          }

          smooth_results <- rbind(smooth_results, data.frame(
            variable = var,
            r_squared = r_sq,
            edf = edf_val,
            f_statistic = f_stat,
            p_value = p_val,
            significant = !is.na(p_val) && p_val < significance_level && r_sq > min_r_squared,
            stringsAsFactors = FALSE
          ))
        },
        error = function(e) {
          message("Warning: Could not fit smooth model for variable '", var, "': ", e$message)
        }
      )
    }
  }

  # Sort results by R-squared (descending)
  linear_results <- linear_results[order(-linear_results$r_squared), ]
  if (smooth_terms) {
    smooth_results <- smooth_results[order(-smooth_results$r_squared), ]
  }

  # Identify significant variables
  significant_linear <- linear_results$variable[linear_results$significant]
  significant_smooth <- if (smooth_terms) smooth_results$variable[smooth_results$significant] else character()
  significant_vars <- unique(c(significant_linear, significant_smooth))

  # Generate recommendations
  recommendations <- generate_recommendations(
    linear_results, smooth_results, significant_vars,
    min_r_squared, significance_level
  )

  # Create diagnostic plots
  plots <- NULL
  if (plot_results && length(significant_vars) > 0) {
    plots <- create_residual_pattern_plots(residual_data, significant_vars, residual_type)
  }

  # Prepare results
  results <- list(
    linear_results = linear_results,
    smooth_results = if (smooth_terms) smooth_results else NULL,
    significant_vars = significant_vars,
    recommendations = recommendations,
    residual_data = residual_data,
    plots = plots,
    analysis_info = list(
      residual_type = residual_type,
      n_observations = nrow(residual_data),
      n_candidate_vars = length(candidate_vars),
      significance_level = significance_level,
      min_r_squared = min_r_squared
    )
  )

  class(results) <- "residual_pattern_analysis"
  return(results)
}

#' @title Generate Recommendations for Model Improvement
#' @description Internal function to generate text recommendations based on residual pattern analysis
#' @param linear_results Data frame of linear model results
#' @param smooth_results Data frame of smooth model results
#' @param significant_vars Character vector of significant variables
#' @param min_r_squared Minimum R-squared threshold
#' @param significance_level Significance level used
#' @return Character vector of recommendations
#' @noRd
generate_recommendations <- function(linear_results, smooth_results, significant_vars,
                                     min_r_squared, significance_level) {
  recommendations <- character()

  if (length(significant_vars) == 0) {
    recommendations <- c(
      paste(
        "[OK] No significant patterns detected in residuals at alpha =", significance_level,
        "and R^2 >", min_r_squared
      ),
      "[OK] Current model appears to capture the main patterns in the data adequately.",
      "Consider checking:",
      "  - Model diagnostics (QQ plots, residual plots)",
      "  - Potential interaction terms between existing variables",
      "  - Temporal or spatial autocorrelation if applicable"
    )
  } else {
    recommendations <- c(
      paste("[WARNING] Significant residual patterns detected for", length(significant_vars), "variable(s):"),
      paste("   -", significant_vars),
      "",
      "[ANALYSIS] Recommended model improvements:"
    )

    for (var in significant_vars) {
      linear_r2 <- linear_results$r_squared[linear_results$variable == var][1]
      smooth_r2 <- if (!is.null(smooth_results)) {
        smooth_results$r_squared[smooth_results$variable == var][1]
      } else {
        NA
      }

      if (!is.na(smooth_r2) && smooth_r2 > linear_r2 + 0.05) {
        edf <- smooth_results$edf[smooth_results$variable == var][1]
        recommendations <- c(
          recommendations,
          paste("  * Add smooth term: s(", var, ") - shows non-linear pattern (EDF=",
            round(edf, 2), ", R^2=", round(smooth_r2, 3), ")",
            sep = ""
          )
        )
      } else if (!is.na(linear_r2)) {
        coef <- linear_results$coefficient[linear_results$variable == var][1]
        direction <- if (coef > 0) "positive" else "negative"
        recommendations <- c(
          recommendations,
          paste(
            "  * Add linear term:", var, "- shows", direction, "relationship (R^2=",
            round(linear_r2, 3), ")"
          )
        )
      }
    }

    recommendations <- c(
      recommendations,
      "",
      "[NEXT STEPS] After adding variables, re-run this analysis to check for remaining patterns.",
      "[TIP] Consider interaction terms if multiple variables are significant."
    )
  }

  return(recommendations)
}

#' @title Create Diagnostic Plots for Residual Patterns
#' @description Internal function to create plots showing residual patterns
#' @param residual_data Data frame with residuals and candidate variables
#' @param significant_vars Character vector of significant variables to plot
#' @param residual_type Type of residuals being analysed
#' @return List of ggplot objects
#' @noRd
create_residual_pattern_plots <- function(residual_data, significant_vars, residual_type) {
  plots <- list()

  max_plots <- min(length(significant_vars), 6)
  if (max_plots > 0) {
    for (i in seq_len(max_plots)) {
      var <- significant_vars[i]
      p <- ggplot2::ggplot(residual_data, ggplot2::aes(x = .data[[var]], y = .data$residuals)) +
        ggplot2::geom_point(alpha = 0.6, colour = "steelblue") +
        ggplot2::geom_smooth(method = "loess", se = TRUE, colour = "royalblue") +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
        ggplot2::labs(
          x = var,
          y = paste(tools::toTitleCase(residual_type), "Residuals"),
          title = paste("Residual Pattern:", var)
        )

      plots[[var]] <- p
    }
  }

  return(plots)
}

#' @title Print Method for Residual Pattern Analysis
#' @description Print method for objects of class 'residual_pattern_analysis'
#' @param x A residual_pattern_analysis object
#' @param ... Additional arguments (unused)
#' @export
print.residual_pattern_analysis <- function(x, ...) {
  cat("Residual Pattern Analysis\n")
  cat("========================\n\n")

  cat("Analysis Summary:\n")
  cat("  Residual type:", x$analysis_info$residual_type, "\n")
  cat("  Observations:", x$analysis_info$n_observations, "\n")
  cat("  Candidate variables tested:", x$analysis_info$n_candidate_vars, "\n")
  cat("  Significance level:", x$analysis_info$significance_level, "\n")
  cat("  Min R-squared threshold:", x$analysis_info$min_r_squared, "\n\n")

  if (nrow(x$linear_results) > 0) {
    cat("Top Linear Relationships (by R-squared):\n")
    print_top <- head(x$linear_results[, c("variable", "r_squared", "p_value", "significant")], 5)
    print(print_top, row.names = FALSE, digits = 4)
    cat("\n")
  }

  if (!is.null(x$smooth_results) && nrow(x$smooth_results) > 0) {
    cat("Top Smooth Relationships (by R-squared):\n")
    print_top <- head(x$smooth_results[, c("variable", "r_squared", "edf", "p_value", "significant")], 5)
    print(print_top, row.names = FALSE, digits = 4)
    cat("\n")
  }

  cat("Recommendations:\n")
  for (rec in x$recommendations) {
    cat(rec, "\n")
  }

  if (!is.null(x$plots)) {
    cat("\nDiagnostic plots generated for", length(x$plots), "significant variable(s).\n")
    cat("Access plots with: obj$plots$variable_name\n")
  }

  invisible(x)
}

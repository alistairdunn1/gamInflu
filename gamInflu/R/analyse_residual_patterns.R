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
#' @param max_missing_percent Numeric; maximum percentage of missing values allowed for candidate variables (default 20).
#' Variables with more missing data will be excluded from analysis.
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
                                      max_missing_percent = 20,
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
    stop("Invalid residual_type. Must be one of: deviance, pearson, response, working", call. = FALSE)
  )

  # Get variables already in the model
  model_vars <- all.vars(formula(model))
  response_var <- model_vars[1]
  predictor_vars <- model_vars[-1]

  # Determine candidate variables
  if (is.null(candidate_vars)) {
    # Use all numeric variables not already in the model
    numeric_vars <- names(data)[sapply(data, is.numeric)]
    candidate_vars <- setdiff(numeric_vars, c(response_var, predictor_vars))
  }

  # Exclude specified variables
  if (!is.null(exclude_vars)) {
    candidate_vars <- setdiff(candidate_vars, exclude_vars)
  }

  # Filter by missing data percentage
  missing_percent <- sapply(candidate_vars, function(var) {
    sum(is.na(data[[var]])) / nrow(data) * 100
  })

  candidate_vars <- candidate_vars[missing_percent <= max_missing_percent]

  if (length(candidate_vars) == 0) {
    stop("No suitable candidate variables found after filtering.", call. = FALSE)
  }

  # Create residual data frame
  residual_data <- data.frame(
    residuals = residuals_vec,
    data[candidate_vars],
    stringsAsFactors = FALSE
  )

  # Remove rows with missing data
  residual_data <- residual_data[complete.cases(residual_data), ]

  # Fit linear models
  linear_results <- data.frame(
    variable = character(),
    r_squared = numeric(),
    p_value = numeric(),
    significant = logical(),
    stringsAsFactors = FALSE
  )

  for (var in candidate_vars) {
    if (var %in% names(residual_data) && sum(!is.na(residual_data[[var]])) > 0) {
      tryCatch(
        {
          fit <- lm(residuals ~ get(var), data = residual_data)
          r_squared <- summary(fit)$r.squared
          p_value <- anova(fit)$"Pr(>F)"[1]

          linear_results <- rbind(linear_results, data.frame(
            variable = var,
            r_squared = r_squared,
            p_value = p_value,
            significant = p_value < significance_level && r_squared > min_r_squared,
            stringsAsFactors = FALSE
          ))
        },
        error = function(e) {
          # Skip variables that cause errors
        }
      )
    }
  }

  # Sort by significance and R-squared
  linear_results <- linear_results[order(-linear_results$significant, -linear_results$r_squared), ]

  # Fit smooth models if requested
  smooth_results <- NULL
  if (smooth_terms && nrow(linear_results) > 0) {
    smooth_results <- data.frame(
      variable = character(),
      r_squared = numeric(),
      p_value = numeric(),
      significant = logical(),
      stringsAsFactors = FALSE
    )

    for (var in candidate_vars) {
      if (var %in% names(residual_data) && sum(!is.na(residual_data[[var]])) > 0) {
        tryCatch(
          {
            fit <- gam(residuals ~ s(get(var)), data = residual_data)
            r_squared <- summary(fit)$r.sq
            p_value <- summary(fit)$s.table[1, 4] # p-value for smooth term

            smooth_results <- rbind(smooth_results, data.frame(
              variable = var,
              r_squared = r_squared,
              p_value = p_value,
              significant = p_value < significance_level && r_squared > min_r_squared,
              stringsAsFactors = FALSE
            ))
          },
          error = function(e) {
            # Skip variables that cause errors
          }
        )
      }
    }

    # Sort by significance and R-squared
    smooth_results <- smooth_results[order(-smooth_results$significant, -smooth_results$r_squared), ]
  }

  # Identify significant variables
  significant_vars <- unique(c(
    linear_results$variable[linear_results$significant],
    if (!is.null(smooth_results)) smooth_results$variable[smooth_results$significant] else character()
  ))

  # Generate recommendations
  recommendations <- character()
  if (length(significant_vars) > 0) {
    recommendations <- c(
      "Consider adding the following variables to your model:",
      paste("-", significant_vars),
      "",
      "These variables show significant patterns in the residuals,",
      "suggesting they may improve model fit if included."
    )
  } else {
    recommendations <- c(
      "No significant patterns found in residuals with tested variables.",
      "This suggests the current model captures the main relationships well."
    )
  }

  # Generate plots if requested
  plots <- list()
  if (plot_results && length(significant_vars) > 0) {
    n_plots <- min(6, length(significant_vars))
    for (var in significant_vars[seq_len(n_plots)]) { # Limit to 6 plots
      if (var %in% names(residual_data)) {
        p <- ggplot(residual_data, aes(x = get(var), y = residuals)) +
          geom_point(alpha = 0.6) +
          geom_smooth(method = "loess", se = TRUE, color = "red") +
          labs(
            x = var,
            y = "Residuals",
            title = paste("Residuals vs", var)
          ) +
          theme_minimal()

        plots[[var]] <- p
      }
    }
  }

  # Create result object
  result <- structure(
    list(
      linear_results = linear_results,
      smooth_results = smooth_results,
      significant_vars = significant_vars,
      recommendations = recommendations,
      residual_data = residual_data,
      plots = plots,
      call = match.call()
    ),
    class = "residual_pattern_analysis"
  )

  return(result)
}

#' @title Print Method for Residual Pattern Analysis
#' @description Print method for objects of class `residual_pattern_analysis`.
#' @param x A `residual_pattern_analysis` object.
#' @param ... Additional arguments passed to print methods.
#' @export
print.residual_pattern_analysis <- function(x, ...) {
  cat("Residual Pattern Analysis\n")
  cat("========================\n\n")

  if (nrow(x$linear_results) > 0) {
    cat("Linear Model Results (top 10):\n")
    cat("------------------------------\n")
    top_linear <- head(x$linear_results, 10)
    print(top_linear, row.names = FALSE)
    cat("\n")
  }

  if (!is.null(x$smooth_results) && nrow(x$smooth_results) > 0) {
    cat("Smooth Model Results (top 10):\n")
    cat("------------------------------\n")
    top_smooth <- head(x$smooth_results, 10)
    print(top_smooth, row.names = FALSE)
    cat("\n")
  }

  if (length(x$significant_vars) > 0) {
    cat("Significant Variables:\n")
    cat("---------------------\n")
    cat(paste(x$significant_vars, collapse = ", "), "\n\n")
  }

  cat("Recommendations:\n")
  cat("---------------\n")
  cat(paste(x$recommendations, collapse = "\n"), "\n")

  if (length(x$plots) > 0) {
    cat("\nDiagnostic plots available for:", paste(names(x$plots), collapse = ", "), "\n")
    cat("Access plots with: result$plots$variable_name\n")
  }

  invisible(x)
}

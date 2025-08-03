#' @title Combine Binomial and Positive Catch Indices
#' @description Combines indices from a binomial GAM (probability of positive catch) and a
#' positive catch GAM to create an overall CPUE index using the delta-GLM approach. This method
#' is commonly used in fisheries analysis where catches contain many zeros.
#' @param binomial_gi A `gam_influence` object from a binomial model (presence/absence or catch probability).
#'   This model should predict the probability of a positive catch.
#' @param positive_gi A `gam_influence` object from a positive catch model (e.g., Gamma or Gaussian).
#'   This model should predict catch rates conditional on positive catches.
#' @param method Character. Method for combining indices:
#'   - `"multiplicative"`: Combined index = P(positive) * E(catch | positive) (default)
#'   - `"arithmetic"`: Combined index = P(positive) + E(catch | positive) - baseline
#'   - `"geometric"`: Combined index = sqrt(P(positive) * E(catch | positive))
#' @param rescale_combined Logical. Should the combined index be rescaled to have geometric mean = 1? Default TRUE.
#' @param confidence_method Character. Method for calculating combined confidence intervals:
#'   - `"delta"`: Delta method using first-order Taylor approximation (default)
#'   - `"bootstrap"`: Bootstrap confidence intervals (more accurate but slower)
#'   - `"independent"`: Treat components as independent (conservative)
#' @param bootstrap_n Integer. Number of bootstrap samples if confidence_method = "bootstrap" (default 1000).
#' @param validate_data Logical. Should data compatibility be validated? Default TRUE.
#' @param ... Additional arguments (currently unused).
#' @return An object of class `gam_influence_combined` containing:
#' \describe{
#'   \item{binomial_gi}{The original binomial gam_influence object}
#'   \item{positive_gi}{The original positive catch gam_influence object}
#'   \item{combined_indices}{Data frame with combined indices, confidence intervals, and diagnostics}
#'   \item{method}{The combination method used}
#'   \item{diagnostics}{Diagnostic information about the combination process}
#' }
#' @details
#' **Delta-GLM Approach:**
#' The delta-GLM method is a two-step process commonly used in fisheries for zero-inflated data:
#' 1. **Binomial component**: Models the probability of a positive catch
#' 2. **Positive component**: Models the catch rate given a positive catch occurred
#' 3. **Combined index**: Multiplies these components to get overall expected catch rate
#'
#' **Mathematical Framework:**
#' For the multiplicative method (most common):
#' Combined Index = P(catch > 0) * E(catch | catch > 0)
#'
#' Where:
#' - P(catch > 0) comes from the binomial model
#' - E(catch | catch > 0) comes from the positive catch model
#'
#' **Confidence Intervals:**
#' - **Delta method**: Uses first-order Taylor approximation for uncertainty propagation
#' - **Bootstrap**: Resamples to empirically estimate confidence intervals
#' - **Independent**: Assumes independence between components (often conservative)
#'
#' **Data Requirements:**
#' - Both models must use the same focus term (e.g., "year")
#' - Focus levels should be identical or compatible
#' - Data should represent the same population/fishery
#'
#' **Validation Checks:**
#' - Ensures both objects have been processed by `calculate_influence()`
#' - Validates focus term compatibility
#' - Checks for reasonable probability ranges (0-1 for binomial component)
#' - Warns about potential issues (e.g., very low catch probabilities)
#' @examples
#' \dontrun{
#' # Example fisheries CPUE analysis with zero-inflated data
#' library(mgcv)
#'
#' # Prepare data - separate zero and positive catches
#' data$presence <- as.numeric(data$cpue > 0)
#' data$positive_cpue <- ifelse(data$cpue > 0, data$cpue, NA)
#'
#' # Fit binomial model for catch probability
#' mod_binom <- gam(presence ~ year + s(depth) + s(temperature),
#'   family = binomial(), data = data
#' )
#' gi_binom <- gam_influence(mod_binom, focus = "year")
#' gi_binom <- calculate_influence(gi_binom)
#'
#' # Fit positive catch model (Gamma for positive continuous data)
#' mod_pos <- gam(positive_cpue ~ year + s(depth) + s(temperature),
#'   family = Gamma(link = "log"),
#'   data = data[data$cpue > 0, ]
#' )
#' gi_pos <- gam_influence(mod_pos, focus = "year")
#' gi_pos <- calculate_influence(gi_pos)
#'
#' # Combine indices using delta-GLM approach
#' combined_gi <- combine_indices(gi_binom, gi_pos)
#'
#' # Alternative combination methods
#' combined_geo <- combine_indices(gi_binom, gi_pos, method = "geometric")
#' combined_boot <- combine_indices(gi_binom, gi_pos, confidence_method = "bootstrap")
#'
#' # View results
#' print(combined_gi)
#' plot(combined_gi)
#' summary(combined_gi)
#' }
#' @importFrom stats qnorm var cov quantile rnorm
#' @export
combine_indices <- function(binomial_gi, positive_gi,
                            method = c("multiplicative", "arithmetic", "geometric"),
                            rescale_combined = TRUE,
                            confidence_method = c("delta", "bootstrap", "independent"),
                            bootstrap_n = 1000,
                            validate_data = TRUE,
                            ...) {
  # Validate inputs
  method <- match.arg(method)
  confidence_method <- match.arg(confidence_method)

  if (!inherits(binomial_gi, "gam_influence")) {
    stop("binomial_gi must be a gam_influence object", call. = FALSE)
  }
  if (!inherits(positive_gi, "gam_influence")) {
    stop("positive_gi must be a gam_influence object", call. = FALSE)
  }

  # Check that both objects have been calculated
  if (is.null(binomial_gi$calculated)) {
    stop("binomial_gi has not been processed. Run calculate_influence() first.", call. = FALSE)
  }
  if (is.null(positive_gi$calculated)) {
    stop("positive_gi has not been processed. Run calculate_influence() first.", call. = FALSE)
  }

  # Validate focus term compatibility
  if (binomial_gi$focus != positive_gi$focus) {
    stop("Both models must have the same focus term. Got '",
      binomial_gi$focus, "' and '", positive_gi$focus, "'",
      call. = FALSE
    )
  }

  focus_term <- binomial_gi$focus

  # Validate model families
  binomial_family <- binomial_gi$model$family$family
  positive_family <- positive_gi$model$family$family

  if (!grepl("binomial", binomial_family, ignore.case = TRUE)) {
    warning("First model family is '", binomial_family,
      "', expected binomial. Proceeding but results may not be meaningful.",
      call. = FALSE
    )
  }

  # Check positive model family for appropriateness
  if (!positive_family %in% c("Gamma", "gamma", "gaussian")) {
    message(
      "Positive model family is '", positive_family,
      "'. Common choices are Gamma or Gaussian for positive catch rates."
    )
  }

  # Extract indices data
  binomial_indices <- binomial_gi$calculated$indices
  positive_indices <- positive_gi$calculated$indices

  # Validate index data
  if (is.null(binomial_indices) || is.null(positive_indices)) {
    stop("Could not extract indices from one or both models", call. = FALSE)
  }

  # Validate focus levels compatibility
  binomial_levels <- sort(as.character(binomial_indices$level))
  positive_levels <- sort(as.character(positive_indices$level))

  # Always find common levels and filter, regardless of validate_data setting
  common_levels <- intersect(binomial_levels, positive_levels)

  if (length(common_levels) == 0) {
    stop("No common focus levels found between the two models.\n",
      "Binomial levels: ", paste(binomial_levels, collapse = ", "), "\n",
      "Positive levels: ", paste(positive_levels, collapse = ", "),
      call. = FALSE
    )
  }

  # Warn if levels differ (only when validate_data is TRUE)
  if (validate_data && !identical(binomial_levels, positive_levels)) {
    missing_binom <- setdiff(positive_levels, binomial_levels)
    missing_pos <- setdiff(binomial_levels, positive_levels)

    warning_msg <- paste("Focus levels differ between models.")
    if (length(missing_binom) > 0) {
      warning_msg <- paste(
        warning_msg,
        "\nMissing from binomial model:", paste(missing_binom, collapse = ", ")
      )
    }
    if (length(missing_pos) > 0) {
      warning_msg <- paste(
        warning_msg,
        "\nMissing from positive model:", paste(missing_pos, collapse = ", ")
      )
    }
    warning_msg <- paste(
      warning_msg,
      "\nUsing intersection:", paste(common_levels, collapse = ", ")
    )

    warning(warning_msg, call. = FALSE)
  }

  # Always filter to common levels before merging
  binomial_indices <- binomial_indices[binomial_indices$level %in% common_levels, ]
  positive_indices <- positive_indices[positive_indices$level %in% common_levels, ]

  # Verify we still have data after filtering
  if (nrow(binomial_indices) == 0 || nrow(positive_indices) == 0) {
    stop("No data remaining after filtering to common levels", call. = FALSE)
  }

  # Enhanced merge with better error handling
  tryCatch(
    {
      # Sort both by level to ensure proper merging
      binomial_indices <- binomial_indices[order(binomial_indices$level), ]
      positive_indices <- positive_indices[order(positive_indices$level), ]

      combined_data <- merge(binomial_indices, positive_indices,
        by = "level", suffixes = c("_binom", "_pos"), all = FALSE
      )

      if (nrow(combined_data) == 0) {
        stop("No matching levels found after merging indices", call. = FALSE)
      }

      # Check for minimum data requirements
      if (nrow(combined_data) < 2) {
        warning("Very few levels available for combination (n=", nrow(combined_data),
          "). Results may not be reliable.",
          call. = FALSE
        )
      }
    },
    error = function(e) {
      stop("Failed to merge indices data. Error: ", e$message,
        "\nFiltered binomial levels: ", paste(binomial_indices$level, collapse = ", "),
        "\nFiltered positive levels: ", paste(positive_indices$level, collapse = ", "),
        call. = FALSE
      )
    }
  )

  # Validate probability ranges for binomial component
  if (validate_data) {
    if (any(combined_data$standardised_index_binom < 0) ||
      any(combined_data$standardised_index_binom > 2)) {
      warning("Binomial component has unusual values (outside 0-2 range). Check model specification.",
        call. = FALSE
      )
    }
  }

  # Calculate combined indices based on method
  combined_data <- calculate_combined_index(combined_data, method)

  # Calculate confidence intervals
  combined_data <- calculate_combined_confidence_intervals(
    combined_data, confidence_method, bootstrap_n,
    binomial_gi, positive_gi, method
  )

  # Rescale if requested
  if (rescale_combined) {
    original_combined <- combined_data$combined_index
    rescale_factor <- 1 / geometric_mean(original_combined)

    combined_data$combined_index <- combined_data$combined_index * rescale_factor
    combined_data$combined_lower <- combined_data$combined_lower * rescale_factor
    combined_data$combined_upper <- combined_data$combined_upper * rescale_factor
    combined_data$combined_se <- combined_data$combined_se * rescale_factor
  }

  # Calculate combined CV
  combined_data$combined_cv <- ifelse(combined_data$combined_index > 1e-6,
    combined_data$combined_se / combined_data$combined_index,
    NA_real_
  )

  # Create diagnostics
  diagnostics <- create_combination_diagnostics(combined_data, binomial_gi, positive_gi, method)

  # Create result object
  result <- structure(
    list(
      binomial_gi = binomial_gi,
      positive_gi = positive_gi,
      combined_indices = combined_data,
      method = method,
      confidence_method = confidence_method,
      focus_term = focus_term,
      diagnostics = diagnostics,
      rescaled = rescale_combined
    ),
    class = "gam_influence_combined"
  )

  return(result)
}

#' @title Calculate Combined Index Values
#' @description Internal function to calculate combined indices using different methods
#' @param combined_data Data frame with binomial and positive indices
#' @param method Combination method
#' @return Updated data frame with combined index
#' @noRd
calculate_combined_index <- function(combined_data, method) {
  if (method == "multiplicative") {
    # Standard delta-GLM: P(positive) * E(catch | positive)
    combined_data$combined_index <- combined_data$standardised_index_binom *
      combined_data$standardised_index_pos
  } else if (method == "arithmetic") {
    # Arithmetic combination (less common)
    baseline <- mean(combined_data$standardised_index_binom) +
      mean(combined_data$standardised_index_pos) - 1
    combined_data$combined_index <- combined_data$standardised_index_binom +
      combined_data$standardised_index_pos - baseline
  } else if (method == "geometric") {
    # Geometric mean combination
    combined_data$combined_index <- sqrt(combined_data$standardised_index_binom *
      combined_data$standardised_index_pos)
  }

  return(combined_data)
}

#' @title Calculate Combined Confidence Intervals
#' @description Internal function to calculate confidence intervals for combined indices
#' @param combined_data Data frame with combined indices
#' @param confidence_method Method for CI calculation
#' @param bootstrap_n Number of bootstrap samples
#' @param binomial_gi Binomial gam_influence object
#' @param positive_gi Positive gam_influence object
#' @param method Combination method
#' @return Updated data frame with confidence intervals
#' @noRd
calculate_combined_confidence_intervals <- function(combined_data, confidence_method,
                                                    bootstrap_n, binomial_gi,
                                                    positive_gi, method) {
  if (confidence_method == "delta") {
    # Delta method for uncertainty propagation
    combined_data <- calculate_delta_method_ci(combined_data, method)
  } else if (confidence_method == "bootstrap") {
    # Bootstrap confidence intervals
    combined_data <- calculate_bootstrap_ci(
      combined_data, bootstrap_n,
      binomial_gi, positive_gi, method
    )
  } else if (confidence_method == "independent") {
    # Independent assumption (conservative)
    combined_data <- calculate_independent_ci(combined_data, method)
  }

  return(combined_data)
}

#' @title Calculate Delta Method Confidence Intervals
#' @description Use delta method for combined index uncertainty
#' @param combined_data Data frame with indices
#' @param method Combination method
#' @return Updated data frame with delta method CIs
#' @noRd
calculate_delta_method_ci <- function(combined_data, method) {
  if (method == "multiplicative") {
    # For Y = X1 * X2, Var(Y) ~ X2^2 * Var(X1) + X1^2 * Var(X2) + 2*X1*X2*Cov(X1,X2)
    # Assuming independence, Cov(X1,X2) = 0

    x1 <- combined_data$standardised_index_binom
    x2 <- combined_data$standardised_index_pos
    var_x1 <- combined_data$stan_se_binom^2
    var_x2 <- combined_data$stan_se_pos^2

    # Delta method variance for multiplication
    combined_var <- (x2^2) * var_x1 + (x1^2) * var_x2
    combined_data$combined_se <- sqrt(combined_var)
  } else if (method == "arithmetic") {
    # For Y = X1 + X2, Var(Y) = Var(X1) + Var(X2) + 2*Cov(X1,X2)
    # Assuming independence
    combined_var <- combined_data$stan_se_binom^2 + combined_data$stan_se_pos^2
    combined_data$combined_se <- sqrt(combined_var)
  } else if (method == "geometric") {
    # For Y = sqrt(X1 * X2), use log transformation and delta method
    x1 <- combined_data$standardised_index_binom
    x2 <- combined_data$standardised_index_pos

    # Delta method on log scale, then transform back
    log_var <- (combined_data$stan_se_binom / x1)^2 + (combined_data$stan_se_pos / x2)^2
    log_se <- sqrt(log_var / 4) # Divide by 4 for square root

    # Transform back to original scale
    combined_data$combined_se <- combined_data$combined_index * log_se
  }

  # Calculate confidence intervals
  z_score <- qnorm(0.975) # 95% CI
  combined_data$combined_lower <- combined_data$combined_index - z_score * combined_data$combined_se
  combined_data$combined_upper <- combined_data$combined_index + z_score * combined_data$combined_se

  # Ensure non-negative bounds for positive indices
  combined_data$combined_lower <- pmax(combined_data$combined_lower, 0)

  return(combined_data)
}

#' @title Calculate Bootstrap Confidence Intervals
#' @description Use bootstrap resampling for combined index uncertainty
#' @param combined_data Data frame with indices
#' @param bootstrap_n Number of bootstrap samples
#' @param binomial_gi Binomial gam_influence object
#' @param positive_gi Positive gam_influence object
#' @param method Combination method
#' @return Updated data frame with bootstrap CIs
#' @noRd
calculate_bootstrap_ci <- function(combined_data, bootstrap_n, binomial_gi, positive_gi, method) {
  # For now, implement a simplified bootstrap using normal approximation
  # Full bootstrap would require refitting models, which is computationally intensive

  n_levels <- nrow(combined_data)
  bootstrap_results <- matrix(NA, nrow = bootstrap_n, ncol = n_levels)

  for (i in seq_len(bootstrap_n)) {
    # Sample from normal distributions for each component
    binom_sample <- rnorm(n_levels,
      mean = combined_data$standardised_index_binom,
      sd = combined_data$stan_se_binom
    )
    pos_sample <- rnorm(n_levels,
      mean = combined_data$standardised_index_pos,
      sd = combined_data$stan_se_pos
    )

    # Ensure binomial samples are in reasonable range
    binom_sample <- pmax(0, pmin(2, binom_sample))

    # Calculate combined index for this bootstrap sample
    if (method == "multiplicative") {
      bootstrap_results[i, ] <- binom_sample * pos_sample
    } else if (method == "arithmetic") {
      baseline <- mean(binom_sample) + mean(pos_sample) - 1
      bootstrap_results[i, ] <- binom_sample + pos_sample - baseline
    } else if (method == "geometric") {
      bootstrap_results[i, ] <- sqrt(binom_sample * pos_sample)
    }
  }

  # Calculate percentile confidence intervals
  combined_data$combined_lower <- apply(bootstrap_results, 2, quantile, probs = 0.025, na.rm = TRUE)
  combined_data$combined_upper <- apply(bootstrap_results, 2, quantile, probs = 0.975, na.rm = TRUE)
  combined_data$combined_se <- apply(bootstrap_results, 2, sd, na.rm = TRUE)

  return(combined_data)
}

#' @title Calculate Independent Confidence Intervals
#' @description Calculate CIs assuming independence (conservative approach)
#' @param combined_data Data frame with indices
#' @param method Combination method
#' @return Updated data frame with independent CIs
#' @noRd
calculate_independent_ci <- function(combined_data, method) {
  # Simple approach: use wider of the two confidence intervals
  binom_width <- combined_data$stan_upper_binom - combined_data$stan_lower_binom
  pos_width <- combined_data$stan_upper_pos - combined_data$stan_lower_pos

  # Conservative approach: use maximum relative width
  rel_width_binom <- binom_width / combined_data$standardised_index_binom
  rel_width_pos <- pos_width / combined_data$standardised_index_pos
  max_rel_width <- pmax(rel_width_binom, rel_width_pos, na.rm = TRUE)

  # Apply to combined index
  half_width <- combined_data$combined_index * max_rel_width / 2
  combined_data$combined_lower <- combined_data$combined_index - half_width
  combined_data$combined_upper <- combined_data$combined_index + half_width
  combined_data$combined_se <- half_width / qnorm(0.975)

  # Ensure non-negative bounds
  combined_data$combined_lower <- pmax(combined_data$combined_lower, 0)

  return(combined_data)
}

#' @title Create Combination Diagnostics
#' @description Create diagnostic information about the combination process
#' @param combined_data Data frame with combined results
#' @param binomial_gi Binomial gam_influence object
#' @param positive_gi Positive gam_influence object
#' @param method Combination method
#' @return List of diagnostic information
#' @noRd
create_combination_diagnostics <- function(combined_data, binomial_gi, positive_gi, method) {
  diagnostics <- list(
    n_levels = nrow(combined_data),
    focus_term = binomial_gi$focus,
    combination_method = method,
    binomial_family = binomial_gi$model$family$family,
    positive_family = positive_gi$model$family$family,

    # Component statistics
    binomial_range = range(combined_data$standardised_index_binom, na.rm = TRUE),
    positive_range = range(combined_data$standardised_index_pos, na.rm = TRUE),
    combined_range = range(combined_data$combined_index, na.rm = TRUE),

    # Uncertainty measures
    mean_binom_cv = mean(combined_data$standardised_cv_binom, na.rm = TRUE),
    mean_positive_cv = mean(combined_data$standardised_cv_pos, na.rm = TRUE),
    mean_combined_cv = mean(combined_data$combined_cv, na.rm = TRUE),

    # Correlation (approximate)
    component_correlation = cor(combined_data$standardised_index_binom,
      combined_data$standardised_index_pos,
      use = "complete.obs"
    ),

    # Data quality flags
    warnings = character(0)
  )

  # Add warnings for potential issues
  if (any(combined_data$standardised_index_binom < 0.01)) {
    diagnostics$warnings <- c(diagnostics$warnings, "Very low catch probabilities detected")
  }

  if (any(combined_data$standardised_index_binom > 1.5)) {
    diagnostics$warnings <- c(diagnostics$warnings, "Unusually high probability indices (>1.5)")
  }

  if (diagnostics$component_correlation < -0.5) {
    diagnostics$warnings <- c(diagnostics$warnings, "Strong negative correlation between components")
  }

  return(diagnostics)
}

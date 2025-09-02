#' @title Perform Influence Calculations
#' @description Perform calculations on a \code{gam_influence} object, including unstandardised and standardised indices, step-wise model building, and influence statistics. Supports subset-based analysis and multiple GLM families.
#' @param obj A \code{gam_influence} object.
#' @param islog Logical. Is the response variable log-transformed? If NULL (default), inferred from response name or model family.
#' @param rescale_method Character. How to rescale the indices. Options: "auto" (default), "geometric_mean", "arithmetic_mean", "raw", "custom".
#' \describe{
#'   \item{auto}{Uses sensible defaults based on model family and response transformation:
#'     \itemize{
#'       \item Binomial: "raw" (preserves probability scale 0-1)
#'       \item Gaussian with log-transformed response: "geometric_mean"
#'       \item Gaussian with original scale response: "arithmetic_mean"
#'       \item Gamma/Poisson: "geometric_mean" (appropriate for positive/count data)
#'     }
#'   }
#'   \item{geometric_mean}{Rescale relative to geometric mean (good for log-normal data)}
#'   \item{arithmetic_mean}{Rescale relative to arithmetic mean (good for normal data)}
#'   \item{raw}{No rescaling (preserves original scale)}
#'   \item{custom}{Rescale to custom value specified by custom_rescale_value}
#' }
#' @param custom_rescale_value Numeric. Custom rescaling value when rescale_method = "custom".
#' @param confidence_level Numeric. Confidence level for standardised index intervals (default 0.95).
#' @param family Character. How to handle different GLM families. Options: "auto" (default), "gaussian", "binomial", "gamma", "poisson".
#' @param subset_var Character. Name of variable to subset on (e.g., "area", "gear_type").
#' @param subset_value Value to subset by (e.g., "North", "Trawl"). If both subset_var and subset_value are provided, analysis is performed only on the subset, using the full model for predictions.
#' @param separate_by_smooth Logical. Whether to separate by-variable smooths (e.g., s(depth, by=gear)) into individual components for influence analysis. If FALSE (default), by-smooths are treated as separate terms for each level. If TRUE, they are grouped together as a single term. This affects how smooth terms with by-variables are handled in influence calculations.
#' @param ... Additional arguments (currently unused).
#' @return The \code{gam_influence} object, now containing a \code{calculated} list with data frames for indices, summary stats, influences, predictions, and s.e. of predictions.
#' @details
#' The function calculates unstandardised and standardised indices, confidence intervals, standard errors, and coefficients of variation (CV) using model predictions. For log-transformed data, calculations are performed in log space and exponentiated. Step-wise model building and influence statistics are also included.
#' @examples
#' \dontrun{
#' data$year <- factor(data$year)
#' mod <- mgcv::gam(cpue ~ year + s(depth), data = data, family = Gamma(link = "log"))
#' gi <- gam_influence(mod, focus = "year")
#' gi <- calculate_influence(gi)
#' }
#' @export
calculate_influence <- function(obj, islog = NULL,
                                rescale_method = c("auto", "geometric_mean", "arithmetic_mean", "raw", "custom"),
                                custom_rescale_value = 1,
                                confidence_level = 0.95,
                                family = c("auto", "gaussian", "binomial", "gamma", "poisson", "weibull"),
                                subset_var = NULL,
                                subset_value = NULL,
                                separate_by_smooth = FALSE,
                                ...) {
  UseMethod("calculate_influence")
}

#' @rdname calculate_influence
#' @method calculate_influence gam_influence
#' @importFrom stats predict aggregate logLik AIC deviance terms update as.formula var cov setNames qnorm
#' @export
calculate_influence.gam_influence <- function(obj, islog = NULL,
                                              rescale_method = c("auto", "geometric_mean", "arithmetic_mean", "raw", "custom"),
                                              custom_rescale_value = 1,
                                              confidence_level = 0.95,
                                              family = c("auto", "gaussian", "binomial", "gamma", "poisson", "weibull"),
                                              subset_var = NULL,
                                              subset_value = NULL,
                                              separate_by_smooth = FALSE,
                                              ...) {
  # Validate inputs
  validate_gam_influence(obj)
  if (confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must be between 0 and 1", call. = FALSE)
  }

  # Extract method preference from object
  use_coeff_method <- obj$use_coeff_method
  if (is.null(use_coeff_method)) {
    use_coeff_method <- TRUE # Default to coefficient method for backwards compatibility
  }
  
  # Check for Gamma with inverse link - override method if needed
  if (is.null(obj$model$family)) {
    # Not a GLM/GAM, so cannot check family
    is_gamma_inverse <- FALSE
  } else {
    is_gamma_inverse <- obj$model$family$family == "Gamma" && 
                        obj$model$family$link == "inverse"
                        
    if (is_gamma_inverse && use_coeff_method) {
      message("Detected Gamma model with inverse link - forcing prediction-based method")
      message("(Coefficient-based method can produce NaN values with inverse link)")
      use_coeff_method <- FALSE
    }
  }

  # --- Subset Analysis Implementation ---
  # Handle subset-based influence analysis if requested
  if (!is.null(subset_var) && !is.null(subset_value)) {
    # Validate subset parameters
    if (!subset_var %in% names(obj$data)) {
      stop("subset_var '", subset_var, "' not found in data", call. = FALSE)
    }

    # Create subset data
    subset_indices <- obj$data[[subset_var]] == subset_value
    if (sum(subset_indices) == 0) {
      stop("No observations found for ", subset_var, " = '", subset_value, "'", call. = FALSE)
    }

    # Update the obj$data to subset for influence calculations
    # Keep original for model context
    original_data <- obj$data
    obj$data <- obj$data[subset_indices, , drop = FALSE]

    message(
      "Subset analysis: Using ", nrow(obj$data), " observations where ",
      subset_var, " = '", subset_value, "' (from ", nrow(original_data), " total)"
    )

    # For subset analysis, we don't need multiple focus levels since we use the full model
    # The focus term validation only applies to full analysis
  }

  # --- Setup and Family Detection ---
  family <- match.arg(family)
  rescale_method <- match.arg(rescale_method)

  # Detect model family
  if (family == "auto") {
    model_family <- obj$model$family$family
    model_link <- obj$model$family$link

    # Map common families, including pattern matching for Tweedie and NB
    if (grepl("^Tweedie\\(", model_family)) {
      family_detected <- "gamma" # Treat Tweedie like Gamma for calculations
    } else if (grepl("^Negative Binomial\\(", model_family)) {
      family_detected <- "poisson" # Treat NB like Poisson for calculations
    } else {
      # Check for Weibull family (custom attribute)
      if (!is.null(attr(obj$model$family, "weibull_approximation")) &&
        attr(obj$model$family, "weibull_approximation") == TRUE) {
        family_detected <- "weibull"
      } else {
        family_detected <- switch(model_family,
          "gaussian" = "gaussian",
          "binomial" = "binomial",
          "Gamma" = "gamma",
          "gamma" = "gamma",
          "poisson" = "poisson",
          "weibull" = "weibull",
          "quasipoisson" = "poisson", # Treat quasi-poisson like poisson
          "quasibinomial" = "binomial", # Treat quasi-binomial like binomial
          "quasi" = "gaussian", # Default for other quasi families
          "gaussian" # Default fallback
        )
      }
    }
    family <- family_detected
    message("Detected model family: ", model_family, " with link: ", model_link)
    message("Using family: ", family)
    
    # Special handling for Gamma with inverse link - force prediction method
    # Note: We already check for this at the beginning of the function,
    # but keep this check for safety with user-specified "auto" family parameter
  }

  # Enhanced islog detection - ONLY based on response name, NOT family link
  if (is.null(islog)) {
    if (is.null(obj$islog)) {
      # Check response name for log transformation (user manually logged response)
      # Look for "log(" or variable names starting with "log_" or "ln_"
      response_log <- substr(obj$response, 1, 4) == "log(" ||
        substr(obj$response, 1, 4) == "log_" ||
        substr(obj$response, 1, 3) == "ln_"

      # NOTE: We do NOT automatically set islog=TRUE for log-link families
      # Log link (e.g., Gamma(link="log")) != Pre-logged response (log(y) ~ ...)
      # Log link: response on original scale, link transforms internally
      # Pre-logged: response on log scale, user transformed before fitting

      islog <- response_log
      if (response_log) {
        message("Detected pre-logged response '", obj$response, "', setting islog = TRUE")
      } else {
        message("Response '", obj$response, "' appears to be on original scale, setting islog = FALSE")
      }
    } else {
      islog <- obj$islog
    }
  } else {
    # User explicitly provided islog parameter
    obj$islog <- islog
  }

  # Auto rescaling method selection based on family and islog
  if (rescale_method == "auto") {
    rescale_method <- switch(family,
      "binomial" = "raw", # Preserve probability scale for binomial
      "gaussian" = if (islog) "geometric_mean" else "arithmetic_mean", # Consider log transformation
      "gamma" = "geometric_mean", # Always geometric for gamma (positive data)
      "poisson" = "geometric_mean", # Geometric for count data
      "weibull" = "geometric_mean", # Geometric for Weibull (positive data)
      "geometric_mean" # Default fallback
    )
    message("Auto-selected rescaling method: ", rescale_method, " (based on family: ", family, ", islog: ", islog, ")")
  }

  # Always ensure obj$islog is set for plotting functions
  obj$islog <- islog

  observed <- obj$data[[obj$response]]

  # Special handling for binomial models with raw rescaling
  # For binomial, raw rescaling should preserve the probability scale (0-1)
  is_binomial <- family == "binomial"
  preserve_probability_scale <- is_binomial && rescale_method == "raw"

  # For binomial models, always preserve probability scale for unstandardised indices
  # but warn if non-raw rescaling is requested
  preserve_unstandardised_probability_scale <- is_binomial
  if (is_binomial && rescale_method != "raw") {
    message("Note: For binomial models, unstandardised indices will preserve probability scale (0-1) regardless of rescaling method.")
    message("Rescaling method '", rescale_method, "' will be applied only to standardised indices.")
  }

  # --- 1. Unstandardised Index ---
  # Calculate the raw, unadjusted index for the focus term.
  # Uses family-appropriate methods for different GLM families.
  indices_df <- calculate_unstandardised_index(observed, obj$data[[obj$focus]], islog, family, preserve_unstandardised_probability_scale)

  # --- 2. Standardised Index (from the full model) ---
  # This is the final index after accounting for all terms in the model.
  # We use type="response" to get full model predictions with proper uncertainty
  # For subset analysis, we predict using the full model but only on subset data
  preds_full_response <- predict(obj$model, newdata = obj$data, type = "response", se.fit = TRUE)

  # Also get terms predictions to extract the focus term's partial effect for influence calculations
  preds_full_terms <- predict(obj$model, newdata = obj$data, type = "terms", se.fit = TRUE)

  # Create data frame with response predictions and focus term levels
  response_pred_df <- data.frame(
    level = obj$data[[obj$focus]],
    pred = preds_full_response$fit,
    se = preds_full_response$se.fit
  )

  # Aggregate by focus term level to get mean predictions and standard errors
  # IMPORTANT: For standard errors, we need SE of the mean, not mean of SEs
  stan_df <- aggregate(cbind(pred, se) ~ level,
    data = response_pred_df,
    FUN = function(x) {
      if (length(x) == 1) {
        return(x)
      }
      c(mean = mean(x), se_of_mean = sqrt(mean(x^2)))
    }
  )

  # Handle the aggregation result format
  if (is.matrix(stan_df$pred)) {
    # Multiple observations per level - use corrected SE calculation
    pred_means <- stan_df$pred[, "mean"]
    se_corrected <- stan_df$se[, "se_of_mean"]
    stan_df <- data.frame(
      level = stan_df$level,
      pred = pred_means,
      se = se_corrected
    )
  } else {
    # Single observation per level - keep original
    stan_df <- data.frame(
      level = stan_df$level,
      pred = stan_df$pred,
      se = stan_df$se
    )
  }

  # Convert to relative index (standardised) using either coefficient or prediction method
  base_pred <- mean(stan_df$pred)

  # Choose calculation method
  if (use_coeff_method) {
    # --- COEFFICIENT-BASED APPROACH (traditional method) ---
    message("Using coefficient-based CI calculation")

    # Extract coefficients for the focus term
    all_coeffs <- coef(obj$model)
    focus_term_pattern <- paste0("^", obj$focus)
    focus_coeff_indices <- grep(focus_term_pattern, names(all_coeffs))

    if (length(focus_coeff_indices) == 0) {
      stop("Could not find coefficients for focus term '", obj$focus, "'", call. = FALSE)
    }

    # Get coefficients for non-reference levels only
    focus_coeffs_nonref <- all_coeffs[focus_coeff_indices]

    # Create full coefficient vector including reference level as 0
    # For subset analysis, use all levels from the original model
    if (!is.null(subset_var) && !is.null(subset_value) && exists("original_data", inherits = FALSE)) {
      # Use levels from original data for coefficient calculations
      focus_levels <- levels(original_data[[obj$focus]])
    } else {
      # Normal analysis - use current data levels
      focus_levels <- levels(obj$data[[obj$focus]])
    }
    n_levels <- length(focus_levels)
    focus_coeffs <- numeric(n_levels)
    names(focus_coeffs) <- focus_levels

    # Reference level (first level) = 0, others from model
    focus_coeffs[1] <- 0 # Reference level
    if (length(focus_coeffs_nonref) > 0) {
      # Extract level names from coefficient names
      coeff_level_names <- gsub(focus_term_pattern, "", names(focus_coeffs_nonref))
      focus_coeffs[coeff_level_names] <- focus_coeffs_nonref
    }

    # Calculate standard errors using Francis method
    # The Francis method calculates SEs for relative differences from the mean
    # This ensures all levels get non-zero confidence intervals
    if (length(focus_coeff_indices) > 0) {
      # Get the variance-covariance matrix including intercept
      all_vcov <- vcov(obj$model)

      # Need intercept + focus term coefficients for Francis method
      intercept_name <- "(Intercept)"
      if (intercept_name %in% rownames(all_vcov)) {
        vcov_names <- c(intercept_name, names(all_coeffs)[focus_coeff_indices])
        focus_vcov_full <- all_vcov[vcov_names, vcov_names, drop = FALSE]

        # Francis method transformation matrix Q
        # Q transforms from absolute to relative-to-mean effects
        n_levels <- length(focus_coeffs)
        Q <- matrix(-1 / n_levels, nrow = n_levels, ncol = ncol(focus_vcov_full))

        # First row: reference level = intercept - mean(all levels)
        Q[1, 1] <- 1 - 1 / n_levels # intercept coefficient
        Q[1, -1] <- -1 / n_levels # other coefficients

        # Subsequent rows: (intercept + coeff_i) - mean(all levels)
        for (i in 2:n_levels) {
          Q[i, 1] <- 1 - 1 / n_levels # intercept coefficient
          Q[i, i] <- 1 - 1 / n_levels # own coefficient
          Q[i, -c(1, i)] <- -1 / n_levels # other coefficients
        }

        # Calculate variance matrix for relative effects: Q * V * Q'
        V_relative <- Q %*% focus_vcov_full %*% t(Q)
        focus_ses <- sqrt(diag(V_relative))
        names(focus_ses) <- focus_levels
      } else {
        # Fallback: no intercept in model (unusual)
        # Use simplified method for coefficients only
        focus_vcov <- all_vcov[focus_coeff_indices, focus_coeff_indices, drop = FALSE]
        focus_ses_nonref <- sqrt(diag(focus_vcov))

        # Create SE vector for all levels
        focus_ses <- numeric(n_levels)
        names(focus_ses) <- focus_levels
        focus_ses[1] <- 0 # Reference level SE = 0 (fallback only)
        focus_ses[coeff_level_names] <- focus_ses_nonref
      }
    } else {
      # All levels are reference level (shouldn't happen for factors)
      focus_ses <- rep(0, n_levels)
      names(focus_ses) <- focus_levels
    }

    # Calculate relative coefficients
    base_coeff <- mean(focus_coeffs)
    relative_coeffs <- focus_coeffs - base_coeff

    # Calculate indices and CIs using coefficient method
    # CI multiplier: use z-score for consistency with confidence_level
    ci_multiplier <- qnorm(1 - (1 - confidence_level) / 2)

    # For subset analysis, filter coefficient results to match levels in stan_df
    if (!is.null(subset_var) && !is.null(subset_value)) {
      # Get levels present in subset data
      subset_levels <- as.character(stan_df$level)

      # Filter coefficients and SEs to match subset levels
      relative_coeffs_subset <- relative_coeffs[subset_levels]
      focus_ses_subset <- focus_ses[subset_levels]

      if (preserve_probability_scale) {
        # For binomial with raw rescaling, convert coefficients back to probability scale
        # but don't apply relative scaling
        stan_df$standardised_index <- plogis(relative_coeffs_subset + base_coeff)
        stan_df$stan_lower <- plogis(relative_coeffs_subset + base_coeff - ci_multiplier * focus_ses_subset)
        stan_df$stan_upper <- plogis(relative_coeffs_subset + base_coeff + ci_multiplier * focus_ses_subset)

        # Calculate SE on probability scale using delta method
        linear_pred <- relative_coeffs_subset + base_coeff
        prob_pred <- plogis(linear_pred)
        stan_df$stan_se <- focus_ses_subset * prob_pred * (1 - prob_pred) # Delta method for logit transform
        stan_df$standardised_cv <- ifelse(stan_df$standardised_index > 1e-6 & stan_df$standardised_index < (1 - 1e-6),
          stan_df$stan_se / stan_df$standardised_index,
          NA_real_
        )
      } else {
        # Standard relative index calculation
        stan_df$standardised_index <- exp(relative_coeffs_subset)
        stan_df$stan_lower <- exp(relative_coeffs_subset - ci_multiplier * focus_ses_subset)
        stan_df$stan_upper <- exp(relative_coeffs_subset + ci_multiplier * focus_ses_subset)

        # Calculate SE and CV on response scale (approximation)
        stan_df$stan_se <- focus_ses_subset * stan_df$standardised_index # Delta method approximation
        stan_df$standardised_cv <- focus_ses_subset # CV ~ SE on log scale for coefficient method
      }
    } else {
      # Normal analysis - use all coefficients
      if (preserve_probability_scale) {
        # For binomial with raw rescaling, convert coefficients back to probability scale
        stan_df$standardised_index <- plogis(relative_coeffs + base_coeff)
        stan_df$stan_lower <- plogis(relative_coeffs + base_coeff - ci_multiplier * focus_ses)
        stan_df$stan_upper <- plogis(relative_coeffs + base_coeff + ci_multiplier * focus_ses)

        # Calculate SE on probability scale using delta method
        linear_pred <- relative_coeffs + base_coeff
        prob_pred <- plogis(linear_pred)
        stan_df$stan_se <- focus_ses * prob_pred * (1 - prob_pred) # Delta method for logit transform
        stan_df$standardised_cv <- ifelse(stan_df$standardised_index > 1e-6 & stan_df$standardised_index < (1 - 1e-6),
          stan_df$stan_se / stan_df$standardised_index,
          NA_real_
        )
      } else {
        # Standard relative index calculation
        stan_df$standardised_index <- exp(relative_coeffs)
        stan_df$stan_lower <- exp(relative_coeffs - ci_multiplier * focus_ses)
        stan_df$stan_upper <- exp(relative_coeffs + ci_multiplier * focus_ses)

        # Calculate SE and CV on response scale (approximation)
        stan_df$stan_se <- focus_ses * stan_df$standardised_index # Delta method approximation
        stan_df$standardised_cv <- focus_ses # CV ~ SE on log scale for coefficient method
      }
    }
  } else {
      # --- PREDICTION-BASED APPROACH (modern method) ---
    message("Using prediction-based CI calculation (modern approach)")
    
    # Special handling for Gamma with inverse link - these can produce NaN values
    is_gamma_inverse <- !is.null(obj$model$family) && 
                        obj$model$family$family == "Gamma" && 
                        obj$model$family$link == "inverse"
    
    if (is_gamma_inverse) {
      message("Detected Gamma model with inverse link - using robust prediction handling")
      # Replace any NaN or Inf values with NA for safer processing
      stan_df$pred[!is.finite(stan_df$pred)] <- NA
      stan_df$se[!is.finite(stan_df$se)] <- NA
      
      # Use median instead of mean for base value calculation if there are NAs
      if (sum(is.na(stan_df$pred)) > 0) {
        message("Found non-finite predictions - using median of finite values for base calculation")
        base_pred <- median(stan_df$pred, na.rm = TRUE)
      }
    }

    if (preserve_probability_scale) {
      # For binomial with raw rescaling, preserve original probability scale
      stan_df$standardised_index <- stan_df$pred      # Calculate confidence intervals directly on probability scale
      alpha <- 1 - confidence_level
      z_score <- qnorm(1 - alpha / 2)

      # For probabilities, ensure CIs stay within [0,1]
      stan_df$stan_lower <- pmax(0, pmin(1, stan_df$pred - z_score * stan_df$se))
      stan_df$stan_upper <- pmax(0, pmin(1, stan_df$pred + z_score * stan_df$se))

      # Calculate SE and CV on probability scale
      stan_df$stan_se <- stan_df$se
      stan_df$standardised_cv <- ifelse(stan_df$standardised_index > 1e-6 & stan_df$standardised_index < (1 - 1e-6),
        stan_df$stan_se / stan_df$standardised_index,
        NA_real_
      )
    } else {
      # Standard relative index calculation
      stan_df$standardised_index <- stan_df$pred / base_pred

      # Calculate confidence intervals on the response scale
      alpha <- 1 - confidence_level
      z_score <- qnorm(1 - alpha / 2)

      # Calculate confidence intervals directly on response scale
      stan_df$stan_lower <- (stan_df$pred - z_score * stan_df$se) / base_pred
      stan_df$stan_upper <- (stan_df$pred + z_score * stan_df$se) / base_pred

      # Ensure confidence intervals are non-negative for positive-valued responses
      if (all(stan_df$pred > 0, na.rm = TRUE)) {
        stan_df$stan_lower <- pmax(stan_df$stan_lower, 0, na.rm = TRUE)
      }

      # Special handling for Gamma with inverse link - handle potential NaN values
      if (is_gamma_inverse) {
        # For Gamma(inverse), ensure all indices are valid by handling NAs
        stan_df$stan_lower[!is.finite(stan_df$stan_lower)] <- NA
        stan_df$stan_upper[!is.finite(stan_df$stan_upper)] <- NA
        
        # Use robust method for handling potential NaN values in standardized indices
        stan_df$standardised_index[!is.finite(stan_df$standardised_index)] <- NA
      }

      # Calculate CV and SE on the response scale
      stan_df$stan_se <- stan_df$se / base_pred

      # Determine if we should use delta method for CV calculation
      # Use delta method for log-link models (regardless of islog setting)
      has_log_link <- !is.null(obj$model$family$link) && obj$model$family$link == "log"

      if (has_log_link) {
        # For log-link models (e.g., Gamma(link="log"), Poisson(link="log"))
        # Use delta method for CV calculation - this is mathematically correct
        # regardless of whether islog=TRUE or FALSE
        log_preds <- predict(obj$model, newdata = obj$data, type = "link", se.fit = TRUE)
        log_pred_df <- data.frame(
          level = obj$data[[obj$focus]],
          log_se = log_preds$se.fit
        )
        log_se_agg <- aggregate(log_se ~ level, data = log_pred_df, FUN = mean)

        # Merge log SE with standardised data
        stan_df <- merge(stan_df, log_se_agg, by = "level", all.x = TRUE)

        # Delta method CV: sqrt(exp(sigma^2) - 1) where sigma is SE on log scale
        stan_df$standardised_cv <- ifelse(stan_df$standardised_index > 1e-6,
          sqrt(exp(stan_df$log_se^2) - 1),
          NA_real_
        )

        # Remove the temporary log_se column
        stan_df$log_se <- NULL
      } else {
        # For non-log-link models (Gaussian, etc.), use standard CV = se/mean
        stan_df$standardised_cv <- ifelse(stan_df$standardised_index > 1e-6,
          stan_df$stan_se / stan_df$standardised_index,
          NA_real_
        )
      }
    }
  } # End of else block for prediction-based approach

  # Apply rescaling to unstandardised and standardised indices
  # For binomial models, always preserve probability scale for unstandardised indices
  if (!preserve_unstandardised_probability_scale) {
    # Apply rescaling to unstandardised indices for non-binomial cases
    indices_df$unstan <- rescale_index(indices_df$unstan, rescale_method, custom_rescale_value)
  }

  # For binomial models with raw rescaling, skip rescaling of standardised indices
  # to preserve probability scale
  if (!preserve_probability_scale) {
    # Store original standardised values for rescaling bounds
    original_stan <- stan_df$standardised_index
    original_lower <- stan_df$stan_lower
    original_upper <- stan_df$stan_upper

    # Rescale the central estimate
    stan_df$standardised_index <- rescale_index(stan_df$standardised_index, rescale_method, custom_rescale_value)

    # For confidence intervals, maintain the same relative distance from the central estimate
    # This preserves the statistical relationship
    if (rescale_method != "raw") {
      rescale_factor <- stan_df$standardised_index / original_stan
      stan_df$stan_lower <- original_lower * rescale_factor
      stan_df$stan_upper <- original_upper * rescale_factor
    } else {
      stan_df$stan_lower <- rescale_index(stan_df$stan_lower, rescale_method, custom_rescale_value)
      stan_df$stan_upper <- rescale_index(stan_df$stan_upper, rescale_method, custom_rescale_value)
    }
  }

  # Handle special case: levels with zero standard error (reference levels)
  zero_se_levels <- abs(stan_df$se) < 1e-10
  if (any(zero_se_levels)) {
    stan_df$stan_lower[zero_se_levels] <- stan_df$standardised_index[zero_se_levels]
    stan_df$stan_upper[zero_se_levels] <- stan_df$standardised_index[zero_se_levels]
    stan_df$stan_se[zero_se_levels] <- 0
    stan_df$standardised_cv[zero_se_levels] <- 0
  }

  # Always ensure CV is correctly calculated as se/index after any rescaling
  stan_df$standardised_cv <- ifelse(stan_df$standardised_index > 1e-6,
    stan_df$stan_se / stan_df$standardised_index,
    NA_real_
  )

  # Merge standardised results with unstandardised indices
  merge_cols <- c("level", "standardised_index", "stan_lower", "stan_upper", "stan_se", "standardised_cv")
  indices_df <- merge(indices_df, stan_df[, merge_cols], by = "level")

  # --- 3. Step-wise Calculations ---
  # Build the model one term at a time to see how the index and model fit change.
  # Note: For subset analysis, stepwise calculations are skipped to avoid
  # model fitting issues with different data
  summary_list <- list()
  step_indices_list <- list()

  # Define all_terms (ordered) for stepwise building.
  # stats::terms() on a gam object only returns parametric-style labels which
  # can collapse multiple smooths on the same variable (e.g. s(depth,by=gear,...) + s(depth,...) -> "depth").
  # To correctly attribute smooth-specific EDF we must rebuild models adding each
  # explicit RHS token (parametric term or full smooth specification) in the original order.
  # We therefore parse the RHS of the model formula manually, splitting on top-level '+'
  # while respecting parentheses.
  parse_formula_terms <- function(frm) {
    rhs <- deparse(formula(frm)[[3]])
    # Collapse multi-line deparse into single string
    rhs <- paste(rhs, collapse = " ")
    chars <- strsplit(rhs, "")[[1]]
    depth_paren <- 0
    buf <- ""
    tokens <- c()
    for (ch in chars) {
      if (ch == "(") depth_paren <- depth_paren + 1
      if (ch == ")") depth_paren <- max(0, depth_paren - 1)
      if (ch == "+" && depth_paren == 0) {
        token <- trimws(buf)
        if (nzchar(token)) tokens <- c(tokens, token)
        buf <- ""
      } else {
        buf <- paste0(buf, ch)
      }
    }
    last_token <- trimws(buf)
    if (nzchar(last_token)) tokens <- c(tokens, last_token)
    # Remove empty tokens and any inline comments (unlikely)
    tokens <- tokens[nzchar(tokens)]
    tokens
  }
  all_terms <- parse_formula_terms(obj$model)

  # Check if this is subset analysis
  is_subset_analysis <- !is.null(subset_var) && !is.null(subset_value)

  if (!is_subset_analysis) {
    # Perform stepwise analysis for full dataset analysis

    # Start with an intercept-only model
    model_int <- update(obj$model, . ~ 1, data = obj$data)
    logLike_int <- as.numeric(logLik(model_int))
    summary_list[["intercept"]] <- data.frame(
      term = "Intercept",
      logLike = logLike_int,
      aic = AIC(model_int),
      r_sq = 0,
      deviance_explained = 0,
      df = 0,
      smooth_edf = 0,
      smooth_edf_cum = 0
    )
    prev_resid_df <- model_int$df.residual
    prev_total_smooth_edf <- 0

    # Sequentially add each term from the original model formula
    for (i in seq_along(all_terms)) {
      current_terms <- all_terms[1:i]
      formula_str <- paste("~", paste(current_terms, collapse = " + "))
      model_step <- update(obj$model, as.formula(formula_str), data = obj$data)

      # Store the step-wise index if the focus term is in the current model
      if (obj$focus %in% current_terms) {
        step_preds <- predict(model_step, newdata = obj$data, type = "terms")
        pred_cols <- colnames(step_preds)
        # Focus columns: exact match for factor parametric term OR smooths containing the focus variable
        # Simpler: gather columns that start with focus term or contain focus term followed by ) for smooths
        focus_cols <- grep(paste0("^", obj$focus, "$|^", obj$focus, "(:|\\)|$)"), pred_cols, value = TRUE)
        if (length(focus_cols) == 0) {
          # Fallback: any column that contains the focus term as standalone word
          focus_cols <- grep(paste0("(^|:)", obj$focus, "($|:)"), pred_cols, value = TRUE)
        }
        if (length(focus_cols) == 0) {
          stop(paste0("Could not find any columns for focus term '", obj$focus, "' in step-wise model predictions (type='terms')."))
        }
        focus_effect_sum <- rowSums(step_preds[, focus_cols, drop = FALSE])
        step_focus_effects <- aggregate(focus_effect_sum ~ obj$data[[obj$focus]], FUN = mean)
        names(step_focus_effects) <- c("level", "effect")
        base_step <- mean(step_focus_effects$effect)
        step_focus_effects$index <- exp(step_focus_effects$effect - base_step)
        term_label <- paste("+", all_terms[i])
        step_indices_list[[term_label]] <- step_focus_effects[, c("level", "index")]
        names(step_indices_list[[term_label]])[2] <- term_label
      }

      # Store summary statistics for the current step
      logLike_step <- as.numeric(logLik(model_step))
      model_summary <- summary(model_step)

      # Get R-squared and Deviance Explained, robust to model type (gam vs. glm)
      r_sq_step <- if ("r.sq" %in% names(model_summary)) model_summary$r.sq else NA
      dev_expl_step <- if ("dev.expl" %in% names(model_summary)) model_summary$dev.expl else (model_step$null.deviance - deviance(model_step)) / model_step$null.deviance

      # Degrees of freedom contributed by this added term (drop in residual df)
      df_added <- prev_resid_df - model_step$df.residual
      if (is.null(df_added) || is.na(df_added)) df_added <- NA_real_

      # Smooth-specific EDF: obtain total smooth edf from the model summary
      s_tab <- try(summary(model_step)$s.table, silent = TRUE)
      total_smooth_edf <- if (!inherits(s_tab, "try-error") && !is.null(s_tab)) {
        # s.table may be matrix or data.frame; edf column usually named 'edf'
        if (is.matrix(s_tab) || is.data.frame(s_tab)) {
          edf_col <- NULL
          if ("edf" %in% colnames(s_tab)) edf_col <- s_tab[, "edf"] else if (ncol(s_tab) >= 1) edf_col <- s_tab[, 1]
          suppressWarnings(sum(as.numeric(edf_col), na.rm = TRUE))
        } else {
          0
        }
      } else {
        0
      }
      smooth_edf_added <- total_smooth_edf - prev_total_smooth_edf
      if (is.null(smooth_edf_added) || is.na(smooth_edf_added)) smooth_edf_added <- NA_real_
      summary_list[[all_terms[i]]] <- data.frame(
        term = all_terms[i],
        logLike = logLike_step,
        aic = AIC(model_step),
        r_sq = r_sq_step,
        deviance_explained = dev_expl_step,
        df = df_added,
        smooth_edf = smooth_edf_added,
        smooth_edf_cum = total_smooth_edf
      )
      prev_resid_df <- model_step$df.residual
      prev_total_smooth_edf <- total_smooth_edf
    }
  } else {
    # For subset analysis, perform stepwise analysis on subset data
    message("Performing stepwise analysis on subset data")

    # Pre-filter terms to exclude those with insufficient factor levels in subset data
    valid_subset_terms <- sapply(all_terms, function(term) {
      # Extract variable names from the term (handles interactions and smooths)
      var_names <- all.vars(as.formula(paste("~", term)))

      # Check each variable in the term
      all(sapply(var_names, function(var_name) {
        if (var_name %in% names(obj$data)) {
          var_data <- obj$data[[var_name]]
          if (is.factor(var_data) || is.character(var_data)) {
            # For factors/characters, need at least 2 levels for contrasts
            length(unique(var_data[!is.na(var_data)])) >= 2
          } else {
            # For numeric variables, just need some variation
            length(unique(var_data[!is.na(var_data)])) >= 2
          }
        } else {
          TRUE # Variable not in data, let the model fitting handle it
        }
      }))
    })

    # Filter out invalid terms and report what was excluded
    invalid_terms <- all_terms[!valid_subset_terms]
    if (length(invalid_terms) > 0) {
      message(
        "Excluding terms from subset stepwise analysis due to insufficient factor levels: ",
        paste(invalid_terms, collapse = ", ")
      )

      # Record the excluded terms in the summary
      for (invalid_term in invalid_terms) {
        summary_list[[paste(invalid_term, "(Subset - Excluded)")]] <- data.frame(
          term = paste(invalid_term, "(Subset - Excluded: insufficient levels)"),
          logLike = NA,
          aic = NA,
          r_sq = NA,
          deviance_explained = NA
        )
      }
    }

    # Use only valid terms for stepwise analysis
    subset_terms <- all_terms[valid_subset_terms]

    # Start with an intercept-only model fitted to subset data
    tryCatch(
      {
        model_int <- update(obj$model, . ~ 1, data = obj$data)
        logLike_int <- as.numeric(logLik(model_int))
        summary_list[["intercept"]] <- data.frame(
          term = "Intercept (Subset)",
          logLike = logLike_int,
          aic = AIC(model_int),
          r_sq = 0,
          deviance_explained = 0,
          df = 0,
          smooth_edf = 0,
          smooth_edf_cum = 0
        )
        prev_resid_df <- model_int$df.residual
        prev_total_smooth_edf <- 0

        # Sequentially add each valid term from the filtered list
        for (i in seq_along(subset_terms)) {
          current_terms <- subset_terms[1:i]

          formula_str <- paste("~", paste(current_terms, collapse = " + "))

          tryCatch(
            {
              # Fit stepwise model to subset data
              model_step <- update(obj$model, as.formula(formula_str), data = obj$data)

              # Store the step-wise index if the focus term is in the current model
              if (obj$focus %in% current_terms) {
                step_preds <- predict(model_step, newdata = obj$data, type = "terms")
                focus_cols <- grep(paste0("^", obj$focus, ""), colnames(step_preds), value = TRUE)
                if (length(focus_cols) > 0) {
                  focus_effect_sum <- rowSums(step_preds[, focus_cols, drop = FALSE])
                  step_focus_effects <- aggregate(focus_effect_sum ~ obj$data[[obj$focus]], FUN = mean)
                  names(step_focus_effects) <- c("level", "effect")

                  base_step <- mean(step_focus_effects$effect)
                  step_focus_effects$index <- exp(step_focus_effects$effect - base_step)

                  # Name the column after the term that was just added
                  term_label <- paste("+", subset_terms[i])
                  step_indices_list[[term_label]] <- step_focus_effects[, c("level", "index")]
                  names(step_indices_list[[term_label]])[2] <- term_label
                }
              }

              # Store summary statistics for the current step
              logLike_step <- as.numeric(logLik(model_step))
              model_summary <- summary(model_step)

              # Get R-squared and Deviance Explained, robust to model type (gam vs. glm)
              r_sq_step <- if ("r.sq" %in% names(model_summary)) model_summary$r.sq else NA
              dev_expl_step <- if ("dev.expl" %in% names(model_summary)) model_summary$dev.expl else (model_step$null.deviance - deviance(model_step)) / model_step$null.deviance

              df_added <- prev_resid_df - model_step$df.residual
              if (is.null(df_added) || is.na(df_added)) df_added <- NA_real_
              s_tab <- try(summary(model_step)$s.table, silent = TRUE)
              total_smooth_edf <- if (!inherits(s_tab, "try-error") && !is.null(s_tab)) {
                if (is.matrix(s_tab) || is.data.frame(s_tab)) {
                  edf_col <- NULL
                  if ("edf" %in% colnames(s_tab)) edf_col <- s_tab[, "edf"] else if (ncol(s_tab) >= 1) edf_col <- s_tab[, 1]
                  suppressWarnings(sum(as.numeric(edf_col), na.rm = TRUE))
                } else {
                  0
                }
              } else {
                0
              }
              smooth_edf_added <- total_smooth_edf - prev_total_smooth_edf
              if (is.null(smooth_edf_added) || is.na(smooth_edf_added)) smooth_edf_added <- NA_real_
              summary_list[[paste(subset_terms[i], "(Subset)")]] <- data.frame(
                term = paste(subset_terms[i], "(Subset)"),
                logLike = logLike_step,
                aic = AIC(model_step),
                r_sq = r_sq_step,
                deviance_explained = dev_expl_step,
                df = df_added,
                smooth_edf = smooth_edf_added,
                smooth_edf_cum = total_smooth_edf
              )
              prev_resid_df <- model_step$df.residual
              prev_total_smooth_edf <- total_smooth_edf
            },
            error = function(e) {
              # If model fitting fails for this step, record the error and continue
              warning(paste(
                "Failed to fit stepwise model with terms:", paste(current_terms, collapse = " + "),
                "Error:", e$message
              ))
              summary_list[[paste(subset_terms[i], "(Subset - Failed)")]] <<- data.frame(
                term = paste(subset_terms[i], "(Subset - Failed)"),
                logLike = NA,
                aic = NA,
                r_sq = NA,
                deviance_explained = NA,
                df = NA,
                smooth_edf = NA,
                smooth_edf_cum = NA
              )
            }
          )
        }
      },
      error = function(e) {
        # If even the intercept model fails, fall back to minimal summary
        warning("Could not perform stepwise analysis on subset data:", e$message)
        summary_list[["full_model"]] <- data.frame(
          term = "Full Model (Subset Analysis)",
          logLike = as.numeric(logLik(obj$model)),
          aic = AIC(obj$model),
          r_sq = summary(obj$model)$r.sq,
          deviance_explained = summary(obj$model)$dev.expl,
          df = NA,
          smooth_edf = NA,
          smooth_edf_cum = NA
        )
      }
    )
  }

  summary_df <- do.call(rbind, summary_list)
  # Calculate successive differences
  # If we performed a full (non-subset) analysis, ensure the final row reflects the full model metrics exactly
  if (!is_subset_analysis && nrow(summary_df) > 1) {
    full_mod_sum <- try(suppressWarnings(summary(obj$model)), silent = TRUE)
    if (!inherits(full_mod_sum, "try-error")) {
      if ("r.sq" %in% names(full_mod_sum)) {
        summary_df$r_sq[nrow(summary_df)] <- full_mod_sum$r.sq
      }
      if ("dev.expl" %in% names(full_mod_sum)) {
        summary_df$deviance_explained[nrow(summary_df)] <- full_mod_sum$dev.expl
      }
    }
  }
  summary_df$r_sq_diff <- c(0, diff(summary_df$r_sq))
  summary_df$deviance_explained_diff <- c(0, diff(summary_df$deviance_explained))
  # Clamp tiny floating point artifacts to zero for readability
  clamp_tol <- 1e-6
  clamp_cols <- intersect(c("df", "smooth_edf", "smooth_edf_cum", "r_sq_diff", "deviance_explained_diff"), names(summary_df))
  for (cc in clamp_cols) {
    bad <- !is.na(summary_df[[cc]]) & abs(summary_df[[cc]]) < clamp_tol
    if (any(bad)) summary_df[[cc]][bad] <- 0
  }
  rownames(summary_df) <- NULL

  # Combine all step-wise indices into a single data frame
  if (length(step_indices_list) > 0) {
    step_indices_df <- Reduce(function(df1, df2) merge(df1, df2, by = "level"), step_indices_list)
    indices_df <- merge(indices_df, step_indices_df, by = "level")
    # Define step_cols as the names of the stepwise index columns (those starting with '+')
    step_cols <- names(step_indices_df)[names(step_indices_df) != "level"]
  } else {
    step_cols <- character(0)
  }

  # --- 4. Influence Calculations ---
  # Calculate how each non-focus term influences the focus term's index.
  # Use terms predictions for influence analysis to isolate individual term effects
  all_preds_df <- as.data.frame(preds_full_terms$fit)
  all_preds_df$level <- obj$data[[obj$focus]]

  # After creating all_preds_df, also create a dataframe for standard errors
  all_preds_se_df <- as.data.frame(preds_full_terms$se.fit)
  all_preds_se_df$level <- obj$data[[obj$focus]]

  influ_list <- list()
  # store setting for downstream functions
  obj$separate_by_smooth <- separate_by_smooth
  non_focus_terms <- setdiff(all_terms, obj$focus)

  for (term in non_focus_terms) {
    # Use the robust helper to find all relevant columns
    term_cols <- find_term_columns(term, colnames(all_preds_df), group_by_bysmooth = !separate_by_smooth)
    if (length(term_cols) == 0) next # Skip if no columns found for this term

    # Sum across all columns for this term (handles smooths and factors)
    term_effect_sum <- rowSums(all_preds_df[, term_cols, drop = FALSE])
    infl <- aggregate(term_effect_sum ~ all_preds_df$level, FUN = mean)
    names(infl) <- c("level", "influence")

    # Overall influence: Average magnitude of the effect
    overall_infl <- exp(mean(abs(infl$influence))) - 1

    # Trend influence: How the effect changes systematically across focus levels
    level_num <- as.numeric(infl$level)
    # Add safety checks for subset analysis
    if (any(is.na(level_num)) || any(is.na(infl$influence)) ||
      length(level_num) < 2 || is.na(var(level_num)) || var(level_num) == 0) {
      trend_infl <- 0
    } else {
      trend_infl <- exp(cov(level_num, infl$influence) / var(level_num)) - 1
    }

    summary_df$overall[summary_df$term == term] <- overall_infl
    summary_df$trend[summary_df$term == term] <- trend_infl

    base_infl <- mean(infl$influence)
    infl$influence <- exp(infl$influence - base_infl)

    infl$term <- term
    influ_list[[term]] <- infl
  }

  influences_df <- do.call(rbind, influ_list)
  if (!is.null(influences_df)) rownames(influences_df) <- NULL

  # Create step_labels: names are step_cols, values are the full term expressions from the model
  step_labels <- tryCatch(
    {
      if (is.null(obj$model$formula)) {
        setNames(step_cols, step_cols)
      } else {
        # Extract terms from model formula
        terms_in_formula <- attr(terms(obj$model), "term.labels")

        # Create mapping from step_cols to full formula terms
        setNames(
          sapply(step_cols, function(col) {
            # Remove leading '+' and find matching term
            term_pattern <- sub("^\\+", "", col)

            # First try exact match
            exact_matches <- terms_in_formula[terms_in_formula == term_pattern]
            if (length(exact_matches) > 0) {
              return(exact_matches[1])
            }

            # Try fixed string match (original approach)
            fixed_matches <- grep(term_pattern, terms_in_formula, value = TRUE, fixed = TRUE)
            if (length(fixed_matches) > 0) {
              return(fixed_matches[1])
            }

            # Try whitespace-normalized matching
            pattern_stripped <- strip_whitespace(term_pattern)
            formula_stripped <- strip_whitespace(terms_in_formula)
            match_idx <- which(formula_stripped == pattern_stripped)
            if (length(match_idx) > 0) {
              return(terms_in_formula[match_idx[1]])
            }

            # Fall back to original column name if no match
            col
          }),
          step_cols
        )
      }
    },
    error = function(e) {
      warning("Could not create step labels, using column names instead")
      setNames(step_cols, step_cols)
    }
  )

  # Rename standardised index columns for consistency with package API
  if ("stan" %in% names(indices_df)) {
    names(indices_df)[names(indices_df) == "stan"] <- "standardised_index"
  }
  if ("stanLower" %in% names(indices_df)) {
    names(indices_df)[names(indices_df) == "stanLower"] <- "stan_lower"
  }
  if ("stanUpper" %in% names(indices_df)) {
    names(indices_df)[names(indices_df) == "stanUpper"] <- "stan_upper"
  }

  # --- Finalize ---
  # Restore original data if subset analysis was performed
  if (!is.null(subset_var) && !is.null(subset_value) && exists("original_data", inherits = FALSE)) {
    obj$data <- original_data
  }

  # Store all calculated data frames in the object's 'calculated' slot
  obj$calculated <- list(
    indices = indices_df,
    summary = summary_df,
    influences = influences_df,
    predictions = all_preds_df,
    prediction_se = all_preds_se_df,
    step_labels = step_labels
  )

  return(obj)
}

#' @title Validate gam_influence object
#' @description Internal function to validate gam_influence objects and check family compatibility
#' @param obj A gam_influence object to validate
#' @return TRUE if valid, otherwise throws an error
#' @noRd
validate_gam_influence <- function(obj) {
  if (!inherits(obj, "gam_influence")) {
    stop("Object must be of class 'gam_influence'", call. = FALSE)
  }
  if (!inherits(obj$model, "gam")) {
    stop("Model must be a GAM object from mgcv package", call. = FALSE)
  }
  if (!obj$focus %in% obj$terms) {
    stop("Focus term must be present in model terms", call. = FALSE)
  }
  if (nrow(obj$data) == 0) {
    stop("Data cannot be empty", call. = FALSE)
  }

  # Check family compatibility
  if (!is.null(obj$model$family)) {
    family_name <- obj$model$family$family
    supported_families <- c("gaussian", "binomial", "Gamma", "gamma", "poisson", "weibull", "quasi", "quasipoisson", "quasibinomial", "Tweedie", "nb")

    # Check for Tweedie family (which reports as "Tweedie(p=X.XXX)")
    is_tweedie <- grepl("^Tweedie\\(", family_name)
    # Check for negative binomial family (which reports as "Negative Binomial(X.XXX)")
    is_nb <- grepl("^Negative Binomial\\(", family_name)

    if (!family_name %in% supported_families && !is_tweedie && !is_nb) {
      warning("Family '", family_name, "' may not be fully supported. Supported families: ",
        paste(c(supported_families, "Tweedie(p=X)", "Negative Binomial(X)"), collapse = ", "),
        call. = FALSE
      )
    }

    # Check for common issues
    observed <- obj$data[[obj$response]]
    if (family_name == "binomial" && !all(observed >= 0 & observed <= 1)) {
      if (!all(observed %in% c(0, 1))) {
        warning("Binomial family detected but response values are not between 0-1 or binary. This may cause issues.", call. = FALSE)
      }
    }

    if (family_name %in% c("Gamma", "gamma") && any(observed <= 0)) {
      warning("Gamma family detected but response contains non-positive values. This may cause issues.", call. = FALSE)
    }

    if (family_name == "weibull" && any(observed <= 0)) {
      warning("Weibull family detected but response contains non-positive values. This may cause issues.", call. = FALSE)
    }

    if (family_name == "poisson" && any(observed < 0)) {
      warning("Poisson family detected but response contains negative values. This may cause issues.", call. = FALSE)
    }

    # Check for Tweedie family requirements
    if (is_tweedie && any(observed < 0)) {
      warning("Tweedie family detected but response contains negative values. Tweedie requires non-negative values.", call. = FALSE)
    }

    # Check for negative binomial family requirements
    if (is_nb && (any(observed < 0) || any(observed != round(observed)))) {
      warning("Negative Binomial family detected but response contains negative or non-integer values. This may cause issues.", call. = FALSE)
    }
  }

  TRUE
}

#' @title Helper to robustly match model terms to prediction columns
#' @description Internal helper function to find column names that correspond to a given term.
#' @param term Character string of the term name to match.
#' @param colnames_vec Character vector of column names to search in.
#' @return Character vector of matching column names.
#' @noRd
find_term_columns <- function(term, colnames_vec, group_by_bysmooth = TRUE) {
  original_term <- term
  # If term is a smooth with arguments, strip spaces
  term_clean <- gsub("[[:space:]]+", "", term)
  # Generic normalisation for smooth constructors s(), te(), ti(), t2()
  normalize_smooth <- function(txt) {
    if (!grepl("^(s|te|ti|t2)\\(", txt)) {
      return(txt)
    }
    # Remove all whitespace for parsing
    no_space <- gsub("[[:space:]]+", "", txt)
    inner <- sub("^[a-z0-9]+\\((.*)\\)$", "\\1", no_space, ignore.case = TRUE)
    parts <- strsplit(inner, ",")[[1]]
    # Collect leading variable tokens up until the first argument containing '=' or 'c('
    # This avoids accidentally capturing fragments of k=c(6,6) etc.
    var_tokens <- c()
    for (p in parts) {
      if (grepl("=", p) || grepl("^c\\(", p)) break
      var_tokens <- c(var_tokens, p)
    }
    if (length(var_tokens) == 0) {
      return(txt)
    }
    builder <- sub("^([a-z0-9]+)\\(.*$", "\\1", no_space, ignore.case = TRUE)
    paste0(builder, "(", paste(var_tokens, collapse = ","), ")")
  }
  normalized_label <- normalize_smooth(term)
  # Direct attempt: if normalized smooth column present
  if (normalized_label %in% colnames_vec) {
    return(normalized_label)
  }
  # For tensor products predict() often drops spaces after commas, so try that variant
  compressed_label <- gsub("[[:space:]]+", "", normalized_label)
  if (compressed_label %in% colnames_vec) {
    return(compressed_label)
  }
  # If still not found and it's a tensor product, attempt partial match (covers by-smooth expansion in future)
  if (grepl("^(te|ti|t2)\\(", term_clean)) {
    # Normalize: drop spaces then rebuild minimal label with just variable list
    tensor_base <- normalize_smooth(term_clean)
    tensor_base <- gsub("[[:space:]]+", "", tensor_base)
    # Direct exact match
    if (tensor_base %in% colnames_vec) {
      return(tensor_base)
    }
    # Escape regex metachars (omit braces to avoid TRE 'Invalid contents of {}') for prefix search
    esc_tensor <- gsub("([().,+*?^$|\\[\\]{}\\\\])", "\\\\\\1", tensor_base)
    tensor_cols <- grep(paste0("^", esc_tensor), colnames_vec, value = TRUE)
    if (length(tensor_cols) > 0) {
      return(tensor_cols)
    }
  }
  # Remove k=... and other arguments inside smooth for matching prediction columns
  if (grepl("^s\\(", term_clean)) {
    # remove ,k=... patterns
    term_base <- sub(",k=[^,)]*", "", term_clean)
    term_base <- sub(",bs=[^,)]*", "", term_base)
    term_base <- sub(",by=[^,)]*", "", term_base)
    term_base <- sub("\\)$", "", term_base)
    # Extract variable list inside s()
    inner <- sub("^s\\((.*)$", "\\1", term_base)
    inner_vars <- strsplit(inner, ",")[[1]]
    if (length(inner_vars) > 0) {
      term_variable <- inner_vars[1]
      # Prediction columns for simple smooths are like s(var)
      simple_pattern <- paste0("^s\\(", term_variable, ".*\\)")
      simple_cols <- grep(simple_pattern, colnames_vec, value = TRUE)
      if (length(simple_cols) > 0) {
        return(simple_cols)
      }
    }
  }
  # Handle interaction terms: convert * to : for prediction column matching
  # In GAM formulas, A*B creates interactions but predict() uses A:B format
  if (grepl("\\*", term) && !grepl("^(s|te|ti|t2)\\(", term)) {
    # Convert interaction term from * format to : format
    interaction_term <- gsub("\\s*\\*\\s*", ":", term)
    if (interaction_term %in% colnames_vec) {
      return(interaction_term)
    }
    # Also try with spaces removed
    interaction_term_clean <- gsub("[[:space:]]+", "", interaction_term)
    if (interaction_term_clean %in% colnames_vec) {
      return(interaction_term_clean)
    }

    # Handle variable order differences: A*B vs B:A
    # Split the interaction term and try both orders
    vars <- strsplit(gsub("\\s*\\*\\s*", "*", term), "\\*")[[1]]
    if (length(vars) == 2) {
      # Try original order A:B
      order1 <- paste(trimws(vars), collapse = ":")
      if (order1 %in% colnames_vec) {
        return(order1)
      }
      # Try reversed order B:A
      order2 <- paste(trimws(rev(vars)), collapse = ":")
      if (order2 %in% colnames_vec) {
        return(order2)
      }
    }
  }

  # Handle interaction terms that are already in colon format
  # Sometimes the formula parsing gives us A:B but prediction columns have B:A
  if (grepl(":", term) && !grepl("^(s|te|ti|t2)\\(", term)) {
    # First try exact match
    if (term %in% colnames_vec) {
      return(term)
    }

    # Try with spaces removed
    term_clean_colon <- gsub("[[:space:]]+", "", term)
    if (term_clean_colon %in% colnames_vec) {
      return(term_clean_colon)
    }

    # Try reversed order for colon-separated terms
    vars <- strsplit(term, ":")[[1]]
    if (length(vars) == 2) {
      # Try reversed order B:A
      reversed_term <- paste(trimws(rev(vars)), collapse = ":")
      if (reversed_term %in% colnames_vec) {
        return(reversed_term)
      }
    }
  }

  # Handle mgcv by-variable smooths (e.g. s(depth,by=gear)) whose term.labels differ
  # from the per-level columns produced by predict(..., type='terms'), e.g. columns like
  #   s(depth):gearLevel1, s(depth):gearLevel2
  # We want to group all those columns back to the single model term 's(depth,by=gear)'.
  term_nospace <- gsub("[[:space:]]+", "", term)
  if (grepl("by=", term_nospace, fixed = TRUE)) {
    if (!group_by_bysmooth) {
      # When separating by-smooths, attempt to return only columns associated with each by-smooth level
      # Keep the core smooth columns distinct from base smooth by returning per-level pattern set
      term_core <- sub(",k=.*$", "", term_nospace)
      core_smooth <- sub(",by=.*$", "", term_core)
      core_escaped <- gsub("([().,+*?^$|\\[\\]{}\\\\])", "\\\\\\1", core_smooth, perl = TRUE)
      pattern_core <- paste0("^", core_escaped, ":")
      candidate_cols <- grep(pattern_core, colnames_vec, value = TRUE)
      if (length(candidate_cols) > 0) {
        return(candidate_cols)
      }
      # fallback to generic flow if none found
    }
    # Remove trailing k= specifications first (e.g. ,k=6)
    term_core <- sub(",k=.*$", "", term_nospace)
    # Extract smooth without by-part
    core_smooth <- sub(",by=.*$", "", term_core)
    by_var <- sub(".*by=([^,)+]+).*$", "\\1", term_core)
    # Build regex for columns: coreSmooth:byVarLevel  (mgcv drops ',by=var')
    core_escaped <- gsub("([().,+*?^$|\\[\\]{}\\\\])", "\\\\\\1", core_smooth, perl = TRUE)
    # Columns start with the core smooth then ':'
    pattern_core <- paste0("^", core_escaped, ":")
    candidate_cols <- grep(pattern_core, colnames_vec, value = TRUE)
    if (length(candidate_cols) == 0) {
      # Fallback: columns containing ':' then by var name (coarse)
      by_pattern <- paste0(":", by_var)
      candidate_cols <- grep(by_pattern, colnames_vec, value = TRUE, fixed = TRUE)
    }
    if (length(candidate_cols) > 0) {
      return(candidate_cols)
    }
  }
  # Try exact match first
  cols <- which(colnames_vec == term)
  if (length(cols) > 0) {
    return(colnames_vec[cols])
  }
  # Try partial match (for factors, by-variables, etc.)
  cols <- grep(paste0("(^|\\(|\\:)", term, "($|\\)|\\:|\\d+)"), colnames_vec, value = TRUE)
  if (length(cols) > 0) {
    return(cols)
  }
  # Try contains (for factor levels, e.g., year1990)
  cols <- grep(term, colnames_vec, value = TRUE, fixed = TRUE)

  # Enhanced error handling
  if (length(cols) == 0) {
    available_terms <- paste(colnames_vec, collapse = ", ")
    stop(sprintf(
      "Could not find columns for term '%s'. Available columns: %s",
      term, available_terms
    ), call. = FALSE)
  }
  return(cols)
}

#' @title Geometric Mean
#' @description Calculate the geometric mean of a numeric vector, with robust handling of missing values and non-positive values.
#' @param x A numeric vector.
#' @param na.rm Logical. Should missing values be removed before calculation? Default is TRUE.
#' @return The geometric mean of the input vector.
#' @title Calculate Group Statistics
#' @description Helper function to calculate mean, standard error, and coefficient of variation for a group
#' @param x Numeric vector
#' @return Named vector with mean, se, cv
#' @noRd
calc_group_stats <- function(x) {
  n <- length(x)
  if (n <= 1) {
    return(c(mean = mean(x), se = NA_real_, cv = NA_real_))
  }
  m <- mean(x)
  se <- sd(x) / sqrt(n)
  cv <- se / abs(m)
  return(c(mean = m, se = se, cv = cv))
}

#' @title Strip Whitespace
#' @description Remove all whitespace from a string
#' @param x Character string
#' @return String with whitespace removed
#' @noRd
strip_whitespace <- function(x) gsub("\\s+", "", x)

#' @title Geometric Mean
#' @description Calculate the geometric mean of a numeric vector
#' @param x Numeric vector of positive values
#' @param na.rm Logical. Should NA values be removed?
#' @return Geometric mean value
#' @export
geometric_mean <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  if (any(x <= 0)) {
    warning("Non-positive values detected, using arithmetic mean")
    return(mean(x))
  }
  exp(mean(log(x)))
}

#' @title Rescale Index
#' @description Rescale an index using different methods.
#' @param index Numeric vector to rescale.
#' @param method Character. Rescaling method.
#' @param custom_value Numeric. Custom rescaling value when method = "custom".
#' @return Rescaled index vector.
#' @noRd
rescale_index <- function(index, method = c("auto", "geometric_mean", "arithmetic_mean", "raw", "custom"),
                          custom_value = 1) {
  method <- match.arg(method)

  # Auto should have been resolved by this point, but provide fallback
  if (method == "auto") {
    method <- "geometric_mean"
    warning("rescale_index received 'auto' method - using 'geometric_mean' as fallback")
  }

  switch(method,
    "geometric_mean" = index / geometric_mean(index),
    "arithmetic_mean" = index / mean(index),
    "raw" = index,
    "custom" = index / geometric_mean(index) * custom_value
  )
}

#' @title Calculate Index by Mean Type
#' @description Standardized helper function to calculate indices using specified mean type
#' @param observed Numeric vector of observed values
#' @param focus_var Factor or numeric vector of focus variable levels
#' @param mean_type Character. Type of mean to use: "arithmetic", "geometric", "raw"
#' @param islog Logical. Whether data is already on log scale
#' @param preserve_probability_scale Logical. Whether to preserve probability scale
#' @return Data frame with level, unstan, unstan_se, unstan_cv columns
calculate_index_by_mean_type <- function(observed, focus_var, mean_type = c("arithmetic", "geometric", "raw"),
                                         islog = FALSE, preserve_probability_scale = FALSE) {
  mean_type <- match.arg(mean_type)

  if (mean_type == "raw") {
    # Raw values - no scaling
    if (preserve_probability_scale) {
      if (islog) {
        # Data is on log scale - exp() first to get back to original scale
        exp_observed <- exp(observed)
        stats_agg <- aggregate(exp_observed, list(level = focus_var), calc_group_stats)
      } else {
        # Data is on original scale - use directly
        stats_agg <- aggregate(observed, list(level = focus_var), calc_group_stats)
      }
      agg_df <- data.frame(
        level = stats_agg$level,
        unstan = stats_agg$x[, 1], # Keep actual values
        unstan_se = stats_agg$x[, 2],
        unstan_cv = stats_agg$x[, 3]
      )
    } else {
      # Convert to relative scale even for raw
      if (islog) {
        # Data is on log scale - exp() first to get back to original scale
        exp_observed <- exp(observed)
        stats_agg <- aggregate(exp_observed, list(level = focus_var), calc_group_stats)
      } else {
        # Data is on original scale - use directly
        stats_agg <- aggregate(observed, list(level = focus_var), calc_group_stats)
      }
      base_mean <- mean(stats_agg$x[, 1])
      agg_df <- data.frame(
        level = stats_agg$level,
        unstan = stats_agg$x[, 1] / base_mean,
        unstan_se = stats_agg$x[, 2] / base_mean,
        unstan_cv = stats_agg$x[, 3]
      )
    }
  } else if (mean_type == "arithmetic") {
    # Arithmetic mean
    if (islog) {
      # Data is on log scale - exp() first to get back to original scale
      exp_observed <- exp(observed)
      stats_agg <- aggregate(exp_observed, list(level = focus_var), calc_group_stats)
      base_mean <- mean(stats_agg$x[, 1])
      agg_df <- data.frame(
        level = stats_agg$level,
        unstan = stats_agg$x[, 1] / base_mean,
        unstan_se = stats_agg$x[, 2] / base_mean,
        unstan_cv = stats_agg$x[, 3]
      )
    } else {
      # Data is on original scale - use directly
      stats_agg <- aggregate(observed, list(level = focus_var), calc_group_stats)
      base_mean <- mean(stats_agg$x[, 1])
      agg_df <- data.frame(
        level = stats_agg$level,
        unstan = stats_agg$x[, 1] / base_mean,
        unstan_se = stats_agg$x[, 2] / base_mean,
        unstan_cv = stats_agg$x[, 3]
      )
    }
  } else if (mean_type == "geometric") {
    # Geometric mean
    if (islog) {
      # Data is already on log scale - use directly
      log_stats <- aggregate(observed, list(level = focus_var), calc_group_stats)
      base_log_mean <- mean(log_stats$x[, 1])
      agg_df <- data.frame(
        level = log_stats$level,
        unstan = exp(log_stats$x[, 1] - base_log_mean),
        unstan_se = NA_real_, # SE not meaningful on linear scale for log data
        unstan_cv = log_stats$x[, 2] # CV is SE in log space
      )
    } else {
      # Data is on original scale - need to log first
      if (any(observed <= 0)) {
        warning("Non-positive values detected for geometric mean. Using arithmetic mean instead.")
        return(calculate_index_by_mean_type(observed, focus_var, "arithmetic", islog, preserve_probability_scale))
      }
      log_stats <- aggregate(log(observed), list(level = focus_var), calc_group_stats)
      base_log_mean <- mean(log_stats$x[, 1])
      agg_df <- data.frame(
        level = log_stats$level,
        unstan = exp(log_stats$x[, 1] - base_log_mean),
        unstan_se = NA_real_, # SE not meaningful on linear scale for log data
        unstan_cv = log_stats$x[, 2] # CV is SE in log space
      )
    }
  }

  return(agg_df)
}

#' @title Calculate Unstandardised Index
#' @description Calculate the unstandardised index with robust zero handling and family-specific methods
#' @param observed Numeric vector of observed values
#' @param focus_var Factor or numeric vector of focus variable levels
#' @param islog Logical. Whether to use log transformation
#' @param family Character. GLM family method to use
#' @return Data frame with level and unstan columns
#' @noRd
calculate_unstandardised_index <- function(observed, focus_var, islog = NULL, family = "gaussian", preserve_probability_scale = FALSE) {
  if (is.null(islog)) {
    islog <- all(observed > 0) && !any(observed == 0)
  }

  # Determine the appropriate mean type based on family and islog
  if (family == "binomial") {
    # For binomial models, use specialized handling for proportions
    mean_type <- ifelse(preserve_probability_scale, "raw", "arithmetic")
  } else if (family %in% c("gamma", "poisson", "weibull")) {
    # For gamma, poisson, and weibull models, use geometric mean (appropriate for positive/count data)
    mean_type <- "geometric"
  } else if (family == "gaussian") {
    # For gaussian models, respect islog parameter
    mean_type <- ifelse(islog, "geometric", "arithmetic")
  } else {
    # Default fallback for other families (quasi, etc.) - respect islog parameter
    mean_type <- ifelse(islog, "geometric", "arithmetic")
  }

  # Handle binomial family special cases
  if (family == "binomial") {
    if (all(observed %in% c(0, 1))) {
      # Binary data - calculate proportions with binomial SE
      prop_stats <- aggregate(observed, list(level = focus_var), function(x) {
        prop <- sum(x) / length(x)
        n <- length(x)
        se <- ifelse(n > 1, sqrt(prop * (1 - prop) / n), NA_real_)
        cv <- ifelse(prop > 0 && prop < 1, se / prop, NA_real_)
        c(prop, se, cv)
      })

      if (preserve_probability_scale) {
        # For raw rescaling, preserve actual proportions (probability scale)
        agg_df <- data.frame(
          level = prop_stats$level,
          unstan = prop_stats$x[, 1], # Keep actual proportions
          unstan_se = prop_stats$x[, 2], # Keep actual SEs
          unstan_cv = prop_stats$x[, 3]
        )
      } else {
        # Convert to relative scale for other rescaling methods
        overall_prop <- mean(observed)
        base_prop <- ifelse(overall_prop > 0 && overall_prop < 1, overall_prop, mean(prop_stats$x[, 1]))

        agg_df <- data.frame(
          level = prop_stats$level,
          unstan = prop_stats$x[, 1] / base_prop,
          unstan_se = prop_stats$x[, 2] / base_prop,
          unstan_cv = prop_stats$x[, 3]
        )
      }
    } else {
      # Already proportions or continuous data between 0-1
      return(calculate_index_by_mean_type(observed, focus_var, mean_type, islog, preserve_probability_scale))
    }
  } else {
    # Use the standardized helper function for all other families
    agg_df <- calculate_index_by_mean_type(observed, focus_var, mean_type, islog, preserve_probability_scale)
  }

  return(agg_df)
}

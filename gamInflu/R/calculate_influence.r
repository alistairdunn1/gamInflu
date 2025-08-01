#' @title Perform Influence Calculations
#' @description Perform calculations on a \code{gam_influence} object, including unstandardised and standardised indices, step-wise model building, and influence statistics. Supports subset-based analysis and multiple GLM families.
#' @param obj A \code{gam_influence} object.
#' @param islog Logical. Is the response variable log-transformed? If NULL (default), inferred from response name or model family.
#' @param rescale_method Character. How to rescale the indices. Options: "geometric_mean" (default), "arithmetic_mean", "raw", "custom".
#' @param custom_rescale_value Numeric. Custom rescaling value when rescale_method = "custom".
#' @param confidence_level Numeric. Confidence level for standardised index intervals (default 0.95).
#' @param family_method Character. How to handle different GLM families. Options: "auto" (default), "gaussian", "binomial", "gamma", "poisson".
#' @param use_coeff_method Logical. If TRUE (default), uses coefficient-based CI calculation (influ.r approach). If FALSE, uses prediction-based CI calculation.
#' @param subset_var Character. Name of variable to subset on (e.g., "area", "gear_type").
#' @param subset_value Value to subset by (e.g., "North", "Trawl"). If both subset_var and subset_value are provided, analysis is performed only on the subset, using the full model for predictions.
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
                                rescale_method = c("geometric_mean", "arithmetic_mean", "raw", "custom"),
                                custom_rescale_value = 1,
                                confidence_level = 0.95,
                                family_method = c("auto", "gaussian", "binomial", "gamma", "poisson"),
                                subset_var = NULL,
                                subset_value = NULL, ...) {
  UseMethod("calculate_influence")
}

#' @rdname calculate_influence
#' @method calculate_influence gam_influence
#' @importFrom stats predict aggregate logLik AIC deviance terms update as.formula var cov setNames qnorm
#' @export
calculate_influence.gam_influence <- function(obj, islog = NULL,
                                              rescale_method = c("geometric_mean", "arithmetic_mean", "raw", "custom"),
                                              custom_rescale_value = 1,
                                              confidence_level = 0.95,
                                              family_method = c("auto", "gaussian", "binomial", "gamma", "poisson"),
                                              subset_var = NULL,
                                              subset_value = NULL, ...) {
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
  family_method <- match.arg(family_method)
  rescale_method <- match.arg(rescale_method)

  # Detect model family
  if (family_method == "auto") {
    model_family <- obj$model$family$family
    model_link <- obj$model$family$link

    # Map common families, including pattern matching for Tweedie and NB
    if (grepl("^Tweedie\\(", model_family)) {
      family_detected <- "gamma" # Treat Tweedie like Gamma for calculations
    } else if (grepl("^Negative Binomial\\(", model_family)) {
      family_detected <- "poisson" # Treat NB like Poisson for calculations
    } else {
      family_detected <- switch(model_family,
        "gaussian" = "gaussian",
        "binomial" = "binomial",
        "Gamma" = "gamma",
        "gamma" = "gamma",
        "poisson" = "poisson",
        "quasi" = "gaussian", # Default for quasi families
        "gaussian" # Default fallback
      )
    }
    family_method <- family_detected
    message("Detected model family: ", model_family, " with link: ", model_link)
    message("Using family_method: ", family_method)
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
      # Log link (e.g., Gamma(link="log")) ≠ Pre-logged response (log(y) ~ ...)
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

  # Always ensure obj$islog is set for plotting functions
  obj$islog <- islog

  observed <- obj$data[[obj$response]]

  # --- 1. Unstandardised Index ---
  # Calculate the raw, unadjusted index for the focus term.
  # Uses family-appropriate methods for different GLM families.
  indices_df <- calculate_unstandardised_index(observed, obj$data[[obj$focus]], islog, family_method)

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
    # --- COEFFICIENT-BASED APPROACH (influ.r method) ---
    message("Using coefficient-based CI calculation (influ.r approach)")

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
    focus_levels <- levels(obj$data[[obj$focus]])
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

        # Calculate variance matrix for relative effects: Q * Σ * Q'
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

    stan_df$standardised_index <- exp(relative_coeffs)
    stan_df$stan_lower <- exp(relative_coeffs - ci_multiplier * focus_ses)
    stan_df$stan_upper <- exp(relative_coeffs + ci_multiplier * focus_ses)

    # Calculate SE and CV on response scale (approximation)
    stan_df$stan_se <- focus_ses * stan_df$standardised_index # Delta method approximation
    stan_df$standardised_cv <- focus_ses # CV ≈ SE on log scale for coefficient method
  } else {
    # --- PREDICTION-BASED APPROACH (modern method) ---
    message("Using prediction-based CI calculation (modern approach)")

    stan_df$standardised_index <- stan_df$pred / base_pred

    # Calculate confidence intervals on the response scale
    alpha <- 1 - confidence_level
    z_score <- qnorm(1 - alpha / 2)

    # Calculate confidence intervals directly on response scale
    stan_df$stan_lower <- (stan_df$pred - z_score * stan_df$se) / base_pred
    stan_df$stan_upper <- (stan_df$pred + z_score * stan_df$se) / base_pred

    # Ensure confidence intervals are non-negative for positive-valued responses
    if (all(stan_df$pred > 0)) {
      stan_df$stan_lower <- pmax(stan_df$stan_lower, 0)
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

      # Delta method CV: sqrt(exp(σ²) - 1) where σ is SE on log scale
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
  } # Apply rescaling to both unstandardised and standardised indices
  indices_df$unstan <- rescale_index(indices_df$unstan, rescale_method, custom_rescale_value)

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

  # Define all_terms for use in both stepwise and influence calculations
  all_terms <- obj$terms

  # Check if this is subset analysis
  is_subset_analysis <- !is.null(subset_var) && !is.null(subset_value)

  if (!is_subset_analysis) {
    # Perform stepwise analysis only for full dataset analysis

    # Start with an intercept-only model
    model_int <- update(obj$model, . ~ 1, data = obj$data)
    logLike_int <- as.numeric(logLik(model_int))
    summary_list[["intercept"]] <- data.frame(
      term = "Intercept",
      logLike = logLike_int,
      aic = AIC(model_int),
      r_sq = 0,
      deviance_explained = 0
    )

    # Sequentially add each term from the original model formula
    for (i in seq_along(all_terms)) {
      current_terms <- all_terms[1:i]
      formula_str <- paste("~", paste(current_terms, collapse = " + "))
      model_step <- update(obj$model, as.formula(formula_str), data = obj$data)

      # Store the step-wise index if the focus term is in the current model
      if (obj$focus %in% current_terms) {
        step_preds <- predict(model_step, newdata = obj$data, type = "terms")
        focus_cols <- grep(paste0("^", obj$focus, ""), colnames(step_preds), value = TRUE)
        if (length(focus_cols) == 0) {
          stop(paste0("Could not find any columns for focus term '", obj$focus, "' in step-wise model predictions (type='terms')."))
        }
        focus_effect_sum <- rowSums(step_preds[, focus_cols, drop = FALSE])
        step_focus_effects <- aggregate(focus_effect_sum ~ obj$data[[obj$focus]], FUN = mean)
        names(step_focus_effects) <- c("level", "effect")

        base_step <- mean(step_focus_effects$effect)
        step_focus_effects$index <- exp(step_focus_effects$effect - base_step)

        # Name the column after the term that was just added
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

      summary_list[[all_terms[i]]] <- data.frame(
        term = all_terms[i],
        logLike = logLike_step,
        aic = AIC(model_step),
        r_sq = r_sq_step,
        deviance_explained = dev_expl_step
      )
    }
  } else {
    # For subset analysis, provide minimal summary stats
    message("Stepwise analysis skipped for subset analysis")
    summary_list[["full_model"]] <- data.frame(
      term = "Full Model (Subset Analysis)",
      logLike = as.numeric(logLik(obj$model)),
      aic = AIC(obj$model),
      r_sq = summary(obj$model)$r.sq,
      deviance_explained = summary(obj$model)$dev.expl
    )
  }

  summary_df <- do.call(rbind, summary_list)
  # Calculate successive differences
  summary_df$r_sq_diff <- c(0, diff(summary_df$r_sq))
  summary_df$deviance_explained_diff <- c(0, diff(summary_df$deviance_explained))
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
  non_focus_terms <- setdiff(all_terms, obj$focus)

  for (term in non_focus_terms) {
    # Use the robust helper to find all relevant columns
    term_cols <- find_term_columns(term, colnames(all_preds_df))
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
            matches <- grep(term_pattern, terms_in_formula, value = TRUE, fixed = TRUE)
            if (length(matches) > 0) matches[1] else col
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
    supported_families <- c("gaussian", "binomial", "Gamma", "gamma", "poisson", "quasi", "quasipoisson", "quasibinomial", "Tweedie", "nb")

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
find_term_columns <- function(term, colnames_vec) {
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
rescale_index <- function(index, method = c("geometric_mean", "arithmetic_mean", "raw", "custom"),
                          custom_value = 1) {
  method <- match.arg(method)
  switch(method,
    "geometric_mean" = index / geometric_mean(index),
    "arithmetic_mean" = index / mean(index),
    "raw" = index,
    "custom" = index / geometric_mean(index) * custom_value
  )
}

#' @title Calculate Unstandardised Index
#' @description Calculate the unstandardised index with robust zero handling and family-specific methods
#' @param observed Numeric vector of observed values
#' @param focus_var Factor or numeric vector of focus variable levels
#' @param islog Logical. Whether to use log transformation
#' @param family_method Character. GLM family method to use
#' @return Data frame with level and unstan columns
#' @noRd
calculate_unstandardised_index <- function(observed, focus_var, islog = NULL, family_method = "gaussian") {
  if (is.null(islog)) {
    islog <- all(observed > 0) && !any(observed == 0)
  }

  # Helper function to calculate statistics for a group
  calc_group_stats <- function(x) {
    n <- length(x)
    if (n <= 1) {
      return(c(mean = mean(x), se = NA_real_, cv = NA_real_))
    }
    mean_val <- mean(x)
    sd_val <- sd(x)
    se_val <- sd_val / sqrt(n)
    cv_val <- ifelse(mean_val != 0, sd_val / abs(mean_val), NA_real_)
    return(c(mean = mean_val, se = se_val, cv = cv_val))
  }

  # Family-specific index calculations
  if (family_method == "binomial") {
    # For binomial models, work with proportions
    if (all(observed %in% c(0, 1))) {
      # Binary data - calculate proportions with binomial SE
      prop_stats <- aggregate(observed, list(level = focus_var), function(x) {
        prop <- sum(x) / length(x)
        n <- length(x)
        se <- ifelse(n > 1, sqrt(prop * (1 - prop) / n), NA_real_)
        cv <- ifelse(prop > 0 && prop < 1, se / prop, NA_real_)
        c(prop, se, cv)
      })

      # Convert to relative scale
      overall_prop <- mean(observed)
      base_prop <- ifelse(overall_prop > 0 && overall_prop < 1, overall_prop, mean(prop_stats$x[, 1]))

      agg_df <- data.frame(
        level = prop_stats$level,
        unstan = prop_stats$x[, 1] / base_prop,
        unstan_se = prop_stats$x[, 2] / base_prop,
        unstan_cv = prop_stats$x[, 3]
      )
    } else {
      # Already proportions or continuous data between 0-1
      stats_agg <- aggregate(observed, list(level = focus_var), calc_group_stats)
      base_mean <- mean(stats_agg$x[, 1])

      agg_df <- data.frame(
        level = stats_agg$level,
        unstan = stats_agg$x[, 1] / base_mean,
        unstan_se = stats_agg$x[, 2] / base_mean,
        unstan_cv = stats_agg$x[, 3]
      )
    }
  } else if (family_method == "gamma") {
    # For gamma models, always use geometric mean (positive data expected)
    if (any(observed <= 0)) {
      warning("Gamma family expects positive values. Using arithmetic mean for non-positive data.")
      stats_agg <- aggregate(observed, list(level = focus_var), calc_group_stats)
      base_mean <- mean(stats_agg$x[, 1])

      agg_df <- data.frame(
        level = stats_agg$level,
        unstan = stats_agg$x[, 1] / base_mean,
        unstan_se = stats_agg$x[, 2] / base_mean,
        unstan_cv = stats_agg$x[, 3]
      )
    } else {
      # Use geometric mean for positive gamma data - work in log space
      log_stats <- aggregate(log(observed), list(level = focus_var), calc_group_stats)
      base_log_mean <- mean(log_stats$x[, 1])

      agg_df <- data.frame(
        level = log_stats$level,
        unstan = exp(log_stats$x[, 1] - base_log_mean),
        unstan_se = NA_real_, # SE not meaningful on linear scale for log data
        unstan_cv = log_stats$x[, 2] # CV is SE in log space
      )
    }
  } else if (family_method == "poisson") {
    # For Poisson models, use geometric mean if positive, arithmetic if zeros present
    if (any(observed < 0)) {
      warning("Poisson family expects non-negative values.")
    }

    if (all(observed > 0)) {
      # Use geometric mean for positive count data - work in log space
      log_stats <- aggregate(log(observed), list(level = focus_var), calc_group_stats)
      base_log_mean <- mean(log_stats$x[, 1])

      agg_df <- data.frame(
        level = log_stats$level,
        unstan = exp(log_stats$x[, 1] - base_log_mean),
        unstan_se = NA_real_, # SE not meaningful on linear scale for log data
        unstan_cv = log_stats$x[, 2] # CV is SE in log space
      )
    } else {
      # Use arithmetic mean when zeros present
      stats_agg <- aggregate(observed, list(level = focus_var), calc_group_stats)
      base_mean <- mean(stats_agg$x[, 1])

      agg_df <- data.frame(
        level = stats_agg$level,
        unstan = stats_agg$x[, 1] / base_mean,
        unstan_se = stats_agg$x[, 2] / base_mean,
        unstan_cv = stats_agg$x[, 3]
      )
    }
  } else {
    # Gaussian or other families - original logic
    if (islog && all(observed > 0)) {
      # Use geometric mean for positive data - work in log space
      log_stats <- aggregate(log(observed), list(level = focus_var), calc_group_stats)
      base_log_mean <- mean(log_stats$x[, 1])

      agg_df <- data.frame(
        level = log_stats$level,
        unstan = exp(log_stats$x[, 1] - base_log_mean),
        unstan_se = NA_real_, # SE not meaningful on linear scale for log data
        unstan_cv = log_stats$x[, 2] # CV is SE in log space
      )
    } else {
      # Fallback to arithmetic mean for data with zeros or non-log
      stats_agg <- aggregate(observed, list(level = focus_var), calc_group_stats)
      base_mean <- mean(stats_agg$x[, 1])

      agg_df <- data.frame(
        level = stats_agg$level,
        unstan = stats_agg$x[, 1] / base_mean,
        unstan_se = stats_agg$x[, 2] / base_mean,
        unstan_cv = stats_agg$x[, 3]
      )
    }
  }

  # Ensure all required columns exist
  if (!"unstan_se" %in% names(agg_df)) agg_df$unstan_se <- NA_real_
  if (!"unstan_cv" %in% names(agg_df)) agg_df$unstan_cv <- NA_real_

  return(agg_df)
}

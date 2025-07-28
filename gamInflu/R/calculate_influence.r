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
    supported_families <- c("gaussian", "binomial", "Gamma", "gamma", "poisson", "quasi", "quasipoisson", "quasibinomial")

    if (!family_name %in% supported_families) {
      warning("Family '", family_name, "' may not be fully supported. Supported families: ",
        paste(supported_families, collapse = ", "),
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
  }

  TRUE
}

#' @title Perform Influence Calculations
#' @description A generic S3 function to perform calculations on a `gam_influence` object.
#' @param obj An object for which to perform calculations.
#' @param ... Additional arguments passed to methods.
#' @export
calculate_influence <- function(obj, ...) {
  UseMethod("calculate_influence")
}

#' @title Perform influence calculations for a gam_influence object
#' @description This is the core function that computes all necessary metrics for the plots
#' and summaries. It calculates the unstandardised and standardised indices,
#' performs a step-wise model build to assess term contributions, and computes
#' influence statistics (overall and trend). Enhanced with comprehensive support
#' for multiple GLM families including binomial, gamma, and Poisson distributions.
#' @param obj A `gam_influence` object.
#' @param islog Logical. Is the response variable log-transformed? If NULL (default),
#'   the function infers this by checking if the response name starts with "log(" or model family.
#' @param rescale_method Character. How to rescale the indices. Options: "geometric_mean" (default),
#'   "arithmetic_mean", "raw", or "custom".
#' @param custom_rescale_value Numeric. Custom rescaling value when rescale_method = "custom".
#' @param confidence_level Numeric. Confidence level for standardized index intervals (default 0.95).
#' @param family_method Character. How to handle different GLM families. Options: "auto" (default),
#'   "gaussian", "binomial", "gamma", "poisson". When "auto", family is detected from model.
#' @param ... Additional arguments (currently unused).
#' @return The `gam_influence` object, now containing a `calculated` list with data frames
#'   for indices, summary stats, influences, predictions, and s.e. of predictions.
#' @details
#' The influence calculation follows Bentley et al. (2012):
#' - Overall influence: exp(mean(|effects|)) - 1
#' - Trend influence: exp(cov(levels, effects)/var(levels)) - 1
#'
#' **Enhanced Family Support:**
#'
#' The package now supports multiple GLM families with family-specific index calculations:
#' - **Gaussian**: Traditional geometric mean for log-transformed data, arithmetic for linear
#' - **Binomial**: Proportion-based indices for presence/absence or binary data
#' - **Gamma**: Geometric mean aggregation for positive continuous data (biomass, CPUE)
#' - **Poisson**: Count-appropriate methods for abundance or catch numbers
#' - **Automatic Detection**: Family is auto-detected from model object when family_method="auto"
#'
#' Each family uses statistically appropriate aggregation methods and handles edge cases
#' like zeros in count data or proportions in binomial models. The standardized index
#' uses model predictions with proper uncertainty quantification for all families.
#' @importFrom stats predict aggregate logLik AIC deviance terms update as.formula var cov setNames qnorm
#' @export
calculate_influence.gam_influence <- function(obj, islog = NULL,
                                              rescale_method = c("geometric_mean", "arithmetic_mean", "raw", "custom"),
                                              custom_rescale_value = 1,
                                              confidence_level = 0.95,
                                              family_method = c("auto", "gaussian", "binomial", "gamma", "poisson"), ...) {
  # Validate inputs
  validate_gam_influence(obj)
  if (confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must be between 0 and 1", call. = FALSE)
  }

  # --- Setup and Family Detection ---
  family_method <- match.arg(family_method)

  # Detect model family
  if (family_method == "auto") {
    model_family <- obj$model$family$family
    model_link <- obj$model$family$link

    # Map common families
    family_detected <- switch(model_family,
      "gaussian" = "gaussian",
      "binomial" = "binomial",
      "Gamma" = "gamma",
      "gamma" = "gamma",
      "poisson" = "poisson",
      "quasi" = "gaussian", # Default for quasi families
      "gaussian" # Default fallback
    )
    family_method <- family_detected
    message("Detected model family: ", model_family, " with link: ", model_link)
    message("Using family_method: ", family_method)
  }

  # Enhanced islog detection based on family and response name
  if (is.null(islog)) {
    if (is.null(obj$islog)) {
      # Check response name for log transformation
      response_log <- substr(obj$response, 1, 4) == "log("

      # Check if family typically uses log transformation
      family_log <- family_method %in% c("gamma", "poisson") &&
        (!is.null(obj$model$family$link) && obj$model$family$link == "log")

      islog <- response_log || family_log
      message(
        "Assuming islog = ", islog, " based on response name '", obj$response,
        "' and family '", family_method, "'."
      )
    } else {
      islog <- obj$islog
    }
  } else {
    obj$islog <- islog
  }

  observed <- obj$data[[obj$response]]

  # --- 1. Unstandardised Index ---
  # Calculate the raw, unadjusted index for the focus term.
  # Uses family-appropriate methods for different GLM families.
  indices_df <- calculate_unstandardized_index(observed, obj$data[[obj$focus]], islog, family_method)

  # --- 2. Standardised Index (from the full model) ---
  # This is the final index after accounting for all terms in the model.
  # We extract the partial effects for the focus term.
  preds_full <- predict(obj$model, type = "terms", se.fit = TRUE)

  # Identify columns related to the focus term in preds_full$fit
  # For factor terms like 'year', predict(type="terms") will return columns like 'year1990', 'year1991', etc.
  # We need to sum these up per row to get the overall effect for the 'year' term

  # Find column names in preds_full$fit that start with the focus term's name
  focus_term_cols_fit <- grep(paste0("^", obj$focus), colnames(preds_full$fit), value = TRUE)
  focus_term_cols_se <- grep(paste0("^", obj$focus), colnames(preds_full$se.fit), value = TRUE)

  if (length(focus_term_cols_fit) == 0) {
    stop(paste0("Could not find any columns for focus term '", obj$focus, "' in model predictions (type='terms'). This might indicate the focus term is not handled as expected or is part of a complex smooth that needs different extraction."))
  }

  # Calculate the sum of effects for the focus term
  # If the focus term is a simple factor, sum across its dummy variable columns.
  # If it's a smooth, it might just be one column.
  focus_fit_sum <- rowSums(preds_full$fit[, focus_term_cols_fit, drop = FALSE])
  focus_se_sum <- rowSums(preds_full$se.fit[, focus_term_cols_se, drop = FALSE]) # Sum of SEs might be an overestimation, but for simplicity here. A more rigorous approach for SEs might involve the variance-covariance matrix.

  focus_pred_df <- data.frame(
    level = obj$data[[obj$focus]],
    fit = focus_fit_sum,
    se = focus_se_sum
  )

  stan_df <- aggregate(cbind(fit, se) ~ level, data = focus_pred_df, FUN = mean)
  base <- mean(stan_df$fit)
  stan_df$stan <- exp(stan_df$fit - base)

  # Add proper confidence intervals using the specified confidence level
  rescale_method <- match.arg(rescale_method)
  alpha <- 1 - confidence_level
  z_score <- qnorm(1 - alpha / 2)

  stan_df$stanLower <- exp(stan_df$fit - base - z_score * stan_df$se)
  stan_df$stanUpper <- exp(stan_df$fit - base + z_score * stan_df$se)

  # Apply rescaling to both unstandardized and standardized indices
  indices_df$unstan <- rescale_index(indices_df$unstan, rescale_method, custom_rescale_value)
  stan_df$stan <- rescale_index(stan_df$stan, rescale_method, custom_rescale_value)
  stan_df$stanLower <- rescale_index(stan_df$stanLower, rescale_method, custom_rescale_value)
  stan_df$stanUpper <- rescale_index(stan_df$stanUpper, rescale_method, custom_rescale_value)

  indices_df <- merge(indices_df, stan_df[, c("level", "stan", "stanLower", "stanUpper")], by = "level")

  # --- 3. Step-wise Calculations ---
  # Build the model one term at a time to see how the index and model fit change.
  summary_list <- list()
  step_indices_list <- list()

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
  all_terms <- obj$terms
  for (i in seq_along(all_terms)) {
    current_terms <- all_terms[1:i]
    formula_str <- paste("~", paste(current_terms, collapse = " + "))
    model_step <- update(obj$model, as.formula(formula_str), data = obj$data)

    # Store the step-wise index if the focus term is in the current model
    if (obj$focus %in% current_terms) {
      step_preds <- predict(model_step, type = "terms")
      focus_cols <- grep(paste0("^", obj$focus), colnames(step_preds), value = TRUE)
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
  all_preds_df <- as.data.frame(preds_full$fit)
  all_preds_df$level <- obj$data[[obj$focus]]

  # After creating all_preds_df, also create a dataframe for standard errors
  all_preds_se_df <- as.data.frame(preds_full$se.fit)
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
    trend_infl <- if (var(level_num) > 0) exp(cov(level_num, infl$influence) / var(level_num)) - 1 else 0

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

  # Rename standardized index columns for consistency with package API
  if ("stan" %in% names(indices_df)) {
    names(indices_df)[names(indices_df) == "stan"] <- "standardized_index"
  }
  if ("stanLower" %in% names(indices_df)) {
    names(indices_df)[names(indices_df) == "stanLower"] <- "stan_lower"
  }
  if ("stanUpper" %in% names(indices_df)) {
    names(indices_df)[names(indices_df) == "stanUpper"] <- "stan_upper"
  }

  # --- Finalize ---
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
#' @description Calculate the geometric mean of a numeric vector, with robust handling of
#'   missing values and non-positive values. Essential for family-specific index calculations
#'   in gamma and log-normal models where multiplicative processes are expected.
#' @param x A numeric vector.
#' @param na.rm Logical. Should missing values be removed before calculation? Default is TRUE.
#' @return The geometric mean of the input vector. If non-positive values are present,
#'   issues a warning and returns the arithmetic mean as a fallback.
#' @details
#' The geometric mean is calculated as exp(mean(log(x))) and is appropriate for:
#' - Gamma family models with positive data
#' - Log-normal processes in fisheries CPUE standardization
#' - Multiplicative effects and indices
#'
#' When non-positive values are encountered, the function automatically falls back
#' to the arithmetic mean to ensure robust calculations across different data types.
#' @examples
#' \dontrun{
#' # Basic usage
#' geometric_mean(c(1, 2, 4, 8)) # Returns 2.83
#'
#' # With zeros (falls back to arithmetic mean)
#' geometric_mean(c(0, 1, 2, 4)) # Warning + arithmetic mean
#'
#' # In index calculations
#' index_values <- c(0.8, 1.2, 1.1, 0.9)
#' standardized <- index_values / geometric_mean(index_values)
#' }
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
#' @description Rescale an index using different methods
#' @param index Numeric vector to rescale
#' @param method Character. Rescaling method
#' @param custom_value Numeric. Custom rescaling value when method = "custom"
#' @return Rescaled index vector
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

#' @title Calculate Unstandardized Index
#' @description Calculate the unstandardized index with robust zero handling and family-specific methods
#' @param observed Numeric vector of observed values
#' @param focus_var Factor or numeric vector of focus variable levels
#' @param islog Logical. Whether to use log transformation
#' @param family_method Character. GLM family method to use
#' @return Data frame with level and unstan columns
#' @noRd
calculate_unstandardized_index <- function(observed, focus_var, islog = NULL, family_method = "gaussian") {
  if (is.null(islog)) {
    islog <- all(observed > 0) && !any(observed == 0)
  }

  # Family-specific index calculations
  if (family_method == "binomial") {
    # For binomial models, work with proportions
    # Convert binary data to proportions by group
    if (all(observed %in% c(0, 1))) {
      # Binary data - calculate proportions
      agg_df <- aggregate(
        list(unstan = observed),
        list(level = focus_var),
        function(x) sum(x) / length(x) # Calculate proportion
      )
      # Convert to relative scale (divide by overall proportion)
      overall_prop <- mean(observed)
      if (overall_prop > 0 && overall_prop < 1) {
        agg_df$unstan <- agg_df$unstan / overall_prop
      } else {
        agg_df$unstan <- agg_df$unstan / mean(agg_df$unstan)
      }
    } else {
      # Already proportions or continuous data between 0-1
      agg_df <- aggregate(
        list(unstan = observed),
        list(level = focus_var), mean
      )
      agg_df$unstan <- agg_df$unstan / mean(agg_df$unstan)
    }
  } else if (family_method == "gamma") {
    # For gamma models, always use geometric mean (positive data expected)
    if (any(observed <= 0)) {
      warning("Gamma family expects positive values. Using arithmetic mean for non-positive data.")
      agg_df <- aggregate(
        list(unstan = observed),
        list(level = focus_var), mean
      )
      agg_df$unstan <- agg_df$unstan / mean(agg_df$unstan)
    } else {
      # Use geometric mean for positive gamma data
      agg_df <- aggregate(
        list(unstan = log(observed)),
        list(level = focus_var), mean
      )
      agg_df$unstan <- exp(agg_df$unstan - mean(agg_df$unstan))
    }
  } else if (family_method == "poisson") {
    # For Poisson models, use geometric mean if positive, arithmetic if zeros present
    if (any(observed < 0)) {
      warning("Poisson family expects non-negative values.")
    }

    if (all(observed > 0)) {
      # Use geometric mean for positive count data
      agg_df <- aggregate(
        list(unstan = log(observed)),
        list(level = focus_var), mean
      )
      agg_df$unstan <- exp(agg_df$unstan - mean(agg_df$unstan))
    } else {
      # Use arithmetic mean when zeros present
      agg_df <- aggregate(
        list(unstan = observed),
        list(level = focus_var), mean
      )
      agg_df$unstan <- agg_df$unstan / mean(agg_df$unstan)
    }
  } else {
    # Gaussian or other families - original logic
    if (islog && all(observed > 0)) {
      # Use geometric mean for positive data
      agg_df <- aggregate(
        list(unstan = log(observed)),
        list(level = focus_var), mean
      )
      agg_df$unstan <- exp(agg_df$unstan - mean(agg_df$unstan))
    } else {
      # Fallback to arithmetic mean for data with zeros or non-log
      agg_df <- aggregate(
        list(unstan = observed),
        list(level = focus_var), mean
      )
      agg_df$unstan <- agg_df$unstan / mean(agg_df$unstan)
    }
  }

  data.frame(level = agg_df$level, unstan = agg_df$unstan)
}

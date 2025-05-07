#' @title calculate_influence
#' @description Performs the core calculations for influence analysis on a pre-created
#' `influence_gam` object. This populates the `indices`, `summary`,
#' `influences`, and `preds` slots of the object.
#'
#' @author Alistair Dunn
#' @param obj An object of class `influence_gam` created by `create_influence_gam`.
#' @param islog Deprecated. Logic now tries to infer from the influence object or the model.
#' @param verbose Logical. Print progress messages?
#' @param ... Additional arguments (currently unused).
#'
#' @return The `influence_gam` object with calculation results added.
#' @export
calculate_influence <- function(obj, islog = NULL, verbose = TRUE, ...) {
  UseMethod("calculate_influence")
}

#' @rdname calculate_influence
#' @export
calculate_influence.influence_gam <- function(obj, islog = NULL, verbose = TRUE, ...) {
  if (!inherits(obj, "influence_gam")) {
    stop("Input 'obj' must be of class 'influence_gam'.")
  }

  model <- obj$model
  data <- obj$data
  response <- obj$response
  focus <- obj$focus
  focus_var <- obj$focus_var # Use the derived variable name
  all_terms <- obj$terms

  if (verbose) message("Starting influence calculations for focus: ", focus)

  # --- Preparations ---
  # Ensure focus variable is a factor in the data used for calculations
  if (!is.factor(data[[focus_var]])) {
    data[[focus_var]] <- factor(data[[focus_var]])
    if (verbose) message("Converted focus variable '", focus_var, "' to factor for calculations.")
  }
  focus_levels <- levels(data[[focus_var]])

  # Get observed values
  observed <- data[[response]]
  # Infer if response is logged (simple check)
  if (is.null(islog)) {
    if (!is.null(obj$islog)) {
      logged <- obj$islog
    } else {
      logged <- grepl("^log\\(.*\\)$", response) || inherits(model$family, "gaussian") && model$family$link == "log"
    }
  } else {
    logged <- islog # Allow override, though discouraged
  }
  # stop("Input 'obj' must be of class 'influence_gam'.")

  # --- 1. Unstandardised Index ---
  if (verbose) message("Calculating unstandardised index...")
  indices_df <- data.frame(level = focus_levels)

  # Use tryCatch for robustness if log is attempted on non-positive values
  try_log <- function(x) tryCatch(log(x), warning = function(w) NA, error = function(e) NA)

  if (logged || all(observed > 0, na.rm = TRUE)) {
    log_observed <- if (logged) observed else try_log(observed)
    if (anyNA(log_observed) && !logged) {
      if (verbose) warning("Non-positive values found in response '", response, "'. Cannot compute geometric mean for unstandardised index. Attempting arithmetic mean.", call. = FALSE)
      # Use stats::aggregate explicitly
      agg_unstan <- stats::aggregate(list(unstan = observed), list(level = data[[focus_var]]), mean, na.rm = TRUE)
      indices_df <- merge(indices_df, agg_unstan, by = "level", all.x = TRUE)
      indices_df$unstan <- indices_df$unstan / mean(indices_df$unstan, na.rm = TRUE)
    } else if (anyNA(log_observed) && logged) {
      if (verbose) warning("NA values found in logged response '", response, "'. Geometric mean might be affected.", call. = FALSE)
      agg_unstan <- stats::aggregate(list(unstan = log_observed), list(level = data[[focus_var]]), mean, na.rm = TRUE)
      indices_df <- merge(indices_df, agg_unstan, by = "level", all.x = TRUE)
      indices_df$unstan <- exp(indices_df$unstan - mean(indices_df$unstan, na.rm = TRUE))
    } else {
      # Standard geometric mean calculation
      agg_unstan <- stats::aggregate(list(unstan = log_observed), list(level = data[[focus_var]]), mean, na.rm = TRUE)
      indices_df <- merge(indices_df, agg_unstan, by = "level", all.x = TRUE)
      indices_df$unstan <- exp(indices_df$unstan - mean(indices_df$unstan, na.rm = TRUE))
    }
  } else {
    if (verbose) message("Response variable '", response, "' contains non-positive values and is not logged. Using arithmetic mean for unstandardised index.")
    agg_unstan <- stats::aggregate(list(unstan = observed), list(level = data[[focus_var]]), mean, na.rm = TRUE)
    indices_df <- merge(indices_df, agg_unstan, by = "level", all.x = TRUE)
    indices_df$unstan <- indices_df$unstan / mean(indices_df$unstan, na.rm = TRUE)
  }
  rownames(indices_df) <- indices_df$level # Easier access later

  # --- 2. Standardised Index (from full model focus term effect) ---
  if (verbose) message("Calculating standardised index...")
  # Predict term contributions for the full model
  # Need data with all levels of focus factor for prediction
  pred_data_focus <- data.frame(level = focus_levels)
  names(pred_data_focus) <- focus_var
  # Add average values for other predictors to predict marginal effect of focus
  other_vars <- setdiff(all.vars(stats::formula(model)[-2]), focus_var)
  for (v in other_vars) {
    if (v %in% names(data)) { # Check if variable exists in data
      if (is.numeric(data[[v]])) {
        pred_data_focus[[v]] <- mean(data[[v]], na.rm = TRUE)
      } else if (is.factor(data[[v]])) {
        # Ensure levels match original data's factor levels
        pred_data_focus[[v]] <- factor(levels(data[[v]])[1], levels = levels(data[[v]]))
      } else {
        # Handle other types if necessary (e.g., character, logical)
        pred_data_focus[[v]] <- data[[v]][1] # Just take the first value
      }
    } else {
      if (verbose) message("  Skipping variable '", v, "' for prediction grid (not found in data, possibly part of smooth spec).")
    }
  }


  # Predict on the link scale
  # Use stats::predict explicitly
  preds_link <- tryCatch(stats::predict(model, newdata = pred_data_focus, type = "link", se.fit = TRUE),
    error = function(e) {
      stop("Prediction failed for standardised index. Check model and prediction data. Error: ", e$message)
    }
  )

  # Calculate relative effect (exponentiated, centred)
  coeffs <- preds_link$fit
  ses <- preds_link$se.fit
  base <- mean(coeffs, na.rm = TRUE)

  # Use the inverse link function
  inv_link <- model$family$linkinv

  indices_df$stan <- inv_link(coeffs - base)
  indices_df$stanLower <- inv_link(coeffs - base - 2 * ses)
  indices_df$stanUpper <- inv_link(coeffs - base + 2 * ses)

  # --- 3. Step Plot Calculations ---
  if (verbose) message("Calculating step indices...")
  summary_list <- list()
  indices_step_list <- list()

  # Original model call components
  orig_call <- model$call
  gam_args <- as.list(orig_call)[-1] # Remove the function name
  gam_args$formula <- NULL # Will be built iteratively
  gam_args$data <- quote(data) # Use the data frame in this function's environment

  # Get full formula components
  full_formula <- stats::formula(model)
  response_str <- deparse(full_formula[[2]])
  terms_rhs <- attr(stats::terms(model), "term.labels")

  # Term 0: Intercept only
  if (verbose) message("  Step 0: Intercept only")
  formula_step0 <- stats::reformulate("1", response = response_str) # Use stats::reformulate
  model_step0 <- tryCatch(do.call(gam, c(list(formula = formula_step0), gam_args)),
    error = function(e) {
      warning("Failed to fit intercept-only model: ", e$message, call. = FALSE)
      NULL
    }
  )

  if (!is.null(model_step0)) {
    logLike0 <- tryCatch(stats::logLik(model_step0), error = function(e) NA) # Use stats::logLik
    aic0 <- tryCatch(stats::AIC(model_step0), error = function(e) NA) # Use stats::AIC
    k0 <- length(stats::coef(model_step0)) # Use stats::coef
    summary_list[["intercept"]] <- data.frame(term = "intercept", k = k0, logLike = logLike0, aic = aic0, r2 = 0, r2Dev = 0, r2Negel = 0)
  } else {
    logLike0 <- NA
    summary_list[["intercept"]] <- data.frame(term = "intercept", k = 1, logLike = NA, aic = NA, r2 = 0, r2Dev = 0, r2Negel = 0)
  }
  logLikeInterceptOnly <- logLike0 # Store for Negelkerke R2
  n_obs <- stats::nobs(model) # Use stats::nobs

  # Loop through terms, adding one at a time
  current_terms <- c() # Start with no terms on RHS
  for (i in 1:length(terms_rhs)) {
    term <- terms_rhs[i]
    current_terms <- c(current_terms, term)
    formula_step <- stats::reformulate(current_terms, response = response_str)
    step_name <- if (i == 1) term else paste("+", term)
    if (verbose) message("  Step ", i, ": Adding '", term, "' (", step_name, ")")

    model_step <- tryCatch(do.call(gam, c(list(formula = formula_step), gam_args)),
      error = function(e) {
        warning("Failed to fit model at step '", step_name, "': ", e$message, call. = FALSE)
        NULL
      }
    )

    if (!is.null(model_step)) {
      # Calculate index for this step (effect of focus term)
      # Use the same prediction approach as for the full model's standardised index
      preds_link_step <- tryCatch(stats::predict(model_step, newdata = pred_data_focus, type = "link", se.fit = FALSE),
        error = function(e) {
          warning("Prediction failed for step '", step_name, "'. Index set to NA. Error: ", e$message, call. = FALSE)
          rep(NA, nrow(pred_data_focus))
        }
      )

      if (all(is.na(preds_link_step))) {
        index_step <- rep(NA, length(focus_levels))
      } else {
        coeffs_step <- preds_link_step
        base_step <- mean(coeffs_step, na.rm = TRUE)
        index_step <- model_step$family$linkinv(coeffs_step - base_step) # Use inv link from this model
      }
      indices_step_list[[step_name]] <- index_step

      # Calculate summary statistics
      logLike_step <- tryCatch(stats::logLik(model_step), error = function(e) NA)
      aic_step <- tryCatch(stats::AIC(model_step), error = function(e) NA)
      k_step <- length(stats::coef(model_step))

      # R-squared variants
      model_summary <- summary(model_step) # Get summary once
      # Deviance R2
      null_deviance <- model_summary$null.deviance
      dev_step <- model_summary$deviance
      if (is.null(null_deviance) || is.null(dev_step) || is.na(null_deviance) || is.na(dev_step)) {
        if (verbose) message("  Warning: Null or Residual Deviance is NULL/NA for step ", step_name, ". Cannot calculate r2Dev.")
        r2Dev_step <- NA
      } else if (null_deviance < 1e-8) { # Avoid division by zero or near-zero
        if (verbose) message("  Warning: Null Deviance is near zero for step ", step_name, ". Setting r2Dev to NA.")
        r2Dev_step <- NA
      } else {
        r2Dev_step <- (null_deviance - dev_step) / null_deviance
      }

      # Negelkerke R2
      r2Negel_step <- NA
      if (!is.na(logLikeInterceptOnly) && !is.na(logLike_step) && !is.na(n_obs) && n_obs > 0) {
        ll0_scaled <- logLikeInterceptOnly * (2 / n_obs)
        llM_scaled <- logLike_step * (2 / n_obs)
        # Avoid potential issues with exp(very large negative number) -> 0
        if (is.finite(ll0_scaled) && is.finite(llM_scaled)) {
          denominator <- 1 - exp(ll0_scaled)
          if (abs(denominator) < .Machine$double.eps^0.5) { # Avoid division by zero if ll0 is large neg
            r2Negel_step <- NA
            if (verbose) message("  Warning: Denominator for Negelkerke R2 is near zero. Setting r2Negel to NA.")
          } else {
            r2Negel_step <- (1 - exp(llM_scaled - ll0_scaled)) / denominator
          }
        } else {
          r2Negel_step <- NA
          if (verbose) message("  Warning: Scaled log-likelihoods non-finite. Cannot calculate r2Negel.")
        }
      } else {
        if (verbose) message("  Warning: Missing values needed for Negelkerke R2 calculation.")
      }


      # Correlation-based R2 (use response scale predictions)
      fitted_vals <- tryCatch(stats::predict(model_step, type = "response"), error = function(e) NULL)
      r2_step <- if (is.null(fitted_vals)) NA else tryCatch(stats::cor(fitted_vals, observed)^2, error = function(e) NA) # Use stats::cor

      summary_list[[step_name]] <- data.frame(term = step_name, k = k_step, logLike = logLike_step, aic = aic_step, r2 = r2_step, r2Dev = r2Dev_step, r2Negel = r2Negel_step)
    } else {
      # Add placeholder if model failed
      summary_list[[step_name]] <- data.frame(term = step_name, k = NA, logLike = NA, aic = NA, r2 = NA, r2Dev = NA, r2Negel = NA)
      indices_step_list[[step_name]] <- rep(NA, length(focus_levels))
    }
  }
  # Combine step indices into the main indices data frame
  indices_df <- cbind(indices_df, do.call(cbind, indices_step_list))
  # Combine summary info
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL

  # --- 4. Influence Calculations ---
  if (verbose) message("Calculating term influences...")
  # Predict term contributions for the full model on the original data
  preds_terms_full <- tryCatch(stats::predict(model, type = "terms", se.fit = FALSE),
    error = function(e) {
      warning("Failed to predict terms for influence calculation: ", e$message, call. = FALSE)
      NULL
    }
  )

  influences_list <- list()
  if (!is.null(preds_terms_full)) {
    # Get names of terms from prediction matrix (handles smooths like s(x0), s(x1):facA etc.)
    pred_term_names <- colnames(preds_terms_full)

    for (term in pred_term_names) {
      # Skip the main effect term(s) associated with the focus variable
      # This logic is complex and heuristic based on naming conventions.
      is_interaction_involving_focus <- grepl(":", term, fixed = TRUE) && grepl(focus_var, term, fixed = TRUE)
      term_no_space <- gsub(" ", "", term)
      focus_var_no_space <- gsub(" ", "", focus_var)
      pattern_by_focus <- paste0("by=", focus_var_no_space)
      is_smooth_by_focus <- grepl(pattern_by_focus, term_no_space, fixed = TRUE)
      is_focus_main_effect <- term == focus_var # Simple case for non-smooth term

      # Skip if it's the focus main effect, UNLESS it's also part of an interaction or a 'by=' smooth
      if (is_focus_main_effect && !is_interaction_involving_focus && !is_smooth_by_focus) {
        if (verbose) message("  Skipping main focus effect term: ", term)
        next
      }
      # Log calculation for interaction/by terms involving focus or other terms
      if (verbose && (is_interaction_involving_focus || is_smooth_by_focus)) {
        message("  Calculating influence for term involving focus: ", term)
      } else if (verbose && !(grepl(focus_var, term, fixed = TRUE))) {
        message("  Calculating influence for term: ", term)
      }

      term_values <- preds_terms_full[, term]

      # Aggregate the effect of this term across levels of the focus variable
      if (is.numeric(term_values)) {
        agg_influ <- stats::aggregate(list(influence = term_values), list(level = data[[focus_var]]), mean, na.rm = TRUE)
        influences_list[[term]] <- agg_influ$influence
      } else {
        if (verbose) warning("Term '", term, "' predictions are not numeric. Skipping influence calculation for this term.", call. = FALSE)
        influences_list[[term]] <- rep(NA, length(focus_levels))
      }
    }
  } else {
    if (verbose) message("  Skipping influence calculation due to term prediction failure.")
  }

  # Combine influences
  influences_df <- data.frame(level = focus_levels)
  if (length(influences_list) > 0) {
    influences_df <- cbind(influences_df, do.call(cbind, influences_list))
  }
  rownames(influences_df) <- influences_df$level

  # --- Store Results ---
  obj$indices <- indices_df
  obj$summary <- summary_df
  obj$influences <- influences_df
  obj$preds <- if (is.null(preds_terms_full)) NULL else as.data.frame(preds_terms_full) # Store term predictions for CDI plot
  obj$calculated <- TRUE

  if (verbose) message("Influence calculations complete.")
  return(obj)
}

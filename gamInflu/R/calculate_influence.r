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
#' influence statistics (overall and trend).
#' @param obj A `gam_influence` object.
#' @param islog Logical. Is the response variable log-transformed? If NULL (default),
#'   the function infers this by checking if the response name starts with "log(".
#' @param ... Additional arguments (currently unused).
#' @return The `gam_influence` object, now containing a `calculated` list with data frames
#'   for indices, summary stats, influences, predictions, and s.e. of predictions.
#' @importFrom stats predict aggregate logLik AIC deviance terms update as.formula var cov setNames
#' @export
calculate_influence.gam_influence <- function(obj, islog = NULL, ...) {
  # --- Setup ---
  if (is.null(islog)) {
    if (is.null(obj$islog)) {
      islog <- substr(obj$response, 1, 4) == "log("
      message("Assuming islog = ", islog, " based on response name '", obj$response, "'.")
    } else {
      islog <- obj$islog
    }
  } else {
    obj$islog <- islog
  }

  observed <- obj$data[[obj$response]]

  # --- 1. Unstandardised Index ---
  # Calculate the raw, unadjusted index for the focus term.
  # Uses geometric mean for log-transformed data, arithmetic mean otherwise.
  if (islog) {
    agg_df <- aggregate(list(unstan = observed), list(level = obj$data[[obj$focus]]), mean)
    agg_df$unstan <- exp(agg_df$unstan - mean(agg_df$unstan))
  } else {
    agg_df <- aggregate(list(unstan = observed), list(level = obj$data[[obj$focus]]), mean)
    agg_df$unstan <- agg_df$unstan / mean(agg_df$unstan)
  }
  indices_df <- data.frame(level = agg_df$level, unstan = agg_df$unstan)

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
  stan_df$stanLower <- exp(stan_df$fit - base - 2 * stan_df$se)
  stan_df$stanUpper <- exp(stan_df$fit - base + 2 * stan_df$se)

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
  return(cols)
}

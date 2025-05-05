# Revised Influence Code for mgcv::gam and ggplot2
# Based on original Influ.r by Nokome Bentley, Trophia Ltd
# Modifications by Gemini (2025)

# Ensure required packages are loaded
# install.packages(c("mgcv", "ggplot2", "dplyr", "tidyr", "rlang"))
library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rlang) # For := assignment in dplyr/ggplot

# --- S3 Generic Functions ---

#' Calculate Influence Metrics for a GAM
#'
#' Performs the core calculations for influence analysis on a GAM model.
#' This populates the object with necessary data frames for plotting and summary.
#'
#' @param x An object of class 'influence_gam'.
#' @param islog Logical, specify if the response variable was log-transformed
#'        before modelling. If NULL (default), attempts to infer from the
#'        response variable name in the formula.
#' @param ... Additional arguments (currently unused).
#' @return The 'influence_gam' object, updated with calculation results
#'         (invisible).
#' @export
calculate <- function(x, ...) {
  UseMethod("calculate")
}

#' Plot Influence Analysis Results
#'
#' Generic plot function for 'influence_gam' objects. Dispatches to specific
#' plot types.
#'
#' @param x An object of class 'influence_gam'.
#' @param type The type of plot to generate. One of "stan", "step", "influ",
#'             "cdi", "step_influ".
#' @param ... Additional arguments passed to specific plot functions
#'            (e.g., `term` for "cdi").
#' @return A ggplot object or plots to the current device.
#' @export
plot.influence_gam <- function(x, type = "step_influ", ...) {
  switch(type,
         stan = plot_stan(x, ...),
         step = plot_step(x, ...),
         influ = plot_influ(x, ...),
         cdi = plot_cdi(x, ...),
         step_influ = plot_step_influ(x, ...),
         stop("Unknown plot type: ", type)
  )
}

#' Summarize Influence Analysis
#'
#' Provides a summary table of model fit statistics and influence metrics
#' for each term added sequentially.
#'
#' @param object An object of class 'influence_gam'.
#' @param ... Additional arguments (currently unused).
#' @return A data frame summarizing the analysis.
#' @export
# --- Calculation Method (Revised for ns()/bs() etc.) ---

#' @export
calculate.influence_gam <- function(x, islog = NULL, ...) {

  model <- x$model
  data <- x$data
  response_col <- x$response_data_col
  response_name <- x$response_name
  focus <- x$focus
  terms <- x$terms # These are the term labels from the model formula

  # Helper function to safely get prediction column names
  get_pred_colnames <- function(pred_obj) {
    if (!is.null(pred_obj) && !is.null(pred_obj$fit)) {
      colnames(pred_obj$fit)
    } else {
      character(0) # Return empty character vector if predictions failed
    }
  }

  # --- Get Observed Data ---
  observed <- data[[response_col]]
  if (inherits(observed, "Surv")) {
    observed <- as.numeric(observed[, 1])
  }

  # --- Check if Response is Logged ---
  if (is.null(islog)) {
    logged <- grepl("^log\\(", response_name) || grepl("^\\s*log\\(", response_name)
  } else {
    logged <- islog
  }
  if (!logged && any(observed <= 0, na.rm = TRUE)) {
     warning("Response variable contains non-positive values, but 'islog' is FALSE or auto-detected as FALSE. Geometric means cannot be calculated for the unstandardized index. Using arithmetic mean instead.")
     use_arithmetic_mean <- TRUE
  } else {
     use_arithmetic_mean <- FALSE
  }

  # --- Prepare Focus Levels ---
  focus_data <- data[[focus]]
  if (!is.factor(focus_data)) {
     if(is.character(focus_data)) {
        focus_data <- factor(focus_data)
        data[[focus]] <- focus_data # Update data frame if character converted
     } else {
        warning("Numeric focus term '", focus, "' detected. Treating levels as distinct categories for aggregation. Consider if binning is appropriate for interpretation.")
     }
     # Use unique sorted values for numeric, or levels for factor
     focus_levels <- sort(unique(focus_data))
     data$.focus_factor <- factor(focus_data, levels = focus_levels)
     focus_col_agg <- ".focus_factor"
  } else {
     focus_levels <- levels(focus_data)
     focus_col_agg <- focus
  }
  indices_df <- data.frame(level = focus_levels)
  if(is.factor(focus_data)) indices_df$level <- factor(indices_df$level, levels = levels(focus_data))


  # --- Calculate Unstandardized Index ---
  # (Logic unchanged from previous version)
  if (logged) {
    log_observed <- observed
  } else if (!use_arithmetic_mean) {
    log_observed <- log(observed)
  }

  if (use_arithmetic_mean) {
      unstan_agg <- aggregate(list(unstan = observed),
                             list(level = data[[focus_col_agg]]),
                             mean, na.rm = TRUE)
      mean_unstan <- mean(unstan_agg$unstan, na.rm=TRUE)
      indices_df <- merge(indices_df, unstan_agg, by = "level", sort = FALSE)
      indices_df$unstan <- indices_df$unstan / mean_unstan
  } else {
      unstan_agg <- aggregate(list(unstan = log_observed),
                              list(level = data[[focus_col_agg]]),
                              mean, na.rm = TRUE)
      mean_log_unstan <- mean(unstan_agg$unstan, na.rm=TRUE)
      indices_df <- merge(indices_df, unstan_agg, by = "level", sort = FALSE)
      indices_df$unstan <- exp(indices_df$unstan - mean_log_unstan)
  }


  # --- Calculate Standardized Index (from full model) ---
  pred_terms_full <- tryCatch({
      predict(model, type = "terms", se.fit = TRUE)
  }, error = function(e) {
      warning("Could not get term predictions for full model. Error: ", e$message)
      NULL # Return NULL on error
  })

  pred_colnames_full <- get_pred_colnames(pred_terms_full)

  # Identify the focus term column(s) in the predictions
  # Use direct matching first, as predict.gam usually names cols by term label
  focus_term_col <- NULL
  if (focus %in% pred_colnames_full) {
      focus_term_col <- focus
  } else {
      # Fallback: Check if focus is part of a term label (e.g., s(x, by=focus))
      # This is less likely needed if 'focus' itself is the primary term label of interest
      potential_matches <- pred_colnames_full[grepl(focus, pred_colnames_full, fixed = TRUE)]
      if (length(potential_matches) == 1) {
          focus_term_col <- potential_matches
          warning("Focus term '", focus, "' matched prediction column '", focus_term_col, "'. Verify this is correct.")
      } else if (length(potential_matches) > 1) {
           warning("Focus term '", focus, "' matched multiple prediction columns: ",
                   paste(potential_matches, collapse=", "), ". Cannot reliably determine standardized index.")
           focus_term_col <- NULL # Ambiguous
      } else {
          # Check if 'focus' is a *variable* used within a term (e.g., ns(focus, df=3))
          # The term label itself (e.g., "ns(focus, df=3)") should be in x$terms
          # We need to find which *term label* contains the focus variable
           term_containing_focus <- NULL
           for(term_lab in terms) {
               if(focus %in% all.vars(stats::reformulate(term_lab)[[2]])) {
                   if(term_lab %in% pred_colnames_full) {
                      term_containing_focus <- term_lab
                      break # Found the term label in predictions
                   }
               }
           }
           if (!is.null(term_containing_focus)) {
               focus_term_col <- term_containing_focus
               warning("Focus variable '", focus, "' found within term '", focus_term_col,
                       "'. Using this term for standardized index calculation.")
           } else {
                warning("Could not find a prediction column directly matching or clearly containing the focus term '", focus, "'. Standardized index calculation may fail.")
           }
      }
  }

  if (!is.null(focus_term_col) && !is.null(pred_terms_full)) {
      focus_effects_link <- pred_terms_full$fit[, focus_term_col]

      # Get Standard Error if available
      focus_effects_se <- NA
      if (!is.null(pred_terms_full$se.fit) && focus_term_col %in% colnames(pred_terms_full$se.fit)) {
           focus_effects_se <- pred_terms_full$se.fit[, focus_term_col]
      } else {
          warning("SE not found for column '", focus_term_col, "' in full model predictions.")
      }

      # Aggregate effects and SEs by focus level
      stan_agg_data <- data.frame(
            level = data[[focus_col_agg]],
            effect = focus_effects_link,
            se = focus_effects_se
      )
      stan_agg <- aggregate(cbind(effect, se) ~ level, data = stan_agg_data, mean, na.rm = TRUE)

      # Center the effects (link scale) and exponentiate
      base_link <- mean(stan_agg$effect, na.rm = TRUE)
      indices_df <- merge(indices_df, stan_agg, by = "level", sort = FALSE)
      indices_df$stan <- exp(indices_df$effect - base_link)

      # Approximate CI on response scale
      indices_df$stanLower <- exp((indices_df$effect - base_link) - 2 * indices_df$se)
      indices_df$stanUpper <- exp((indices_df$effect - base_link) + 2 * indices_df$se)

  } else {
      warning("Could not reliably identify focus term '", focus, "' effects in full model predictions. Standardized index/CI may be missing or NA.")
      indices_df$stan <- NA
      indices_df$stanLower <- NA
      indices_df$stanUpper <- NA
  }


  # --- Iteratively Add Terms and Calculate Indices/Metrics ---
  summary_list <- list()
  current_formula_terms <- character(0) # Start with empty terms for intercept model
  iter_terms_info <- c(list(list(label="intercept", actual="(Intercept)")), # Info for intercept step
                       lapply(terms, function(t) list(label=paste("+", t), actual=t))) # Info for subsequent terms

  # Initial null model (intercept only)
  null_formula <- reformulate("1", response = response_name)
  null_model <- tryCatch(
       gam(null_formula, data = data, family = family(model)),
       error = function(e) {warning("Failed to fit null model: ", e$message); NULL}
     )
   if(is.null(null_model)) {
       stop("Cannot proceed without fitting the null model.")
   }

  logLikeInterceptOnly <- if (!is.null(logLik(null_model))) logLik(null_model) else NA
  null_deviance <- if(!is.null(null_model$deviance)) null_model$deviance else NA
  n_obs <- nobs(model)
  target_y <- if (logged) log(observed) else observed
  target_y_finite <- is.finite(target_y)


  for (i in seq_along(iter_terms_info)) {
    step_info <- iter_terms_info[[i]]
    term_label_for_step <- step_info$label   # Label for this step (e.g., "+ year")
    term_actual_added <- step_info$actual  # Actual term ("year") or "(Intercept)"

    if (i == 1) { # Intercept model
        step_model <- null_model
        current_formula_terms <- character(0)
    } else {
        # Add the next term to the formula
        current_formula_terms <- c(current_formula_terms, term_actual_added)
        step_formula <- reformulate(current_formula_terms, response = response_name)

        # Fit model for this step
        step_model <- tryCatch(
            gam(step_formula, data = data, family = family(model)),
             error = function(e) {warning("Failed to fit model at step ", i, " (term: ", term_actual_added,"): ", e$message); NULL}
        )
    }

    # --- Calculate Metrics for step_model ---
    k_params <- NA; logLike <- NA; aic_val <- NA; deviance <- NA
    r2 <- NA; r2Dev <- NA; r2Negel <- NA
    index_col_name <- term_label_for_step # Use "+ term" or "intercept" as column name

    if (!is.null(step_model)) {
        k_params <- length(coef(step_model))
        logLike <- if (!is.null(logLik(step_model))) logLik(step_model) else NA
        aic_val <- if (!is.null(logLik(step_model))) AIC(step_model) else NA
        deviance <- if (!is.null(step_model$deviance)) step_model$deviance else NA

        # R-squared variants (logic unchanged)
        fitted_vals_link <- predict(step_model, newdata=data)
        if (logged) {
             if (length(na.omit(target_y)) > 2 && length(na.omit(fitted_vals_link)) > 2) {
                 finite_idx <- target_y_finite & is.finite(fitted_vals_link)
                 if(sum(finite_idx) > 2) r2 <- cor(target_y[finite_idx], fitted_vals_link[finite_idx])^2
             }
         } else {
              fitted_vals_response <- family(step_model)$linkinv(fitted_vals_link)
              if (length(na.omit(target_y)) > 2 && length(na.omit(fitted_vals_response)) > 2) {
                  finite_idx <- target_y_finite & is.finite(fitted_vals_response)
                  if(sum(finite_idx) > 2) r2 <- cor(target_y[finite_idx], fitted_vals_response[finite_idx])^2
              }
         }
        if (!is.na(deviance) && !is.na(null_deviance) && null_deviance > 0) r2Dev <- (null_deviance - deviance) / null_deviance
        if (!is.na(logLike) && !is.na(logLikeInterceptOnly) && !is.na(n_obs) && n_obs > 0) {
          term1_exp <- suppressWarnings(exp((logLikeInterceptOnly - logLike) * (2 / n_obs)))
          term2_exp <- suppressWarnings(exp(logLikeInterceptOnly * (2 / n_obs)))
          if(is.finite(term1_exp) && is.finite(term2_exp) && (1 - term2_exp) != 0) {
             r2Negel <- (1 - term1_exp) / (1 - term2_exp)
          } else { r2Negel <- NA }
        }

         # --- Calculate effects/index for this model step (relative to focus term) ---
         step_pred_terms <- tryCatch(
             predict(step_model, type = "terms"),
             error = function(e){ warning("Cannot get term preds for step ", i, ". Error: ", e$message); NULL}
         )
         step_pred_colnames <- if(!is.null(step_pred_terms)) colnames(step_pred_terms) else character(0)

         # Find focus term column in this step_model's predictions
         step_focus_term_col <- NULL
         if (focus %in% step_pred_colnames) {
             step_focus_term_col <- focus
         } else {
              # Check if focus var is within a term label present in this step model
              term_containing_focus_step <- NULL
              step_term_labels <- attr(terms(step_model), "term.labels") # Terms in *this* model
               for(term_lab in step_term_labels) {
                   if(focus %in% all.vars(stats::reformulate(term_lab)[[2]])) {
                       if(term_lab %in% step_pred_colnames) {
                          term_containing_focus_step <- term_lab
                          break
                       }
                   }
               }
              if (!is.null(term_containing_focus_step)) {
                 step_focus_term_col <- term_containing_focus_step
              }
         }


         if (!is.null(step_focus_term_col) && !is.null(step_pred_terms)) {
             step_focus_effects_link <- step_pred_terms[, step_focus_term_col]
             step_stan_agg_data <- data.frame(level = data[[focus_col_agg]], effect = step_focus_effects_link)
             step_stan_agg <- aggregate(effect ~ level, data = step_stan_agg_data, mean, na.rm = TRUE)

             step_base_link <- mean(step_stan_agg$effect, na.rm = TRUE)
             step_stan_agg[[index_col_name]] <- exp(step_stan_agg$effect - step_base_link)

             indices_df <- merge(indices_df, step_stan_agg[, c("level", index_col_name)], by = "level", all.x = TRUE, sort = FALSE)

         } else {
             # Focus term not in model yet or not found in predictions, index is flat (1)
             # Or prediction failed
              indices_df[[index_col_name]] <- if(is.null(step_pred_terms)) NA else 1.0
              if(!is.null(step_pred_terms)) { # Only warn if preds exist but term wasn't found
                 # This is expected if focus term hasn't been added yet
                 # warning("Focus term '", focus, "' not found in step ", i, " predictions. Index set to 1.")
              } else {
                  warning("Index calculation failed for step: ", index_col_name, " due to prediction error.")
              }
         }

    } else { # step_model fitting failed
        indices_df[[index_col_name]] <- NA # Mark index as NA
        warning("Model fitting failed for step ", i, ". Metrics and index set to NA.")
    }

     # Store summary row
     summary_list[[i]] <- data.frame(
          term_added = term_label_for_step,
          term_actual = term_actual_added,
          k = k_params,
          logLik = logLike,
          AIC = aic_val,
          R2_cor = r2,
          R2_dev = r2Dev,
          R2_nagel = r2Negel
        )

  } # End loop through terms

  summary_df <- bind_rows(summary_list)

  # Calculate differences for summary table
  summary_df <- summary_df %>%
      mutate(
          k_diff = c(NA, diff(k)),
          R2_cor_diff = c(NA, diff(R2_cor)),
          R2_dev_diff = c(NA, diff(R2_dev)),
          R2_nagel_diff = c(NA, diff(R2_nagel))
      ) %>%
      select(term_added, term_actual, k, k_diff, logLik, AIC,
             R2_cor, R2_cor_diff, R2_dev, R2_dev_diff, R2_nagel, R2_nagel_diff)


  # --- Calculate Term Predictions and Influences (from full model) ---
  influences_df <- data.frame(level = focus_levels)
  if(is.factor(focus_data)) influences_df$level <- factor(influences_df$level, levels = levels(focus_data))

  overall_influ <- numeric(length(terms))
  trend_influ <- numeric(length(terms))
  names(overall_influ) <- terms
  names(trend_influ) <- terms

  # Reuse full model predictions from earlier
  if (!is.null(pred_terms_full)) {
      pred_df <- as.data.frame(pred_terms_full$fit)
      pred_se_df <- if(!is.null(pred_terms_full$se.fit)) as.data.frame(pred_terms_full$se.fit) else NULL
      pred_colnames <- colnames(pred_df) # Use colnames from the prediction matrix

      # Combine data with predictions, handling potential NAs
      full_pred_data <- NULL
      if(nrow(pred_df) == nrow(data)) {
          full_pred_data <- bind_cols(data, pred_df)
          if(!is.null(pred_se_df)) {
              names(pred_se_df) <- paste0(pred_colnames, "_se") # Ensure SE names match fit names
              full_pred_data <- bind_cols(full_pred_data, pred_se_df)
          }
      } else if (!is.null(model$na.action)) {
           # (NA handling logic unchanged)
           warning("NA actions detected. Aligning predictions with original data. Check results carefully.")
           na_rows <- as.integer(model$na.action)
           padded_pred_df <- matrix(NA, nrow = nrow(data), ncol = ncol(pred_df)); colnames(padded_pred_df) <- pred_colnames
           padded_pred_df[-na_rows, ] <- as.matrix(pred_df)
           full_pred_data <- bind_cols(data, as.data.frame(padded_pred_df))
           if(!is.null(pred_se_df)) {
               padded_pred_se_df <- matrix(NA, nrow = nrow(data), ncol = ncol(pred_se_df)); colnames(padded_pred_se_df) <- paste0(pred_colnames, "_se")
               padded_pred_se_df[-na_rows, ] <- as.matrix(pred_se_df)
               full_pred_data <- bind_cols(full_pred_data, as.data.frame(padded_pred_se_df))
           }
       } else {
          warning("Prediction rows (", nrow(pred_df), ") do not match data rows (", nrow(data), ") and no NA action found. Influence calculations may be incorrect.")
          full_pred_data <- data # Placeholder
       }


      x$predictions <- full_pred_data # Store predictions with data

      # Calculate influences per term (excluding focus term itself)
      for (term in terms) { # Iterate through term labels from the original model
          if (term == focus) { # Skip influence calculation for the focus term itself
              influences_df[[term]] <- NA
              overall_influ[term] <- NA
              trend_influ[term] <- NA
              next # Move to the next term
          }

          # Check if this term exists as a column in the prediction matrix
          if (term %in% pred_colnames && !is.null(full_pred_data)) {
              term_effect_link <- full_pred_data[[term]] # Directly use the column

              # Aggregate influence by focus level
              infl_agg_data <- data.frame(level = full_pred_data[[focus_col_agg]], value = term_effect_link)
              infl_agg <- aggregate(value ~ level, data = infl_agg_data, mean, na.rm = TRUE)
              names(infl_agg)[names(infl_agg)=="value"] <- term # Rename col to term name

              # Calculate overall and trend influence stats (logic unchanged)
              if (nrow(infl_agg) > 1) {
                 overall_influ[term] <- exp(mean(abs(infl_agg[[term]]), na.rm = TRUE)) - 1
                 level_num <- if (is.factor(infl_agg$level)) as.numeric(infl_agg$level) else infl_agg$level
                 valid_trend_data <- is.finite(level_num) & is.finite(infl_agg[[term]])
                 if(sum(valid_trend_data) > 1) {
                    trend_cov <- cov(level_num[valid_trend_data], infl_agg[[term]][valid_trend_data])
                    trend_var <- var(level_num[valid_trend_data])
                    if (is.finite(trend_cov) && is.finite(trend_var) && trend_var != 0) {
                       trend_influ[term] <- exp(trend_cov / trend_var) - 1
                    } else { trend_influ[term] <- NA }
                 } else { trend_influ[term] <- NA }
              } else { overall_influ[term] <- NA; trend_influ[term] <- NA }

              influences_df <- merge(influences_df, infl_agg[, c("level", term)], by = "level", all.x = TRUE, sort=FALSE)

          } else { # Term not found in prediction columns or prediction failed
              warning("Influence calculation skipped for term '", term, "' as it was not found in prediction columns or predictions failed.")
               influences_df[[term]] <- NA # Add NA column
               overall_influ[term] <- NA
               trend_influ[term] <- NA
          }
      } # End loop terms for influence

  } else { # Full model prediction failed earlier
      warning("Cannot calculate influences due to full model prediction failure.")
      x$predictions <- data # Store original data only
      # Fill influence stats with NA
      overall_influ[] <- NA
      trend_influ[] <- NA
      # Add NA columns to influences_df
      for(term in terms) if(term != focus) influences_df[[term]] <- NA
  }

  # Add influence stats to summary_df
  influ_stats <- data.frame(term_actual = terms, overall = overall_influ, trend = trend_influ)
  influ_stats <- bind_rows(data.frame(term_actual="(Intercept)", overall=NA, trend=NA), influ_stats)
  summary_df$term_actual <- as.character(summary_df$term_actual)
  influ_stats$term_actual <- as.character(influ_stats$term_actual)
  summary_df <- left_join(summary_df, influ_stats, by = "term_actual")

  # --- Store Results ---
  x$indices <- indices_df
  x$summary <- summary_df
  x$influences <- influences_df

  invisible(x) # Return updated object invisibly
}
  if (is.null(object$summary)) {
    stop("Calculations not performed yet. Run calculate() first.")
  }
  print(object$summary)
  invisible(object$summary)
}
# --- Calculation Method (Revised for ns()/bs() etc.) ---

#' @export
calculate.influence_gam <- function(x, islog = NULL, ...) {

  model <- x$model
  data <- x$data
  response_col <- x$response_data_col
  response_name <- x$response_name
  focus <- x$focus
  terms <- x$terms # These are the term labels from the model formula

  # Helper function to safely get prediction column names
  get_pred_colnames <- function(pred_obj) {
    if (!is.null(pred_obj) && !is.null(pred_obj$fit)) {
      colnames(pred_obj$fit)
    } else {
      character(0) # Return empty character vector if predictions failed
    }
  }

  # --- Get Observed Data ---
  observed <- data[[response_col]]
  if (inherits(observed, "Surv")) {
    observed <- as.numeric(observed[, 1])
  }

  # --- Check if Response is Logged ---
  if (is.null(islog)) {
    logged <- grepl("^log\\(", response_name) || grepl("^\\s*log\\(", response_name)
  } else {
    logged <- islog
  }
  if (!logged && any(observed <= 0, na.rm = TRUE)) {
     warning("Response variable contains non-positive values, but 'islog' is FALSE or auto-detected as FALSE. Geometric means cannot be calculated for the unstandardized index. Using arithmetic mean instead.")
     use_arithmetic_mean <- TRUE
  } else {
     use_arithmetic_mean <- FALSE
  }

  # --- Prepare Focus Levels ---
  focus_data <- data[[focus]]
  if (!is.factor(focus_data)) {
     if(is.character(focus_data)) {
        focus_data <- factor(focus_data)
        data[[focus]] <- focus_data # Update data frame if character converted
     } else {
        warning("Numeric focus term '", focus, "' detected. Treating levels as distinct categories for aggregation. Consider if binning is appropriate for interpretation.")
     }
     # Use unique sorted values for numeric, or levels for factor
     focus_levels <- sort(unique(focus_data))
     data$.focus_factor <- factor(focus_data, levels = focus_levels)
     focus_col_agg <- ".focus_factor"
  } else {
     focus_levels <- levels(focus_data)
     focus_col_agg <- focus
  }
  indices_df <- data.frame(level = focus_levels)
  if(is.factor(focus_data)) indices_df$level <- factor(indices_df$level, levels = levels(focus_data))


  # --- Calculate Unstandardized Index ---
  # (Logic unchanged from previous version)
  if (logged) {
    log_observed <- observed
  } else if (!use_arithmetic_mean) {
    log_observed <- log(observed)
  }

  if (use_arithmetic_mean) {
      unstan_agg <- aggregate(list(unstan = observed),
                             list(level = data[[focus_col_agg]]),
                             mean, na.rm = TRUE)
      mean_unstan <- mean(unstan_agg$unstan, na.rm=TRUE)
      indices_df <- merge(indices_df, unstan_agg, by = "level", sort = FALSE)
      indices_df$unstan <- indices_df$unstan / mean_unstan
  } else {
      unstan_agg <- aggregate(list(unstan = log_observed),
                              list(level = data[[focus_col_agg]]),
                              mean, na.rm = TRUE)
      mean_log_unstan <- mean(unstan_agg$unstan, na.rm=TRUE)
      indices_df <- merge(indices_df, unstan_agg, by = "level", sort = FALSE)
      indices_df$unstan <- exp(indices_df$unstan - mean_log_unstan)
  }


  # --- Calculate Standardized Index (from full model) ---
  pred_terms_full <- tryCatch({
      predict(model, type = "terms", se.fit = TRUE)
  }, error = function(e) {
      warning("Could not get term predictions for full model. Error: ", e$message)
      NULL # Return NULL on error
  })

  pred_colnames_full <- get_pred_colnames(pred_terms_full)

  # Identify the focus term column(s) in the predictions
  # Use direct matching first, as predict.gam usually names cols by term label
  focus_term_col <- NULL
  if (focus %in% pred_colnames_full) {
      focus_term_col <- focus
  } else {
      # Fallback: Check if focus is part of a term label (e.g., s(x, by=focus))
      # This is less likely needed if 'focus' itself is the primary term label of interest
      potential_matches <- pred_colnames_full[grepl(focus, pred_colnames_full, fixed = TRUE)]
      if (length(potential_matches) == 1) {
          focus_term_col <- potential_matches
          warning("Focus term '", focus, "' matched prediction column '", focus_term_col, "'. Verify this is correct.")
      } else if (length(potential_matches) > 1) {
           warning("Focus term '", focus, "' matched multiple prediction columns: ",
                   paste(potential_matches, collapse=", "), ". Cannot reliably determine standardized index.")
           focus_term_col <- NULL # Ambiguous
      } else {
          # Check if 'focus' is a *variable* used within a term (e.g., ns(focus, df=3))
          # The term label itself (e.g., "ns(focus, df=3)") should be in x$terms
          # We need to find which *term label* contains the focus variable
           term_containing_focus <- NULL
           for(term_lab in terms) {
               if(focus %in% all.vars(stats::reformulate(term_lab)[[2]])) {
                   if(term_lab %in% pred_colnames_full) {
                      term_containing_focus <- term_lab
                      break # Found the term label in predictions
                   }
               }
           }
           if (!is.null(term_containing_focus)) {
               focus_term_col <- term_containing_focus
               warning("Focus variable '", focus, "' found within term '", focus_term_col,
                       "'. Using this term for standardized index calculation.")
           } else {
                warning("Could not find a prediction column directly matching or clearly containing the focus term '", focus, "'. Standardized index calculation may fail.")
           }
      }
  }

  if (!is.null(focus_term_col) && !is.null(pred_terms_full)) {
      focus_effects_link <- pred_terms_full$fit[, focus_term_col]

      # Get Standard Error if available
      focus_effects_se <- NA
      if (!is.null(pred_terms_full$se.fit) && focus_term_col %in% colnames(pred_terms_full$se.fit)) {
           focus_effects_se <- pred_terms_full$se.fit[, focus_term_col]
      } else {
          warning("SE not found for column '", focus_term_col, "' in full model predictions.")
      }

      # Aggregate effects and SEs by focus level
      stan_agg_data <- data.frame(
            level = data[[focus_col_agg]],
            effect = focus_effects_link,
            se = focus_effects_se
      )
      stan_agg <- aggregate(cbind(effect, se) ~ level, data = stan_agg_data, mean, na.rm = TRUE)

      # Center the effects (link scale) and exponentiate
      base_link <- mean(stan_agg$effect, na.rm = TRUE)
      indices_df <- merge(indices_df, stan_agg, by = "level", sort = FALSE)
      indices_df$stan <- exp(indices_df$effect - base_link)

      # Approximate CI on response scale
      indices_df$stanLower <- exp((indices_df$effect - base_link) - 2 * indices_df$se)
      indices_df$stanUpper <- exp((indices_df$effect - base_link) + 2 * indices_df$se)

  } else {
      warning("Could not reliably identify focus term '", focus, "' effects in full model predictions. Standardized index/CI may be missing or NA.")
      indices_df$stan <- NA
      indices_df$stanLower <- NA
      indices_df$stanUpper <- NA
  }


  # --- Iteratively Add Terms and Calculate Indices/Metrics ---
  summary_list <- list()
  current_formula_terms <- character(0) # Start with empty terms for intercept model
  iter_terms_info <- c(list(list(label="intercept", actual="(Intercept)")), # Info for intercept step
                       lapply(terms, function(t) list(label=paste("+", t), actual=t))) # Info for subsequent terms

  # Initial null model (intercept only)
  null_formula <- reformulate("1", response = response_name)
  null_model <- tryCatch(
       gam(null_formula, data = data, family = family(model)),
       error = function(e) {warning("Failed to fit null model: ", e$message); NULL}
     )
   if(is.null(null_model)) {
       stop("Cannot proceed without fitting the null model.")
   }

  logLikeInterceptOnly <- if (!is.null(logLik(null_model))) logLik(null_model) else NA
  null_deviance <- if(!is.null(null_model$deviance)) null_model$deviance else NA
  n_obs <- nobs(model)
  target_y <- if (logged) log(observed) else observed
  target_y_finite <- is.finite(target_y)


  for (i in seq_along(iter_terms_info)) {
    step_info <- iter_terms_info[[i]]
    term_label_for_step <- step_info$label   # Label for this step (e.g., "+ year")
    term_actual_added <- step_info$actual  # Actual term ("year") or "(Intercept)"

    if (i == 1) { # Intercept model
        step_model <- null_model
        current_formula_terms <- character(0)
    } else {
        # Add the next term to the formula
        current_formula_terms <- c(current_formula_terms, term_actual_added)
        step_formula <- reformulate(current_formula_terms, response = response_name)

        # Fit model for this step
        step_model <- tryCatch(
            gam(step_formula, data = data, family = family(model)),
             error = function(e) {warning("Failed to fit model at step ", i, " (term: ", term_actual_added,"): ", e$message); NULL}
        )
    }

    # --- Calculate Metrics for step_model ---
    k_params <- NA; logLike <- NA; aic_val <- NA; deviance <- NA
    r2 <- NA; r2Dev <- NA; r2Negel <- NA
    index_col_name <- term_label_for_step # Use "+ term" or "intercept" as column name

    if (!is.null(step_model)) {
        k_params <- length(coef(step_model))
        logLike <- if (!is.null(logLik(step_model))) logLik(step_model) else NA
        aic_val <- if (!is.null(logLik(step_model))) AIC(step_model) else NA
        deviance <- if (!is.null(step_model$deviance)) step_model$deviance else NA

        # R-squared variants (logic unchanged)
        fitted_vals_link <- predict(step_model, newdata=data)
        if (logged) {
             if (length(na.omit(target_y)) > 2 && length(na.omit(fitted_vals_link)) > 2) {
                 finite_idx <- target_y_finite & is.finite(fitted_vals_link)
                 if(sum(finite_idx) > 2) r2 <- cor(target_y[finite_idx], fitted_vals_link[finite_idx])^2
             }
         } else {
              fitted_vals_response <- family(step_model)$linkinv(fitted_vals_link)
              if (length(na.omit(target_y)) > 2 && length(na.omit(fitted_vals_response)) > 2) {
                  finite_idx <- target_y_finite & is.finite(fitted_vals_response)
                  if(sum(finite_idx) > 2) r2 <- cor(target_y[finite_idx], fitted_vals_response[finite_idx])^2
              }
         }
        if (!is.na(deviance) && !is.na(null_deviance) && null_deviance > 0) r2Dev <- (null_deviance - deviance) / null_deviance
        if (!is.na(logLike) && !is.na(logLikeInterceptOnly) && !is.na(n_obs) && n_obs > 0) {
          term1_exp <- suppressWarnings(exp((logLikeInterceptOnly - logLike) * (2 / n_obs)))
          term2_exp <- suppressWarnings(exp(logLikeInterceptOnly * (2 / n_obs)))
          if(is.finite(term1_exp) && is.finite(term2_exp) && (1 - term2_exp) != 0) {
             r2Negel <- (1 - term1_exp) / (1 - term2_exp)
          } else { r2Negel <- NA }
        }

         # --- Calculate effects/index for this model step (relative to focus term) ---
         step_pred_terms <- tryCatch(
             predict(step_model, type = "terms"),
             error = function(e){ warning("Cannot get term preds for step ", i, ". Error: ", e$message); NULL}
         )
         step_pred_colnames <- if(!is.null(step_pred_terms)) colnames(step_pred_terms) else character(0)

         # Find focus term column in this step_model's predictions
         step_focus_term_col <- NULL
         if (focus %in% step_pred_colnames) {
             step_focus_term_col <- focus
         } else {
              # Check if focus var is within a term label present in this step model
              term_containing_focus_step <- NULL
              step_term_labels <- attr(terms(step_model), "term.labels") # Terms in *this* model
               for(term_lab in step_term_labels) {
                   if(focus %in% all.vars(stats::reformulate(term_lab)[[2]])) {
                       if(term_lab %in% step_pred_colnames) {
                          term_containing_focus_step <- term_lab
                          break
                       }
                   }
               }
              if (!is.null(term_containing_focus_step)) {
                 step_focus_term_col <- term_containing_focus_step
              }
         }


         if (!is.null(step_focus_term_col) && !is.null(step_pred_terms)) {
             step_focus_effects_link <- step_pred_terms[, step_focus_term_col]
             step_stan_agg_data <- data.frame(level = data[[focus_col_agg]], effect = step_focus_effects_link)
             step_stan_agg <- aggregate(effect ~ level, data = step_stan_agg_data, mean, na.rm = TRUE)

             step_base_link <- mean(step_stan_agg$effect, na.rm = TRUE)
             step_stan_agg[[index_col_name]] <- exp(step_stan_agg$effect - step_base_link)

             indices_df <- merge(indices_df, step_stan_agg[, c("level", index_col_name)], by = "level", all.x = TRUE, sort = FALSE)

         } else {
             # Focus term not in model yet or not found in predictions, index is flat (1)
             # Or prediction failed
              indices_df[[index_col_name]] <- if(is.null(step_pred_terms)) NA else 1.0
              if(!is.null(step_pred_terms)) { # Only warn if preds exist but term wasn't found
                 # This is expected if focus term hasn't been added yet
                 # warning("Focus term '", focus, "' not found in step ", i, " predictions. Index set to 1.")
              } else {
                  warning("Index calculation failed for step: ", index_col_name, " due to prediction error.")
              }
         }

    } else { # step_model fitting failed
        indices_df[[index_col_name]] <- NA # Mark index as NA
        warning("Model fitting failed for step ", i, ". Metrics and index set to NA.")
    }

     # Store summary row
     summary_list[[i]] <- data.frame(
          term_added = term_label_for_step,
          term_actual = term_actual_added,
          k = k_params,
          logLik = logLike,
          AIC = aic_val,
          R2_cor = r2,
          R2_dev = r2Dev,
          R2_nagel = r2Negel
        )

  } # End loop through terms

  summary_df <- bind_rows(summary_list)

  # Calculate differences for summary table
  summary_df <- summary_df %>%
      mutate(
          k_diff = c(NA, diff(k)),
          R2_cor_diff = c(NA, diff(R2_cor)),
          R2_dev_diff = c(NA, diff(R2_dev)),
          R2_nagel_diff = c(NA, diff(R2_nagel))
      ) %>%
      select(term_added, term_actual, k, k_diff, logLik, AIC,
             R2_cor, R2_cor_diff, R2_dev, R2_dev_diff, R2_nagel, R2_nagel_diff)


  # --- Calculate Term Predictions and Influences (from full model) ---
  influences_df <- data.frame(level = focus_levels)
  if(is.factor(focus_data)) influences_df$level <- factor(influences_df$level, levels = levels(focus_data))

  overall_influ <- numeric(length(terms))
  trend_influ <- numeric(length(terms))
  names(overall_influ) <- terms
  names(trend_influ) <- terms

  # Reuse full model predictions from earlier
  if (!is.null(pred_terms_full)) {
      pred_df <- as.data.frame(pred_terms_full$fit)
      pred_se_df <- if(!is.null(pred_terms_full$se.fit)) as.data.frame(pred_terms_full$se.fit) else NULL
      pred_colnames <- colnames(pred_df) # Use colnames from the prediction matrix

      # Combine data with predictions, handling potential NAs
      full_pred_data <- NULL
      if(nrow(pred_df) == nrow(data)) {
          full_pred_data <- bind_cols(data, pred_df)
          if(!is.null(pred_se_df)) {
              names(pred_se_df) <- paste0(pred_colnames, "_se") # Ensure SE names match fit names
              full_pred_data <- bind_cols(full_pred_data, pred_se_df)
          }
      } else if (!is.null(model$na.action)) {
           # (NA handling logic unchanged)
           warning("NA actions detected. Aligning predictions with original data. Check results carefully.")
           na_rows <- as.integer(model$na.action)
           padded_pred_df <- matrix(NA, nrow = nrow(data), ncol = ncol(pred_df)); colnames(padded_pred_df) <- pred_colnames
           padded_pred_df[-na_rows, ] <- as.matrix(pred_df)
           full_pred_data <- bind_cols(data, as.data.frame(padded_pred_df))
           if(!is.null(pred_se_df)) {
               padded_pred_se_df <- matrix(NA, nrow = nrow(data), ncol = ncol(pred_se_df)); colnames(padded_pred_se_df) <- paste0(pred_colnames, "_se")
               padded_pred_se_df[-na_rows, ] <- as.matrix(pred_se_df)
               full_pred_data <- bind_cols(full_pred_data, as.data.frame(padded_pred_se_df))
           }
       } else {
          warning("Prediction rows (", nrow(pred_df), ") do not match data rows (", nrow(data), ") and no NA action found. Influence calculations may be incorrect.")
          full_pred_data <- data # Placeholder
       }


      x$predictions <- full_pred_data # Store predictions with data

      # Calculate influences per term (excluding focus term itself)
      for (term in terms) { # Iterate through term labels from the original model
          if (term == focus) { # Skip influence calculation for the focus term itself
              influences_df[[term]] <- NA
              overall_influ[term] <- NA
              trend_influ[term] <- NA
              next # Move to the next term
          }

          # Check if this term exists as a column in the prediction matrix
          if (term %in% pred_colnames && !is.null(full_pred_data)) {
              term_effect_link <- full_pred_data[[term]] # Directly use the column

              # Aggregate influence by focus level
              infl_agg_data <- data.frame(level = full_pred_data[[focus_col_agg]], value = term_effect_link)
              infl_agg <- aggregate(value ~ level, data = infl_agg_data, mean, na.rm = TRUE)
              names(infl_agg)[names(infl_agg)=="value"] <- term # Rename col to term name

              # Calculate overall and trend influence stats (logic unchanged)
              if (nrow(infl_agg) > 1) {
                 overall_influ[term] <- exp(mean(abs(infl_agg[[term]]), na.rm = TRUE)) - 1
                 level_num <- if (is.factor(infl_agg$level)) as.numeric(infl_agg$level) else infl_agg$level
                 valid_trend_data <- is.finite(level_num) & is.finite(infl_agg[[term]])
                 if(sum(valid_trend_data) > 1) {
                    trend_cov <- cov(level_num[valid_trend_data], infl_agg[[term]][valid_trend_data])
                    trend_var <- var(level_num[valid_trend_data])
                    if (is.finite(trend_cov) && is.finite(trend_var) && trend_var != 0) {
                       trend_influ[term] <- exp(trend_cov / trend_var) - 1
                    } else { trend_influ[term] <- NA }
                 } else { trend_influ[term] <- NA }
              } else { overall_influ[term] <- NA; trend_influ[term] <- NA }

              influences_df <- merge(influences_df, infl_agg[, c("level", term)], by = "level", all.x = TRUE, sort=FALSE)

          } else { # Term not found in prediction columns or prediction failed
              warning("Influence calculation skipped for term '", term, "' as it was not found in prediction columns or predictions failed.")
               influences_df[[term]] <- NA # Add NA column
               overall_influ[term] <- NA
               trend_influ[term] <- NA
          }
      } # End loop terms for influence

  } else { # Full model prediction failed earlier
      warning("Cannot calculate influences due to full model prediction failure.")
      x$predictions <- data # Store original data only
      # Fill influence stats with NA
      overall_influ[] <- NA
      trend_influ[] <- NA
      # Add NA columns to influences_df
      for(term in terms) if(term != focus) influences_df[[term]] <- NA
  }

  # Add influence stats to summary_df
  influ_stats <- data.frame(term_actual = terms, overall = overall_influ, trend = trend_influ)
  influ_stats <- bind_rows(data.frame(term_actual="(Intercept)", overall=NA, trend=NA), influ_stats)
  summary_df$term_actual <- as.character(summary_df$term_actual)
  influ_stats$term_actual <- as.character(influ_stats$term_actual)
  summary_df <- left_join(summary_df, influ_stats, by = "term_actual")

  # --- Store Results ---
  x$indices <- indices_df
  x$summary <- summary_df
  x$influences <- influences_df

  invisible(x) # Return updated object invisibly
}

#' Print Influence GAM Object
#'
#' Prints basic information about the influence_gam object.
#'
#' @param x An object of class 'influence_gam'.
#' @param ... Additional arguments (currently unused).
#' @export
print.influence_gam <- function(x, ...) {
  cat("Influence Analysis Object (mgcv)\n")
  cat("---------------------------------\n")
  cat("Model Call:", deparse(x$model$call), "\n")
  cat("Focus Term:", x$focus, "\n")
  cat("Response:", x$response_name, "\n")
  cat("Terms:", paste(x$terms, collapse = ", "), "\n")
  if (!is.null(x$summary)) {
    cat("Status: Calculations complete. Use summary() or plot().\n")
  } else {
    cat("Status: Calculations pending. Run calculate().\n")
  }
  invisible(x)
}


# --- S3 Constructor ---

#' Create an Influence Object for GAMs
#'
#' Initializes an object to analyze the influence of terms in a GAM model
#' fitted with `mgcv::gam`.
#'
#' @param model A fitted GAM model object from `mgcv::gam`.
#' @param focus The name (character string) of the term in the model for which
#'              influence is the primary focus (e.g., "year"). This term's
#'              levels/values will form the x-axis of many plots. It must be
#'              present in `model`.
#' @param data The data frame used to fit the model. Often extractable from
#'             the model object, but required if `model$data` is NULL.
#' @param response Optional: The name of the response variable. If NULL, it's
#'                 inferred from the model formula.
#' @param term_labels Optional: A named list or vector to provide custom labels
#'                   for terms in plots (e.g., `list(fac_A = "Factor A")`).
#' @param term_orders Optional: A named list or vector specifying the order for
#'                   levels in CDI plots ("asis" or "coef"). Default is "asis".
#' @return An object of class 'influence_gam'. Run `calculate()` on this object
#'         before plotting or summarizing.
#' @export
#' @examples
#' \dontrun{
#' library(mgcv)
#' library(ggplot2)
#' set.seed(123)
#' n <- 200
#' dat <- gamSim(1, n = n, scale = 2)
#' dat$year <- factor(sample(2001:2010, n, replace = TRUE))
#' dat$fac <- factor(sample(letters[1:4], n, replace = TRUE))
#'
#' # Add a log-normal response
#' dat$log_catch <- rnorm(n, mean = (dat$y / 5) + as.numeric(dat$year)*0.05, sd = 0.5)
#' dat$catch <- exp(dat$log_catch)
#'
#' # Fit a GAM
#' m1 <- gam(log_catch ~ year + s(x0) + s(x1) + s(x2) + fac, data = dat)
#'
#' # Create influence object
#' influ_obj <- influence_gam(m1, focus = "year", data = dat)
#'
#' # Calculate metrics
#' influ_obj <- calculate(influ_obj)
#'
#' # Explore results
#' summary(influ_obj)
#' plot(influ_obj, type = "stan")
#' plot(influ_obj, type = "step")
#' plot(influ_obj, type = "influ")
#' plot(influ_obj, type = "cdi", term = "s(x0)")
#' plot(influ_obj, type = "cdi", term = "fac")
#' plot(influ_obj, type = "step_influ")
#' }
influence_gam <- function(model, focus, data = NULL, response = NULL,
                          term_labels = NULL, term_orders = NULL) {

  if (!inherits(model, "gam")) {
    stop("Model must be of class 'gam' (from mgcv package).")
  }
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Please install and load the 'mgcv' package.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install and load the 'ggplot2' package.")
  }

  # --- Data Extraction ---
  model_data <- model$model # This is the model frame
  if (is.null(model_data)) {
      if (is.null(data)) {
        stop("Model frame not found in gam object and 'data' argument not provided.")
      }
      # Need to reconstruct model frame if not stored and data is provided
      # This is simplified; might need terms handling for complex formulas
      # model_vars <- all.vars(formula(model))
      # model_data <- data[, model_vars, drop = FALSE] # Basic attempt
      model_data <- data # Assume provided data has all necessary columns
      # It's generally better to ensure gam stores the data or model frame
      warning("Using provided 'data'. Ensure it matches the data used for fitting exactly.")
  } else {
      data <- model$model # Use model frame if available
  }


  # --- Term and Response Identification ---
  model_terms <- attr(terms(model), "term.labels")
  if (!focus %in% model_terms && !focus %in% names(model$var.summary) ) { # Check if focus might be a smooth term variable
       # Check if focus is *part* of a term (e.g., a variable in s(var))
       var_in_term <- vapply(model_terms, function(tm) focus %in% all.vars(parse(text = tm)[[1]]), logical(1))
       if (!any(var_in_term)) {
            stop("Focus term '", focus, "' not found in model terms: ", paste(model_terms, collapse=", "))
       }
       # If focus is a variable within a smooth/interaction, use the variable name.
       # The code below relies on the *column name* in the data.
       warning("Focus '", focus, "' seems to be a variable within a term. Proceeding using the variable name.")
       # No change needed here as downstream code uses focus to index data[[focus]]
  }


  if (is.null(response)) {
    response_name <- as.character(formula(model)[[2]])
    # Handle response transformations like log(y) in formula
    response_data_col <- all.vars(formula(model)[[2]])[1]
  } else {
    response_name <- response
    response_data_col <- response
  }
   if(!response_data_col %in% names(data)) {
       stop("Response column '", response_data_col, "' not found in the provided data.")
   }
   # Check if focus column exists in data
    if(!focus %in% names(data)) {
        stop("Focus column '", focus, "' not found in the provided data.")
    }
   # Ensure focus term is a factor if expected to be categorical
   if (!is.factor(data[[focus]]) && !is.character(data[[focus]])) {
        warning("Focus term '", focus, "' is numeric. Plots will treat it as continuous or may require binning. Ensure this is intended.")
        # Consider converting character to factor
        if(is.character(data[[focus]])) data[[focus]] <- factor(data[[focus]])
   } else if (is.character(data[[focus]])) {
        data[[focus]] <- factor(data[[focus]])
   }


  # --- Initialize Labels and Orders ---
  labels <- as.list(model_terms)
  names(labels) <- model_terms
  if (!is.null(term_labels)) {
    labels[names(term_labels)] <- term_labels
  }
  # Add label for focus term variable if not already a term label
  if (!focus %in% names(labels)) labels[[focus]] <- focus


  orders <- as.list(rep("asis", length(model_terms)))
  names(orders) <- model_terms
  if (!is.null(term_orders)) {
      # Ensure provided orders are valid
      valid_orders <- names(term_orders) %in% model_terms
      if(!all(valid_orders)) {
          warning("Ignoring term_orders for terms not in the model: ",
                  paste(names(term_orders)[!valid_orders], collapse = ", "))
          term_orders <- term_orders[valid_orders]
      }
      invalid_values <- !sapply(term_orders, function(x) x %in% c("asis", "coef"))
      if (any(invalid_values)) {
           warning("Ignoring invalid term_orders values (must be 'asis' or 'coef'): ",
                   paste(term_orders[invalid_values], collapse = ", "))
           term_orders <- term_orders[!invalid_values]
      }
    orders[names(term_orders)] <- term_orders
  }
   # Add order for focus term variable if not already a term label
  if (!focus %in% names(orders)) orders[[focus]] <- "asis"


  # --- Create Structure ---
  structure(
    list(
      model = model,
      data = data,
      focus = focus,
      response_name = response_name, # Name potentially including log()
      response_data_col = response_data_col, # Name of column in data
      terms = model_terms,
      labels = labels,
      orders = orders,
      # Placeholders for results
      indices = NULL,
      summary = NULL,
      predictions = NULL,
      influences = NULL
    ),
    class = "influence_gam"
  )
}


# --- Calculation Method ---

#' @export
calculate.influence_gam <- function(x, islog = NULL, ...) {

  model <- x$model
  data <- x$data
  response_col <- x$response_data_col
  response_name <- x$response_name
  focus <- x$focus
  terms <- x$terms

  # --- Get Observed Data ---
  observed <- data[[response_col]]
  # Basic check for Surv object (though less common with GAMs for this purpose)
  if (inherits(observed, "Surv")) {
    observed <- as.numeric(observed[, 1])
  }

  # --- Check if Response is Logged ---
  if (is.null(islog)) {
    # Infer from response name in formula
    logged <- grepl("^log\\(", response_name) || grepl("^\\s*log\\(", response_name)
  } else {
    logged <- islog
  }

  if (!logged && any(observed <= 0, na.rm = TRUE)) {
     warning("Response variable contains non-positive values, but 'islog' is FALSE or auto-detected as FALSE. Geometric means cannot be calculated for the unstandardized index. Using arithmetic mean instead.")
     use_arithmetic_mean <- TRUE
  } else {
     use_arithmetic_mean <- FALSE
  }

  # --- Prepare Focus Levels ---
  focus_data <- data[[focus]]
  if (!is.factor(focus_data)) {
     # Convert numeric focus to factor for aggregation, keeping original numeric for potential later use
     # Or handle numeric focus differently depending on desired output
     warning("Numeric focus term '", focus, "' detected. Treating levels as distinct categories for aggregation. Consider if binning is appropriate for interpretation.")
     focus_levels <- sort(unique(focus_data))
     data$.focus_factor <- factor(focus_data, levels = focus_levels)
     focus_col_agg <- ".focus_factor" # Use the factor for aggregation
  } else {
     focus_levels <- levels(focus_data)
     focus_col_agg <- focus # Use original factor name
  }

  # Create initial data frame for indices
  indices_df <- data.frame(level = focus_levels)
  # Keep original order for plotting if factor, otherwise sort numeric
  if(is.factor(focus_data)) indices_df$level <- factor(indices_df$level, levels = levels(focus_data))


  # --- Calculate Unstandardized Index ---
  # Need log-transformed observed for geometric mean if applicable
  if (logged) {
    log_observed <- observed # Already logged
  } else if (!use_arithmetic_mean) {
    log_observed <- log(observed)
  }

  if (use_arithmetic_mean) {
      unstan_agg <- aggregate(list(unstan = observed),
                             list(level = data[[focus_col_agg]]),
                             mean, na.rm = TRUE)
      # Turn into a relative index
      mean_unstan <- mean(unstan_agg$unstan, na.rm=TRUE)
      indices_df <- merge(indices_df, unstan_agg, by = "level", sort = FALSE)
      indices_df$unstan <- indices_df$unstan / mean_unstan

  } else {
      unstan_agg <- aggregate(list(unstan = log_observed),
                              list(level = data[[focus_col_agg]]),
                              mean, na.rm = TRUE)
       # Calculate geometric mean index (relative)
      mean_log_unstan <- mean(unstan_agg$unstan, na.rm=TRUE)
      indices_df <- merge(indices_df, unstan_agg, by = "level", sort = FALSE)
      indices_df$unstan <- exp(indices_df$unstan - mean_log_unstan)
  }


  # --- Calculate Standardized Index (from full model) ---
  # Use predict.gam to get term predictions for the focus term
  # Need newdata where other variables are set to reference levels (or mean for numeric)
  # Simpler: get predictions for the focus term and average by level.

  pred_terms <- tryCatch({
      predict(model, type = "terms", se.fit = TRUE)
  }, error = function(e) {
      warning("Could not get term predictions with standard errors directly. Trying without SE. Error: ", e$message)
      # Fallback without SE if needed for basic index
      list(fit = predict(model, type = "terms"), se.fit = NULL)
  })

  if (!is.null(pred_terms$fit)) {
      focus_term_cols <- grep(paste0("^(s|te|ti|t2)?\\(?.*?", focus, ".*\\)?$"), colnames(pred_terms$fit), value=TRUE)
      if(length(focus_term_cols) == 0) {
          # Try exact match if grep failed (e.g., simple factor)
          if(focus %in% colnames(pred_terms$fit)) {
              focus_term_cols <- focus
          } else {
              stop("Could not find the focus term '", focus, "' in the predicted terms: ", paste(colnames(pred_terms$fit), collapse=", "))
          }
      }

      # Sum up contributions if focus term is involved in multiple smooths/interactions
      focus_effects_link <- rowSums(pred_terms$fit[, focus_term_cols, drop = FALSE])

      # Calculate standard errors (approximate if SEs not available per term)
      # V <- vcov(model) # Full covariance matrix
      # Need to map focus term levels to coefficients - complex for smooths!
      # Alternative: Use the SE from predict(type='terms') if available and focus term is simple factor

      focus_effects_se <- NA # Placeholder, SE calculation is complex

      # If se.fit was returned and the focus term is a simple factor
      if (!is.null(pred_terms$se.fit) && focus %in% colnames(pred_terms$se.fit)) {
         focus_effects_se <- pred_terms$se.fit[, focus]

         # Aggregate effects and SEs by focus level
         stan_agg <- aggregate(list(effect = focus_effects_link, se = focus_effects_se),
                              list(level = data[[focus_col_agg]]),
                              mean, na.rm = TRUE) # Using mean effect/SE per level

         # Center the effects (link scale) and exponentiate
         base_link <- mean(stan_agg$effect, na.rm = TRUE)
         indices_df <- merge(indices_df, stan_agg, by = "level", sort = FALSE)

         indices_df$stan <- exp(indices_df$effect - base_link)
         # Approximate CI on response scale
         indices_df$stanLower <- exp((indices_df$effect - base_link) - 2 * indices_df$se)
         indices_df$stanUpper <- exp((indices_df$effect - base_link) + 2 * indices_df$se)

      } else {
          warning("Standard errors for the focus term '", focus, "' could not be reliably extracted from predict(type='terms'). Standardized index CIs will be missing.")
           stan_agg <- aggregate(list(effect = focus_effects_link),
                              list(level = data[[focus_col_agg]]),
                              mean, na.rm = TRUE) # Using mean effect per level
          # Center the effects (link scale) and exponentiate
         base_link <- mean(stan_agg$effect, na.rm = TRUE)
         indices_df <- merge(indices_df, stan_agg, by = "level", sort = FALSE)
         indices_df$stan <- exp(indices_df$effect - base_link)
         indices_df$stanLower <- NA
         indices_df$stanUpper <- NA
      }

  } else {
      warning("Could not obtain term predictions ('predict(type=\"terms\")'). Standardized index cannot be calculated.")
      indices_df$stan <- NA
      indices_df$stanLower <- NA
      indices_df$stanUpper <- NA
  }


  # --- Iteratively Add Terms and Calculate Indices/Metrics ---
  summary_list <- list()
  current_formula <- reformulate("1", response = response_name)
  iter_terms <- c("intercept", terms) # Include intercept step

  # Initial null model (intercept only)
  null_model <- tryCatch(
       gam(current_formula, data = data, family = family(model)),
       error = function(e) {warning("Failed to fit null model: ", e$message); NULL}
     )
   if(is.null(null_model)) {
       stop("Cannot proceed without fitting the null model.")
   }

  logLikeInterceptOnly <- if (!is.null(logLik(null_model))) logLik(null_model) else NA
  null_deviance <- if(!is.null(null_model$deviance)) null_model$deviance else NA # Store null deviance
  n_obs <- nobs(model) # Use nobs from full model for consistency

  # Store fitted values for R2 calc later if needed
  # For log response, comparison is often done on log scale
  target_y <- if (logged) log(observed) else observed
  target_y_finite <- is.finite(target_y)


  for (i in seq_along(iter_terms)) {
    term_label <- iter_terms[i] # Label for this step

    if (i == 1) { # Intercept model
        step_model <- null_model
        formula_string <- "~ 1" # For index naming
    } else {
        term_to_add <- terms[i-1] # Actual model term
        # Update formula
        current_terms <- terms[1:(i-1)]
        current_formula <- reformulate(current_terms, response = response_name)
        formula_string <- paste("~", paste(current_terms, collapse = " + ")) # For index naming

        # Fit model for this step
        step_model <- tryCatch(
            gam(current_formula, data = data, family = family(model)),
             error = function(e) {warning("Failed to fit model at step ", i, " (term: ", term_to_add,"): ", e$message); NULL}
        )
    }

    if (!is.null(step_model)) {
        k_params <- length(coef(step_model)) # Number of coefficients
        logLike <- if (!is.null(logLik(step_model))) logLik(step_model) else NA
        aic_val <- if (!is.null(logLik(step_model))) AIC(step_model) else NA
        deviance <- if (!is.null(step_model$deviance)) step_model$deviance else NA

         # Calculate R-squared variants
         r2 <- NA; r2Dev <- NA; r2Negel <- NA

         # Fitted values on the link scale, then transformed if necessary
         # Use predict to ensure consistency, esp. with offset terms etc.
         fitted_vals_link <- predict(step_model, newdata=data) # Link scale by default

         # Need predictions on response scale for some R2
         # fitted_vals_response <- family(step_model)$linkinv(fitted_vals_link)

         # SS-based R2 (Correlation based - use log scale if original response was logged)
         if (logged) {
             # Compare observed log(y) vs fitted log(y) (which is link scale for log link)
             if (length(na.omit(target_y)) > 2 && length(na.omit(fitted_vals_link)) > 2) {
                 finite_idx <- target_y_finite & is.finite(fitted_vals_link)
                 if(sum(finite_idx) > 2) {
                    r2 <- cor(target_y[finite_idx], fitted_vals_link[finite_idx])^2
                 }
             }
         } else {
             # Compare observed y vs fitted y (response scale)
              fitted_vals_response <- family(step_model)$linkinv(fitted_vals_link)
              if (length(na.omit(target_y)) > 2 && length(na.omit(fitted_vals_response)) > 2) {
                  finite_idx <- target_y_finite & is.finite(fitted_vals_response)
                  if(sum(finite_idx) > 2) {
                    r2 <- cor(target_y[finite_idx], fitted_vals_response[finite_idx])^2
                  }
             }
         }

        # Deviance R2
        if (!is.na(deviance) && !is.na(null_deviance) && null_deviance > 0) {
          r2Dev <- (null_deviance - deviance) / null_deviance
        }

        # Nagelkerke R2
        if (!is.na(logLike) && !is.na(logLikeInterceptOnly) && !is.na(n_obs) && n_obs > 0) {
          # Handle potential Inf from exp()
          term1_exp <- suppressWarnings(exp((logLikeInterceptOnly - logLike) * (2 / n_obs)))
          term2_exp <- suppressWarnings(exp(logLikeInterceptOnly * (2 / n_obs)))
          if(is.finite(term1_exp) && is.finite(term2_exp) && (1 - term2_exp) != 0){
             r2Negel <- (1 - term1_exp) / (1 - term2_exp)
          } else {
             r2Negel <- NA # Assign NA if calculation resulted in Inf/NaN or division by zero
          }
        }

         # Calculate effects/index for this model step (relative to focus term)
         # Need predictions specific to the focus term from this *step_model*
         step_pred_terms <- tryCatch(
             predict(step_model, type = "terms"),
             error = function(e){ warning("Cannot get term preds for step ", i, ". Error: ", e$message); NULL}
         )

         index_col_name <- if (i == 1) "intercept" else paste("+", terms[i-1]) # Consistent naming with original

         if (!is.null(step_pred_terms)) {
             step_focus_cols <- grep(paste0("^(s|te|ti|t2)?\\(?.*?", focus, ".*\\)?$"), colnames(step_pred_terms), value=TRUE)
             if(length(step_focus_cols) == 0 && focus %in% colnames(step_pred_terms)) step_focus_cols <- focus

             if (length(step_focus_cols) > 0) {
                 step_focus_effects_link <- rowSums(step_pred_terms[, step_focus_cols, drop = FALSE])
                 step_stan_agg <- aggregate(list(effect = step_focus_effects_link),
                                            list(level = data[[focus_col_agg]]),
                                            mean, na.rm = TRUE)
                 step_base_link <- mean(step_stan_agg$effect, na.rm = TRUE)
                 step_stan_agg[[index_col_name]] <- exp(step_stan_agg$effect - step_base_link)

                 # Merge into main indices_df
                 indices_df <- merge(indices_df, step_stan_agg[, c("level", index_col_name)], by = "level", all.x = TRUE, sort = FALSE)

             } else {
                 # Focus term not in model yet, index is flat (1)
                  indices_df[[index_col_name]] <- 1.0
             }
         } else {
             # Prediction failed, cannot calculate index for this step
             indices_df[[index_col_name]] <- NA
             warning("Index calculation failed for step: ", index_col_name)
         }


    } else {
        # Model fitting failed for this step
        k_params <- NA; logLike <- NA; aic_val <- NA; r2 <- NA; r2Dev <- NA; r2Negel <- NA
        index_col_name <- if (i == 1) "intercept" else paste("+", terms[i-1])
        indices_df[[index_col_name]] <- NA # Mark index as NA for this failed step
    }

     summary_list[[i]] <- data.frame(
          term_added = term_label, # Use the step label (intercept, + term1, ...)
          term_actual = ifelse(i==1, "(Intercept)", terms[i-1]), # Underlying term added
          k = k_params,
          logLik = logLike,
          AIC = aic_val,
          R2_cor = r2,
          R2_dev = r2Dev,
          R2_nagel = r2Negel
        )

  } # End loop through terms

  summary_df <- bind_rows(summary_list)

  # Calculate differences for summary table (like original)
  summary_df <- summary_df %>%
      mutate(
          k_diff = c(NA, diff(k)), # Diff doesn't make sense for first row (intercept)
          R2_cor_diff = c(NA, diff(R2_cor)),
          R2_dev_diff = c(NA, diff(R2_dev)),
          R2_nagel_diff = c(NA, diff(R2_nagel))
      ) %>%
      select(term_added, term_actual, k, k_diff, logLik, AIC,
             R2_cor, R2_cor_diff, R2_dev, R2_dev_diff, R2_nagel, R2_nagel_diff)


  # --- Calculate Term Predictions and Influences (from full model) ---
  # Reuse predictions from earlier if available
  if (is.null(pred_terms$fit)) {
      pred_terms <- tryCatch({
          predict(model, type = "terms", se.fit = TRUE)
      }, error = function(e) {
          warning("Final predictions failed: ", e$message)
          NULL
      })
  }

  influences_df <- data.frame(level = focus_levels)
  if(is.factor(focus_data)) influences_df$level <- factor(influences_df$level, levels = levels(focus_data))

  overall_influ <- numeric(length(terms))
  trend_influ <- numeric(length(terms))
  names(overall_influ) <- terms
  names(trend_influ) <- terms

  if (!is.null(pred_terms$fit)) {
      pred_df <- as.data.frame(pred_terms$fit)
      pred_se_df <- if(!is.null(pred_terms$se.fit)) as.data.frame(pred_terms$se.fit) else NULL

      # Combine data with predictions
      # Ensure row alignment if model frame was different from original data due to NAs
      if(nrow(pred_df) == nrow(data)) {
          full_pred_data <- bind_cols(data, pred_df)
          if(!is.null(pred_se_df)) {
              names(pred_se_df) <- paste0(names(pred_df), "_se")
              full_pred_data <- bind_cols(full_pred_data, pred_se_df)
          }
      } else if (!is.null(model$na.action)) {
           # Handle NAs - predict only returns values for non-NA rows
           warning("NA actions detected. Aligning predictions with original data. Check results carefully.")
           # Pad predictions with NAs
           na_rows <- as.integer(model$na.action)
           padded_pred_df <- matrix(NA, nrow = nrow(data), ncol = ncol(pred_df))
           padded_pred_df[-na_rows, ] <- as.matrix(pred_df)
           padded_pred_df <- as.data.frame(padded_pred_df)
           names(padded_pred_df) <- names(pred_df)
           full_pred_data <- bind_cols(data, padded_pred_df)

           if(!is.null(pred_se_df)) {
               padded_pred_se_df <- matrix(NA, nrow = nrow(data), ncol = ncol(pred_se_df))
               padded_pred_se_df[-na_rows, ] <- as.matrix(pred_se_df)
               padded_pred_se_df <- as.data.frame(padded_pred_se_df)
               names(padded_pred_se_df) <- paste0(names(pred_df), "_se")
               full_pred_data <- bind_cols(full_pred_data, padded_pred_se_df)
           }
       } else {
          warning("Prediction rows (", nrow(pred_df), ") do not match data rows (", nrow(data), ") and no NA action found. Influence calculations may be incorrect.")
          # Attempt merge based on row names or identifiers if possible, otherwise skip
          full_pred_data <- data # Placeholder, influences will fail
       }


      x$predictions <- full_pred_data # Store predictions with data

      # Calculate influences per term (excluding focus term itself)
      focus_term_pattern <- paste0("^(s|te|ti|t2)?\\(?.*?", focus, ".*\\)?$") # Pattern to identify focus term cols

      for (term in terms) {
          # Identify prediction columns related to this term
          term_cols <- grep(paste0("^(s|te|ti|t2)?\\(?.*?", term, ".*\\)?$"), colnames(pred_df), value=TRUE)
          if(length(term_cols) == 0 && term %in% colnames(pred_df)) term_cols <- term # Exact match for simple factors

           # Exclude columns that *also* match the focus term if term != focus
           # (handles interactions like s(year, by=fac) when focus='year', term='fac')
           if (term != focus) {
                term_cols <- term_cols[!grepl(focus_term_pattern, term_cols)]
           }


          if (length(term_cols) > 0 && term != focus) {
              # Sum contributions if term is split (e.g., interaction, tensor smooth)
              term_effect_link <- rowSums(full_pred_data[, term_cols, drop = FALSE], na.rm = TRUE)

              # Aggregate influence by focus level
              infl_agg <- aggregate(list(value = term_effect_link),
                                   list(level = full_pred_data[[focus_col_agg]]),
                                   mean, na.rm = TRUE)
              names(infl_agg)[names(infl_agg)=="value"] <- term # Rename col to term name

              # Calculate overall and trend influence statistics (on link scale then exp)
              if (nrow(infl_agg) > 1) {
                 overall_influ[term] <- exp(mean(abs(infl_agg[[term]]), na.rm = TRUE)) - 1

                 # Trend influence (simple linear trend over ordered levels)
                 # Requires numeric representation of levels
                 if (is.factor(infl_agg$level)) {
                    level_num <- as.numeric(infl_agg$level)
                 } else {
                    level_num <- infl_agg$level # Assume already numeric or comparable
                 }

                 valid_trend_data <- is.finite(level_num) & is.finite(infl_agg[[term]])
                 if(sum(valid_trend_data) > 1) {
                    trend_cov <- cov(level_num[valid_trend_data], infl_agg[[term]][valid_trend_data])
                    trend_var <- var(level_num[valid_trend_data])
                    if (is.finite(trend_cov) && is.finite(trend_var) && trend_var != 0) {
                       trend_influ[term] <- exp(trend_cov / trend_var) - 1
                    } else {
                        trend_influ[term] <- NA # Assign NA if cov/var computation failed
                    }
                 } else {
                      trend_influ[term] <- NA # Assign NA if not enough data for trend
                 }
              } else {
                 overall_influ[term] <- NA
                 trend_influ[term] <- NA
              }

              # Merge into main influences_df
              influences_df <- merge(influences_df, infl_agg[, c("level", term)], by = "level", all.x = TRUE, sort=FALSE)

          } else if (term == focus) {
              # Add NA column for focus term itself
              influences_df[[term]] <- NA
              overall_influ[term] <- NA
              trend_influ[term] <- NA
          } else {
              # Term not found in predictions or is focus term
              warning("Influence calculation skipped for term: ", term)
               overall_influ[term] <- NA
               trend_influ[term] <- NA
          }
      } # End loop terms for influence

  } else {
      warning("Cannot calculate influences due to prediction failure.")
      x$predictions <- data # Store original data only
  }

  # Add influence stats to summary_df (match by actual term name)
  influ_stats <- data.frame(term_actual = terms, overall = overall_influ, trend = trend_influ)

  # Add NA row for intercept to allow merging
  influ_stats <- bind_rows(data.frame(term_actual="(Intercept)", overall=NA, trend=NA), influ_stats)

  # Ensure term_actual is character for merging
  summary_df$term_actual <- as.character(summary_df$term_actual)
  influ_stats$term_actual <- as.character(influ_stats$term_actual)


  summary_df <- left_join(summary_df, influ_stats, by = "term_actual")

  # --- Store Results ---
  x$indices <- indices_df
  x$summary <- summary_df
  x$influences <- influences_df

  invisible(x) # Return updated object invisibly
}


# --- Plotting Functions ---

#' Standardization Plot (ggplot2)
#'
#' Compares unstandardized and standardized indices for the focus term.
#'
#' @param x An object of class 'influence_gam' where `calculate()` has been run.
#' @param show_ci Logical, whether to show confidence intervals for the standardized index.
#' @return A ggplot object.
#' @export
plot_stan <- function(x, show_ci = TRUE) {
  if (is.null(x$indices)) stop("Run calculate() first.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")

  df <- x$indices
  focus_lab <- x$labels[[x$focus]] %||% x$focus # Use label or name
  is_numeric_focus <- is.numeric(df$level)

  # Need to pivot for ggplot
  df_long <- df %>%
    select(level, unstan, stan, stanLower, stanUpper) %>%
    pivot_longer(cols = c("unstan", "stan"), names_to = "Index Type", values_to = "Index Value") %>%
    mutate(`Index Type` = ifelse(`Index Type` == "stan", "Standardized", "Unstandardized")) %>%
    # Manually add CI info back for standardized points
    left_join(select(df, level, stanLower, stanUpper), by = "level") %>%
    # Ensure level retains its original type (factor or numeric)
    mutate(level = if(is.factor(df$level)) factor(level, levels=levels(df$level)) else as.numeric(level))


  gg <- ggplot(df_long, aes(x = level, y = `Index Value`)) +
    geom_line(aes(group = `Index Type`, color = `Index Type`)) +
    geom_point(aes(color = `Index Type`, shape = `Index Type`), size = 2.5)

  if (show_ci && "stanLower" %in% names(df_long) && "stanUpper" %in% names(df_long)) {
    # Add CI only for Standardized points
    ci_data <- filter(df_long, `Index Type` == "Standardized")
    # Check if CI data is available and finite
    if(nrow(ci_data) > 0 && any(is.finite(ci_data$stanLower) & is.finite(ci_data$stanUpper))) {
       gg <- gg + geom_errorbar(data = ci_data,
                              aes(ymin = stanLower, ymax = stanUpper),
                              width = 0.1, color = "black") # Adjust color/width as needed
    } else {
       warning("Confidence interval data for standardized index is missing or non-finite.")
    }
  }

  gg <- gg +
    scale_color_manual(values = c("Standardized" = "black", "Unstandardized" = "grey50")) +
    scale_shape_manual(values = c("Standardized" = 16, "Unstandardized" = 1)) +
    labs(
      title = "Standardization Plot",
      subtitle = paste("Focus Term:", focus_lab),
      x = focus_lab,
      y = "Relative Index",
      color = "Index Type",
      shape = "Index Type"
    ) +
    theme_bw() +
    theme(legend.position = "top")

  # Handle axis for factor vs numeric focus
  if (!is_numeric_focus) {
    gg <- gg + scale_x_discrete() # Ensures factor order is respected
  }

  return(gg)
}


#' Step Plot (ggplot2)
#'
#' Shows how the standardized index for the focus term changes as each model
#' term is added sequentially.
#'
#' @param x An object of class 'influence_gam' where `calculate()` has been run.
#' @param panels Logical, if TRUE (default), creates faceted plot. If FALSE,
#'               plots all steps on a single panel.
#' @param base_color Color for the current step's line/points.
#' @param previous_color Color for the previous step's line (in panel view).
#' @param other_color Color for other preceding steps (in panel view).
#' @return A ggplot object.
#' @export
plot_step <- function(x, panels = TRUE, base_color = "blue",
                     previous_color = "darkgrey", other_color = "lightgrey") {
  if (is.null(x$indices)) stop("Run calculate() first.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required.")

  df <- x$indices
  focus_lab <- x$labels[[x$focus]] %||% x$focus
  is_numeric_focus <- is.numeric(df$level)

  # Identify index columns (start after stanUpper usually)
  index_cols <- setdiff(names(df), c("level", "unstan", "effect", "se", "stan", "stanLower", "stanUpper"))
  if (length(index_cols) == 0) stop("No step index columns found in indices dataframe.")

  # Pivot longer
  df_long <- df %>%
    select(level, all_of(index_cols)) %>%
    pivot_longer(cols = all_of(index_cols), names_to = "Step", values_to = "Index") %>%
    # Ensure level retains its original type
    mutate(level = if(is.factor(df$level)) factor(level, levels=levels(df$level)) else as.numeric(level))

  # Ensure steps are ordered correctly (based on column order in original df)
  df_long$Step <- factor(df_long$Step, levels = index_cols)

  if (panels) {
    # Create data for previous/other steps for faceting
     df_panel <- df_long %>%
       mutate(step_idx = as.integer(Step)) %>%
       group_by(level) %>%
       # Lagged index for previous step
       mutate(Index_prev = lag(Index, order_by = step_idx),
             Step_prev = lag(Step, order_by = step_idx)) %>%
       ungroup()

    # Base plot for panels
    gg <- ggplot(df_panel, aes(x = level, y = Index)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey30") +
      # Lines for steps *before* the previous one (greyed out)
      geom_line(data = . %>% mutate(current_step_idx = as.integer(Step)) %>%
                       tidyr::crossing(prev_step_idx = unique(df_panel$step_idx)) %>% # Show all previous lines on each panel
                       filter(prev_step_idx < current_step_idx -1) %>%
                       left_join(df_long %>% rename(Index_hist=Index, Step_hist=Step) %>% mutate(step_idx=as.integer(Step_hist)), by=c("level", "prev_step_idx"="step_idx"))
                       , aes(y = Index_hist, group = prev_step_idx), color = other_color, linewidth = 0.8) +

      # Line for the immediately preceding step (dashed)
      geom_line(aes(y = Index_prev, group = 1), color = previous_color, linetype = "dashed", linewidth = 1) +
      # Line and points for the current step
      geom_line(aes(group = 1), color = base_color, linewidth = 1.2) +
      geom_point(color = base_color, size = 2.5) +
      facet_wrap(~Step, ncol = 1, strip.position = "left") + # Use strip.position
      labs(
        title = "Step Plot: Index Evolution by Term Addition",
        subtitle = paste("Focus Term:", focus_lab),
        x = focus_lab,
        y = "Relative Index" # Facet strips indicate the step
      ) +
      theme_bw() +
      theme(strip.background = element_blank(),
            strip.placement = "outside", # Ensures labels are outside plot area
            strip.text.y.left = element_text(angle = 0, hjust=1)) # Adjust strip text

     # Handle axis for factor vs numeric focus
      if (!is_numeric_focus) {
        gg <- gg + scale_x_discrete()
      }


  } else {
    # Single panel plot
    gg <- ggplot(df_long, aes(x = level, y = Index, color = Step, shape = Step)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey30") +
      geom_line(aes(group = Step), linewidth=0.8) +
      geom_point(size = 2.5) +
      scale_color_viridis_d() + # Or other discrete scale
      scale_shape_manual(values = seq_along(levels(df_long$Step))) + # Different shapes
      labs(
        title = "Step Plot: Index Evolution by Term Addition",
        subtitle = paste("Focus Term:", focus_lab),
        x = focus_lab,
        y = "Relative Index",
        color = "Step Added",
        shape = "Step Added"
      ) +
      theme_bw() +
      theme(legend.position = "top")

     # Handle axis for factor vs numeric focus
      if (!is_numeric_focus) {
        gg <- gg + scale_x_discrete()
      }
  }

  return(gg)
}


#' Influence Plot (ggplot2)
#'
#' Shows the influence of each non-focus term on the focus term index across
#' the levels of the focus term. Influence is plotted on the response scale (exp(effect)).
#'
#' @param x An object of class 'influence_gam' where `calculate()` has been run.
#' @param panels Logical, if TRUE (default), creates faceted plot per influencing term.
#'               If FALSE, plots all influences on a single panel.
#' @param base_color Color for lines/points.
#' @return A ggplot object.
#' @export
plot_influ <- function(x, panels = TRUE, base_color = "red") {
  if (is.null(x$influences)) stop("Run calculate() or ensure predictions were successful first.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required.")

  df <- x$influences
  focus_lab <- x$labels[[x$focus]] %||% x$focus
  is_numeric_focus <- is.numeric(df$level)

  # Identify influence columns (all except 'level' and maybe the focus term itself)
  influ_cols <- setdiff(names(df), c("level", x$focus)) # Exclude level and focus col
  if (length(influ_cols) == 0) stop("No influence columns found (excluding focus term).")

  # Pivot longer
  df_long <- df %>%
    select(level, all_of(influ_cols)) %>%
    pivot_longer(cols = all_of(influ_cols), names_to = "Influencing Term", values_to = "Influence Link Scale") %>%
    # Transform to response scale (exponentiate)
    mutate(Influence = exp(`Influence Link Scale`)) %>%
    # Ensure level retains its original type
    mutate(level = if(is.factor(df$level)) factor(level, levels=levels(df$level)) else as.numeric(level)) %>%
     # Use custom labels if available
     mutate(`Influencing Term` = factor(
         `Influencing Term`,
         levels = influ_cols, # Original order
         labels = sapply(influ_cols, function(term) x$labels[[term]] %||% term)
     ))

  if (panels) {
    gg <- ggplot(df_long, aes(x = level, y = Influence)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey30") +
      geom_line(aes(group = 1), color = base_color, linewidth = 1) +
      geom_point(color = base_color, size = 2.5) +
      facet_wrap(~`Influencing Term`, ncol = 1, strip.position = "left") +
      labs(
        title = "Influence Plot: Effect of Terms on Focus Term Index",
        subtitle = paste("Focus Term:", focus_lab),
        x = focus_lab,
        y = "Influence Multiplier" # Facet strips indicate the term
      ) +
      scale_y_log10() + # Often useful for multiplicative influence
      theme_bw() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            strip.text.y.left = element_text(angle = 0, hjust=1))

      # Handle axis for factor vs numeric focus
      if (!is_numeric_focus) {
        gg <- gg + scale_x_discrete()
      }

  } else {
    gg <- ggplot(df_long, aes(x = level, y = Influence, color = `Influencing Term`, shape = `Influencing Term`)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey30") +
      geom_line(aes(group = `Influencing Term`), linewidth = 0.8) +
      geom_point(size = 2.5) +
      scale_color_viridis_d() + # Or other discrete scale
      scale_shape_manual(values = seq_along(levels(df_long$`Influencing Term`))) +
      scale_y_log10() + # Often useful for multiplicative influence
      labs(
        title = "Influence Plot: Effect of Terms on Focus Term Index",
        subtitle = paste("Focus Term:", focus_lab),
        x = focus_lab,
        y = "Influence Multiplier",
        color = "Influencing Term",
        shape = "Influencing Term"
      ) +
      theme_bw() +
      theme(legend.position = "top")

      # Handle axis for factor vs numeric focus
      if (!is_numeric_focus) {
        gg <- gg + scale_x_discrete()
      }
  }

  return(gg)
}


#' Combined Step and Influence Plot (ggplot2)
#'
#' Arranges the faceted step plot and faceted influence plot side-by-side.
#' Requires the 'patchwork' package.
#'
#' @param x An object of class 'influence_gam' where `calculate()` has been run.
#' @param ... Arguments passed to `plot_step` and `plot_influ`.
#' @return A patchwork object combining the plots.
#' @export
plot_step_influ <- function(x, ...) {
   if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Please install the 'patchwork' package to use plot_step_influ.")
   }

   # Generate faceted plots
   p_step <- plot_step(x, panels = TRUE, ...)
   p_influ <- plot_influ(x, panels = TRUE, ...)

   # Combine plots
   # Note: Alignment might not be perfect if number of terms differs
   # or if focus term included in influence plot implicitly.
   # Adjust layout as needed.

   # Remove y-axis title from influence plot for cleaner look
   p_influ <- p_influ + theme(axis.title.y = element_blank(),
                             axis.text.y = element_blank(),
                             axis.ticks.y = element_blank(),
                             strip.text.y.left = element_text(angle = 0, hjust=0) # Adjust label position
                             ) +
                        scale_y_log10(labels=scales::label_log()) # Keep log scale but maybe show labels

   # Adjust strip position if needed
   # p_influ <- p_influ + theme(strip.position="right")


   combined_plot <- patchwork::wrap_plots(p_step, p_influ, ncol = 2) +
                    patchwork::plot_layout(guides = 'collect') & # Collect legends if any
                    theme(legend.position = 'top')

   print(combined_plot) # Ensure plot is displayed
   invisible(combined_plot)
}


#' Coefficient-Distribution-Influence (CDI) Plot (ggplot2)
#'
#' Creates a three-panel plot for a specific model term showing:
#' 1. Its coefficient/effect across its levels/values.
#' 2. The distribution of the focus term across the levels/values of this term.
#' 3. The influence of this term on the focus term index.
#' Requires the 'patchwork' package.
#'
#' @param x An object of class 'influence_gam' where `calculate()` has been run.
#' @param term The model term (character string, e.g., "s(x0)", "fac") for which
#'             to generate the CDI plot.
#' @param variable Optional: The specific variable name in the data corresponding
#'                 to the term, needed if it cannot be automatically inferred
#'                 (e.g., for `log(var)` term, variable is `"var"`).
#' @param bin_numeric Logical or Integer: If the term's variable is numeric,
#'                    should it be binned for plotting? If TRUE, uses `pretty`
#'                    breaks. If an integer, uses that many bins (`cut_number`).
#'                    Default is TRUE. Set to FALSE to plot against raw numeric values.
#' @param base_color Color used for points/lines in the plots.
#' @return A patchwork object combining the three plots.
#' @export
plot_cdi <- function(x, term, variable = NULL, bin_numeric = TRUE, base_color = "purple") {
  if (is.null(x$predictions) || is.null(x$influences)) {
      stop("Run calculate() first and ensure predictions/influences were successful.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork package required for CDI plots.")
  }
  if (!requireNamespace("rlang", quietly = TRUE)) {
      stop("rlang package required for CDI plots.") # For :=
  }

  preds <- x$predictions
  influences <- x$influences
  focus_var <- x$focus
  focus_lab <- x$labels[[focus_var]] %||% focus_var
  term_lab <- x$labels[[term]] %||% term
  term_order <- x$orders[[term]] %||% "asis"

  # --- Identify Variable and Term Columns ---
  term_cols_fit <- grep(paste0("^", gsub("\\(", "\\\\(", gsub("\\)", "\\\\)", term)), "_se$"), names(preds), invert = TRUE, value = TRUE) # Avoid SE cols first
  term_cols_fit <- grep(paste0("^(s|te|ti|t2)?\\(?.*?", term, ".*\\)?$"), term_cols_fit, value=TRUE) # Filter for the term

  # Try to find exact match if grep failed
  if(length(term_cols_fit) == 0 && term %in% names(preds)) term_cols_fit <- term

  if (length(term_cols_fit) == 0) stop("Could not find prediction column for term: ", term)
  # For simplicity, sum if multiple columns match (e.g., tensor smooth components)
  preds$.term_effect <- rowSums(preds[, term_cols_fit, drop = FALSE], na.rm = TRUE)

  # Find SE column if exists
  term_col_se <- paste0(term_cols_fit[1], "_se") # Assume SE matches first fit col name convention
  if(length(term_cols_fit) == 1 && term_col_se %in% names(preds)) {
      preds$.term_se <- preds[[term_col_se]]
  } else {
      preds$.term_se <- NA
      if(length(term_cols_fit) == 1) warning("SE column '", term_col_se, "' not found for term '", term, "'. CIs omitted.")
      if(length(term_cols_fit) > 1) warning("SEs for multi-column terms like '", term, "' not directly supported for CDI plot CIs. CIs omitted.")
  }


  # --- Identify the Variable associated with the term ---
  if (is.null(variable)) {
    # Attempt to extract variable name(s) from the term string
    term_vars <- all.vars(stats::reformulate(term)[[2]]) # Extract vars from RHS of dummy formula
    # Heuristic: if only one var, use it. If multiple, prefer one present in data.
    if (length(term_vars) == 1 && term_vars %in% names(preds)) {
        variable <- term_vars
    } else {
        possible_vars <- term_vars[term_vars %in% names(preds)]
        if(length(possible_vars) == 1) {
            variable <- possible_vars
        } else if (length(possible_vars) > 1) {
            # If multiple variables from term are in data (e.g., interaction), which one to plot against?
            # Defaulting to the first one, but user might need to specify.
            variable <- possible_vars[1]
             warning("Term '", term, "' involves multiple variables found in data: ",
                     paste(possible_vars, collapse=", "), ". Plotting against '", variable,
                     "'. Specify 'variable' argument if different axis is desired.")
        } else {
             stop("Could not automatically determine the variable for term '", term,
                  "'. Please specify using the 'variable' argument.")
        }
    }
  }
  if (!variable %in% names(preds)) stop("Specified variable '", variable, "' not found in prediction data.")

  var_lab <- x$labels[[variable]] %||% variable
  preds$.term_variable <- preds[[variable]] # Use specific name for clarity


  # --- Handle Numeric Variable Binning ---
  is_numeric_term <- is.numeric(preds$.term_variable)
  if (is_numeric_term && !isFALSE(bin_numeric)) {
    if (isTRUE(bin_numeric)) {
        # Use pretty breaks
        breaks <- pretty(preds$.term_variable, n = 20) # More breaks than original
        preds$.term_variable_bin <- cut(preds$.term_variable, breaks = breaks, include.lowest = TRUE)
        plot_var <- ".term_variable_bin" # Use binned version for plotting axes
        term_axis_lab <- paste(var_lab, "(Binned)")
    } else if (is.numeric(bin_numeric) && bin_numeric > 0) {
        # Use cut_number for specified bins
        preds$.term_variable_bin <- cut_number(preds$.term_variable, n = bin_numeric)
        plot_var <- ".term_variable_bin"
        term_axis_lab <- paste(var_lab, "(Binned)")
    } else {
        warning("Invalid 'bin_numeric' value. Plotting against raw numeric values.")
        plot_var <- ".term_variable" # Plot raw numeric
        term_axis_lab <- var_lab
    }

  } else {
     # Categorical variable or no binning requested
     preds$.term_variable_bin <- factor(preds$.term_variable) # Ensure factor if not numeric
     plot_var <- ".term_variable_bin"
     term_axis_lab <- var_lab
  }

  # --- Aggregate for Plots ---
  # Use the binned/factored variable for aggregation
  cdi_agg <- preds %>%
    group_by(!!sym(plot_var)) %>%
    summarise(
      coeff = mean(.term_effect, na.rm = TRUE),
      se = if (all(is.na(.term_se))) NA else mean(.term_se, na.rm = TRUE), # Mean SE if available
      .groups = 'drop'
    ) %>%
    mutate(
      lower = coeff - 1.96 * se,
      upper = coeff + 1.96 * se
    ) %>%
    # Ensure levels are ordered correctly (asis or by coeff)
    mutate(!!sym(plot_var) := factor(!!sym(plot_var), levels = if(term_order == "coef") {
                                                              levels(!!sym(plot_var))[order(coeff)]
                                                          } else {
                                                              levels(!!sym(plot_var)) # Default factor order
                                                          }
                                                          ))

  # Ensure plot_var exists after summarise
  if (!plot_var %in% names(cdi_agg)) {
     stop("Plotting variable '", plot_var, "' lost during aggregation.")
  }

  # --- 1. Coefficient Plot ---
  p_coeff <- ggplot(cdi_agg, aes(x = !!sym(plot_var), y = exp(coeff))) + # Exponentiate effect
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_point(color = base_color, size = 2.5) +
    scale_y_log10() + # Use log scale for multiplicative effects
    labs(
        #title = paste("Term Effect:", term_lab),
        x = NULL, # Remove x-axis label/ticks
        y = "Coefficient\n(Effect Multiplier)" # Add line break
    ) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 5, r = 5, b = 0, l = 5)) # Adjust margins

  if (any(!is.na(cdi_agg$se))) {
    p_coeff <- p_coeff + geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper)), width = 0.2, color = base_color)
  }


  # --- 2. Distribution Plot ---
  # Aggregate counts of focus variable within each bin of term variable
  distr_agg <- preds %>%
    # Ensure focus variable is factor for discrete axis
    mutate(.focus_factor = factor(!!sym(focus_var))) %>%
    count(!!sym(plot_var), .focus_factor, name = "count") %>%
    group_by(.focus_factor) %>%
    mutate(total = sum(count), prop = count / total) %>%
    ungroup() %>%
     # Ensure factor levels match coefficient plot
     mutate(!!sym(plot_var) := factor(!!sym(plot_var), levels = levels(cdi_agg[[plot_var]])))


  p_distr <- ggplot(distr_agg, aes(x = !!sym(plot_var), y = .focus_factor)) +
    geom_point(aes(size = prop), shape = 16, color = base_color, show.legend = TRUE) +
    scale_size_area(max_size = 10, name = "Proportion within\nFocus Level") + # Adjust max_size
    labs(
      x = term_axis_lab, # Label for the term's variable
      y = focus_lab # Label for the focus variable
    ) +
    theme_bw() +
    theme(legend.position = "right",
         axis.text.x = element_text(angle = 45, hjust = 1, size=9), # Rotate if needed
         plot.margin = margin(t = 0, r = 5, b = 5, l = 5))


  # --- 3. Influence Plot ---
  # Merge influence values with aggregated coefficient data levels for alignment
  influ_data <- influences %>%
     select(level, !!sym(term)) %>% # Select focus level and the term's influence column
     filter(!is.na(!!sym(term))) %>% # Remove NA influences if any
     # Ensure focus level is factor for y-axis
     mutate(level = factor(level, levels=levels(preds$.focus_factor))) %>% # Use focus factor levels from distr plot
     rename(InfluenceLink = !!sym(term)) %>%
     mutate(Influence = exp(InfluenceLink))

  p_influ <- ggplot(influ_data, aes(x = Influence, y = level)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_point(color = base_color, size = 2.5) +
    scale_x_log10() +
    labs(
        x = "Influence\nMultiplier",
        y = NULL # Remove y-axis label/ticks
    ) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(t = 0, r = 5, b = 5, l = 0)) # Adjust margins


  # --- Combine Plots using Patchwork ---
  # Layout: Coeff plot top, Distr bottom-left, Influ bottom-right
  layout <- "
  AAAA
  B#CC
  "
  # Simpler layout: Coeff above, Dist + Influ below
  layout2 <- "
  AA
  BC
  "
  # Layout like original: Coeff wide top, Distr left, Influ right narrow
  layout3 <- "
  AAA#
  BDDC
  BDDC
  "
  # Try layout 2
  combined_plot <- p_coeff + p_distr + p_influ +
                    patchwork::plot_layout(ncol = 2, nrow=2, design=layout2, widths=c(1, 0.6), heights=c(0.6, 1)) + # Adjust relative widths/heights
                    patchwork::plot_annotation(title = paste("CDI Plot for Term:", term_lab))


  print(combined_plot)
  invisible(combined_plot)
}

# Helper for default values (like %||% from rlang but base R)
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

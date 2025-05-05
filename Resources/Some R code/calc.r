# Updated calc function to handle splines package smoothers (ns, bs, etc.)
calc <- function(obj, islog = NULL) {
  # Get observed values
  if (inherits(obj$model, "gam")) {
    observed <- obj$model$model[[obj$response]]
  } else {
    observed <- obj$model$model[, obj$response]
  }
  
  if (inherits(observed, "Surv")) {
    observed <- as.numeric(observed[, 1])
  }
  
  if (is.null(islog)) {
    logged <- startsWith(obj$response, "log(")
  } else {
    logged <- islog
  }
  
  # Extract focus variable levels
  if (inherits(obj$model, "gam")) {
    focus_var <- obj$model$model[[obj$focus]]
  } else {
    focus_var <- obj$model$model[, obj$focus]
  }
  
  # Ensure focus_var is a factor
  if (!is.factor(focus_var)) {
    focus_var <- factor(focus_var)
  }
  
  # Create indices dataframe
  obj$indices <- data.frame(level = levels(focus_var))
  
  # Add unStandardised index
  if (logged || sum(observed <= 0) == 0) {
    if (logged) {
      log_observed <- observed
    } else {
      log_observed <- log(observed)
    }
    # Calculate geometric mean
    unstan_df <- aggregate(list(unstan = log_observed), list(level = focus_var), mean)
    obj$indices <- merge(obj$indices, unstan_df)
    # Turn into relative index
    obj$indices$unstan <- with(obj$indices, exp(unstan - mean(unstan)))
  } else {
    # Use arithmetic mean instead
    unstan_df <- aggregate(list(unstan = observed), list(level = focus_var), mean)
    obj$indices <- merge(obj$indices, unstan_df)
    # Turn into relative index
    obj$indices$unstan <- with(obj$indices, unstan / mean(unstan))
  }
  
  # Add Standardised index
  coeffs <- get_coeffs(obj)
  ses <- get_ses(obj)
  base <- mean(coeffs)
  obj$indices$stan <- exp(coeffs - base)
  obj$indices$stanLower <- exp(coeffs - base - 2 * ses)
  obj$indices$stanUpper <- exp(coeffs - base + 2 * ses)
  
  # Create models with terms successively added
  obj$summary <- NULL
  
  # Calculate influence statistics
  for (termCount in 0:length(obj$terms)) {
    if (termCount > 0) {
      term <- obj$terms[termCount]
      
      # Check if the term is a spline smoother from splines package
      is_spline_smoother <- term %in% names(obj$smoothers) && 
                           obj$smoothers[[term]]$type %in% c("ns", "bs", "poly", "lo", "pspline", "rcs")
      
      # Handle models with interactions and spline smoothers differently
      if (inherits(obj$model, "gam") && any(sapply(obj$model$smooth, function(s) s$label == term))) {
        # For smooth terms in GAMs
        terms_to_include <- obj$terms[1:termCount]
        smooth_terms <- sapply(obj$model$smooth, function(s) s$label)
        smooth_to_include <- smooth_terms[smooth_terms %in% terms_to_include]
        
        # Build formula for GAM, preserving interactions
        formula_terms <- character(0)
        for (t in terms_to_include) {
          if (t %in% smooth_to_include) {
            # Find the original specification of the smooth term
            for (s in obj$model$smooth) {
              if (s$label == t) {
                # Handle factor year in smooth terms
                if ("year" %in% s$term && is.factor(obj$data[["year"]])) {
                  # For factor year in smooths, convert to parametric
                  if (inherits(s, "tensor.smooth") || inherits(s, "ti")) {
                    # For interactions with year as factor, modify the term
                    vars <- s$term
                    year_idx <- which(vars == "year")
                    other_vars <- vars[-year_idx]
                    
                    if (length(other_vars) > 0) {
                      # For interactions with other variables
                      other_var_terms <- paste0("s(", other_vars, ")")
                      year_term <- "year"
                      # Create interaction between year and smooth
                      formula_terms <- c(formula_terms, 
                                         other_var_terms,
                                         year_term,
                                         paste0(other_var_terms, ":year"))
                    } else {
                      # Just year by itself
                      formula_terms <- c(formula_terms, "year")
                    }
                  } else {
                    # Single term with year as factor
                    formula_terms <- c(formula_terms, "year")
                  }
                } else {
                  # Regular smooth terms without factor year
                  if (inherits(s, "tensor.smooth") || inherits(s, "ti")) {
                    # For tensor interaction smooth terms (ti)
                    vars <- paste(s$term, collapse = ",")
                    formula_terms <- c(formula_terms, paste0("ti(", vars, ")"))
                  } else {
                    # For regular smooth terms (s)
                    formula_terms <- c(formula_terms, paste0("s(", s$term, ")"))
                  }
                }
                break
              }
            }
          } else if (t %in% names(obj$interactions) && obj$interactions[[t]]$type == "parametric") {
            # For parametric interactions
            formula_terms <- c(formula_terms, t)
          } else if (is_spline_smoother && t == term) {
            # For splines package smoothers (ns, bs, etc.)
            var_name <- obj$smoothers[[t]]$vars
            smoother_type <- obj$smoothers[[t]]$type
            
            # Extract parameters from the term
            term_params <- ""
            param_match <- regexpr("\\([^)]*,([^)]+)\\)", t)
            if (param_match > 0) {
              match_text <- regmatches(t, param_match)
              param_part <- sub("^\\([^,]*,", "", match_text)
              param_part <- sub("\\)$", "", param_part)
              term_params <- paste0(", ", param_part)
            }
            
            # Add the smoother term to the formula
            formula_terms <- c(formula_terms, paste0(smoother_type, "(", var_name, term_params, ")"))
          } else {
            # For regular parametric terms
            formula_terms <- c(formula_terms, t)
          }
        }
        
        # Create the updated formula
        formula_rhs <- paste(formula_terms, collapse = "+")
        formula_str <- paste(obj$response, "~", formula_rhs)
        new_formula <- as.formula(formula_str)
        
        # Update model with new formula
        model <- update(obj$model, formula = new_formula)
      } else {
        # For GLMs and parametric terms, preserving interactions
        terms_to_include <- obj$terms[1:termCount]
        
        # Check for interaction terms involving any of the included terms
        interaction_terms <- names(obj$interactions)[sapply(obj$interactions, function(i) {
          any(i$vars %in% terms_to_include)
        })]
        
        # Only include interaction terms if all their variables are in terms_to_include
        interaction_terms_to_include <- character(0)
        for (i_term in interaction_terms) {
          if (all(obj$interactions[[i_term]]$vars %in% terms_to_include)) {
            interaction_terms_to_include <- c(interaction_terms_to_include, i_term)
          }
        }
        
        # Handle spline smoothers from splines package
        spline_terms_to_include <- character(0)
        for (t in terms_to_include) {
          if (t %in% names(obj$smoothers) && 
              obj$smoothers[[t]]$type %in% c("ns", "bs", "poly", "lo", "pspline", "rcs")) {
            
            # Extract smoother information
            var_name <- obj$smoothers[[t]]$vars
            smoother_type <- obj$smoothers[[t]]$type
            
            # Extract parameters from the term
            term_params <- ""
            param_match <- regexpr("\\([^)]*,([^)]+)\\)", t)
            if (param_match > 0) {
              match_text <- regmatches(t, param_match)
              param_part <- sub("^\\([^,]*,", "", match_text)
              param_part <- sub("\\)$", "", param_part)
              term_params <- paste0(", ", param_part)
            }
            
            # Create term string
            spline_term <- paste0(smoother_type, "(", var_name, term_params, ")")
            spline_terms_to_include <- c(spline_terms_to_include, spline_term)
          }
        }
        
        # Combine main effects and allowed interaction terms
        regular_terms <- terms_to_include[!terms_to_include %in% c(names(obj$interactions), names(obj$smoothers))]
        all_terms_to_include <- unique(c(
          regular_terms,
          interaction_terms_to_include,
          spline_terms_to_include
        ))
        
        # Create formula
        formula_rhs <- paste(all_terms_to_include, collapse = "+")
        formula_str <- paste(obj$response, "~", formula_rhs)
        new_formula <- as.formula(formula_str)
        
        # Update model
        model <- update(obj$model, formula = new_formula)
      }
      
      # Get index for this model
      index <- get_effects(obj, model)
      
      # Add column to indices
      obj$indices <- cbind(obj$indices, index)
      
      # Give index column the right hand side of formula as name
      names(obj$indices)[ncol(obj$indices)] <- if (termCount == 1) term else paste("+", term)
    } else {
      term <- "intercept"
      
      # Create intercept-only model
      if (inherits(obj$model, "gam")) {
        formula_str <- paste(obj$response, "~ 1")
        model <- gam(as.formula(formula_str), data = obj$data, family = obj$model$family)
      } else {
        model <- update(obj$model, . ~ 1)
      }
    }
    
    # Extract model statistics
    if (inherits(model, "survreg")) {
      logLike <- model$loglik[2]
      fitted <- predict(model, type = "response")
    } else {
      logLike <- as.numeric(logLik(model))
      fitted <- fitted(model)
    }
    
    # Calculate R-squared values
    if (termCount == 0) {
      r2 <- 0
    } else {
      if (!logged) {
        r2 <- cor(log(observed), log(fitted))^2
      } else {
        r2 <- cor(observed, fitted)^2
      }
    }
    
    # Deviance pseudo-R2
    if (inherits(model, c("glm", "gam"))) {
      r2Dev <- (model$null.deviance - model$deviance) / model$null.deviance
    } else {
      r2Dev <- NA
    }
    
    # Negelkerke pseudo-R2
    if (termCount == 0) {
      logLikeInterceptOnly <- logLike
    }
    n <- length(observed)
    r2Negel <- (1 - exp((logLikeInterceptOnly - logLike) * (2 / n))) / (1 - exp(logLikeInterceptOnly * (2 / n)))
    
    # Add to summary dataframe
    obj$summary <- rbind(obj$summary, data.frame(
      term = term,
      k = length(coef(model)),
      logLike = logLike,
      aic = extractAIC(model)[2],
      r2 = r2,
      r2Dev = r2Dev,
      r2Negel = r2Negel
    ))
  }
  
  # R2 values presented as differences
  obj$summary <- within(obj$summary, {
    k <- c(1, diff(k))
    r2 <- c(NA, diff(r2))
    r2Dev <- c(NA, diff(r2Dev))
    r2Negel <- c(NA, diff(r2Negel))
  })
  
  # Calculate predicted values and handle interactions
  if (inherits(obj$model, "gam")) {
    # For GAMs, get predictions for each term separately
    preds_fit <- list()
    preds_se <- list()
    
    # Handle parametric terms first
    if (length(obj$model$pterms) > 0) {
      param_terms <- attr(obj$model$pterms, "term.labels")
      for (term in param_terms) {
        if (term != "") {  # Skip empty terms
          # Check if this is part of an interaction
          if (term %in% names(obj$interactions)) {
            # For interaction terms, we need a prediction grid
            pred_data <- create_prediction_grid(obj, term)
            pred <- predict(obj$model, newdata = pred_data, type = "terms", terms = term, se.fit = TRUE)
            
            # Map back to the original data points
            # This is a simplification - in practice we'd need to interpolate or match values
            all_preds <- predict(obj$model, type = "terms", terms = term, se.fit = TRUE)
            preds_fit[[term]] <- all_preds$fit
            preds_se[[term]] <- all_preds$se.fit
          } else {
            # For regular terms
            pred <- predict(obj$model, type = "terms", terms = term, se.fit = TRUE)
            preds_fit[[term]] <- pred$fit
            preds_se[[term]] <- pred$se.fit
          }
        }
      }
    }
    
    # Handle smooth terms
    for (i in seq_along(obj$model$smooth)) {
      term <- obj$model$smooth[[i]]$label
      
      # Check if this is an interaction term involving the focus or year
      if (term %in% names(obj$interactions) && 
          (obj$focus %in% obj$interactions[[term]]$vars || "year" %in% obj$interactions[[term]]$vars)) {
        # For interaction terms involving focus or year
        pred_data <- create_prediction_grid(obj, term)
        pred <- predict(obj$model, newdata = pred_data, type = "terms", terms = term, se.fit = TRUE)
        
        # Map back to original data - simplified approach
        all_preds <- predict(obj$model, type = "terms", terms = term, se.fit = TRUE)
        preds_fit[[term]] <- all_preds$fit
        preds_se[[term]] <- all_preds$se.fit
      } else {
        # For regular smooth terms
        pred <- predict(obj$model, type = "terms", terms = term, se.fit = TRUE)
        preds_fit[[term]] <- pred$fit
        preds_se[[term]] <- pred$se.fit
      }
    }
    
    # Handle spline smoothers from splines package if present in GAM
    for (term_name in names(obj$smoothers)) {
      if (obj$smoothers[[term_name]]$type %in% c("ns", "bs", "poly", "lo", "pspline", "rcs")) {
        # Create prediction grid for this spline
        pred_data <- create_prediction_grid(obj, term_name)
        
        # Get the variable and smoother type
        var_name <- obj$smoothers[[term_name]]$vars
        smoother_type <- obj$smoothers[[term_name]]$type
        
        # Extract parameters
        term_params <- ""
        param_match <- regexpr("\\([^)]*,([^)]+)\\)", term_name)
        if (param_match > 0) {
          match_text <- regmatches(term_name, param_match)
          param_part <- sub("^\\([^,]*,", "", match_text)
          param_part <- sub("\\)$", "", param_part)
          term_params <- paste0(", ", param_part)
        }
        
        # Create model formula for just this term
        formula_text <- paste0(obj$response, " ~ ", smoother_type, "(", var_name, term_params, ")")
        
        # Fit a temporary model to get predictions
        if (inherits(obj$model, "gam")) {
          temp_model <- gam(as.formula(formula_text), data = obj$data, family = obj$model$family)
        } else {
          temp_model <- glm(as.formula(formula_text), data = obj$data, family = obj$model$family)
        }
        
        # Get predictions
        preds <- predict(temp_model, newdata = pred_data, type = "terms", se.fit = TRUE)
        
        # Store predictions
        term_col <- paste0(smoother_type, "(", var_name)
        if (term_col %in% colnames(preds$fit)) {
          preds_fit[[term_name]] <- preds$fit[, term_col]
          preds_se[[term_name]] <- preds$se.fit[, term_col]
        }
      }
    }
    
    # Create dataframe with all predictions
    preds_df <- data.frame(obj$data)
    for (term in names(preds_fit)) {
      preds_df[[paste0("fit.", term)]] <- preds_fit[[term]]
      preds_df[[paste0("se.fit.", term)]] <- preds_se[[term]]
    }
    obj$preds <- preds_df
  } else {
    # For GLMs, handle interactions and splines specifically
    has_interactions <- length(obj$interactions) > 0
    has_splines <- any(sapply(obj$smoothers, function(s) s$type %in% c("ns", "bs", "poly", "lo", "pspline", "rcs")))
    
    if (has_interactions || has_splines) {
      # Get standard predictions for all terms
      preds <- predict(obj$model, type = "terms", se.fit = TRUE)
      fit <- as.data.frame(preds$fit)
      se.fit <- as.data.frame(preds$se.fit)
      
      # Process interactions involving the focus or year
      for (term in names(obj$interactions)) {
        interaction_info <- obj$interactions[[term]]
        if (obj$focus %in% interaction_info$vars || "year" %in% interaction_info$vars) {
          # For interactions involving focus or year, we need special handling
          if (term %in% colnames(fit)) {
            # The term exists in the predictions, no special handling needed
            next
          }
          
          # Create grid for this interaction
          pred_data <- create_prediction_grid(obj, term)
          
          # Get predictions on the grid
          pred_grid <- predict(obj$model, newdata = pred_data, type = "terms", se.fit = TRUE)
          
          # Map back to original data points - simplified approach
          if (term %in% colnames(pred_grid$fit)) {
            # Add as new column in the fit and se.fit dataframes
            grid_values <- pred_grid$fit[, term]
            grid_se <- pred_grid$se.fit[, term]
            
            # Map values to original data by matching focus or year levels
            mapped_values <- numeric(nrow(obj$data))
            mapped_se <- numeric(nrow(obj$data))
            
            # Determine which variable to use for mapping
            map_var <- if (obj$focus %in% interaction_info$vars) obj$focus else "year"
            
            for (level in levels(factor(obj$data[[map_var]]))) {
              level_rows <- obj$data[[map_var]] == level
              grid_rows <- pred_data[[map_var]] == level
              
              if (any(level_rows) && any(grid_rows)) {
                # Take the mean value for this level
                mapped_values[level_rows] <- mean(grid_values[grid_rows])
                mapped_se[level_rows] <- mean(grid_se[grid_rows])
              }
            }
            
            # Add to fit and se.fit
            fit[[term]] <- mapped_values
            se.fit[[term]] <- mapped_se
          }
        }
      }
      
      # Process spline smoothers
      for (term_name in names(obj$smoothers)) {
        smoother_info <- obj$smoothers[[term_name]]
        if (smoother_info$type %in% c("ns", "bs", "poly", "lo", "pspline", "rcs")) {
          var_name <- smoother_info$vars
          smoother_type <- smoother_info$type
          
          # Create necessary column name pattern to find in predictions
          pattern <- paste0("^", smoother_type, "\\(", var_name)
          
          # Find relevant columns in predictions
          term_columns <- grep(pattern, colnames(fit), value = TRUE)
          
          if (length(term_columns) > 0) {
            # Calculate combined effect for all basis functions of this smoother
            combined_fit <- rowSums(fit[, term_columns, drop = FALSE])
            
            # For standard errors, we need to account for covariance between basis functions
            # This is a simplification - full calculation would use vcov matrix
            combined_se <- sqrt(rowSums(se.fit[, term_columns, drop = FALSE]^2))
            
            # Add combined effect to predictions with the term name
            fit[[term_name]] <- combined_fit
            se.fit[[term_name]] <- combined_se
          } else {
            # If the spline term doesn't appear directly in predictions,
            # create a prediction grid and calculate manually
            pred_data <- create_prediction_grid(obj, term_name)
            
            # Extract parameters
            term_params <- ""
            param_match <- regexpr("\\([^)]*,([^)]+)\\)", term_name)
            if (param_match > 0) {
              match_text <- regmatches(term_name, param_match)
              param_part <- sub("^\\([^,]*,", "", match_text)
              param_part <- sub("\\)$", "", param_part)
              term_params <- paste0(", ", param_part)
            }
            
            # Create model formula for just this term
            formula_text <- paste0(obj$response, " ~ ", smoother_type, "(", var_name, term_params, ")")
            
            # Fit a temporary model
            temp_model <- glm(as.formula(formula_text), data = obj$data, family = obj$model$family)
            
            # Get predictions
            preds <- predict(temp_model, newdata = pred_data, type = "terms", se.fit = TRUE)
            
            # Map back to original data
            term_col <- paste0(smoother_type, "(", var_name)
            if (term_col %in% colnames(preds$fit)) {
              # Create mapping from variable levels to predictions
              var_levels <- pred_data[[var_name]]
              pred_values <- preds$fit[, term_col]
              pred_se <- preds$se.fit[, term_col]
              
              # For continuous variables, use approx to interpolate
              if (is.numeric(var_levels) && is.numeric(obj$data[[var_name]])) {
                mapped_values <- approx(var_levels, pred_values, obj$data[[var_name]], rule = 2)$y
                mapped_se <- approx(var_levels, pred_se, obj$data[[var_name]], rule = 2)$y
              } else {
                # For factors, match exactly
                mapped_values <- numeric(nrow(obj$data))
                mapped_se <- numeric(nrow(obj$data))
                
                for (level in unique(var_levels)) {
                  level_rows <- obj$data[[var_name]] == level
                  grid_rows <- var_levels == level
                  
                  if (any(level_rows) && any(grid_rows)) {
                    mapped_values[level_rows] <- mean(pred_values[grid_rows])
                    mapped_se[level_rows] <- mean(pred_se[grid_rows])
                  }
                }
              }
              
              # Add to fit and se.fit
              fit[[term_name]] <- mapped_values
              se.fit[[term_name]] <- mapped_se
            }
          }
        }
      }
      
      # Combine into preds dataframe
      preds <- cbind(fit, se.fit)
      names(preds) <- c(paste0("fit.", names(fit)), paste0("se.fit.", names(se.fit)))
      obj$preds <- cbind(obj$data, preds)
    } else {
      # Standard GLM without interactions or splines
      preds <- predict(obj$model, type = "terms", se.fit = TRUE)
      fit <- as.data.frame(preds$fit)
      se.fit <- as.data.frame(preds$se.fit)
      preds <- cbind(fit, se.fit)
      names(preds) <- c(paste0("fit.", names(fit)), paste0("se.fit.", names(se.fit)))
      obj$preds <- cbind(obj$data, preds)
    }
  }
  
  # Calculate influences and statistics for each term
  focus_var <- obj$data[[obj$focus]]
  if (!is.factor(focus_var)) {
    focus_var <- factor(focus_var)
  }
  
  obj$influences <- data.frame(level = levels(focus_var))
  overall <- c(NA, NA)  # NAs for null model and for focus term
  trend <- c(NA, NA)
  
  # Handle influences for all terms (including splines)
  for (term in obj$terms) {
    if (term != obj$focus) {
      # Check if this is an interaction term involving focus or year
      is_interaction <- term %in% names(obj$interactions) && 
                       (obj$focus %in% obj$interactions[[term]]$vars || "year" %in% obj$interactions[[term]]$vars)
      
      # Check if this is a spline smoother
      is_spline_smoother <- term %in% names(obj$smoothers) && 
                           obj$smoothers[[term]]$type %in% c("ns", "bs", "poly", "lo", "pspline", "rcs")
      
      if (is_interaction || is_spline_smoother) {
        # For interactions or splines, we need to use the prediction grid
        pred_data <- create_prediction_grid(obj, term)
        
        if (is_interaction) {
          # For interaction terms
          pred <- predict(obj$model, newdata = pred_data, type = "terms", terms = term)
        } else {
          # For spline terms, extract from model matrix and coefficients
          var_name <- obj$smoothers[[term]]$vars
          smoother_type <- obj$smoothers[[term]]$type
          
          # Extract parameters
          term_params <- ""
          param_match <- regexpr("\\([^)]*,([^)]+)\\)", term)
          if (param_match > 0) {
            match_text <- regmatches(term, param_match)
            param_part <- sub("^\\([^,]*,", "", match_text)
            param_part <- sub("\\)$", "", param_part)
            term_params <- paste0(", ", param_part)
          }
          
          # Build formula with just this term
          formula_text <- paste0("~", smoother_type, "(", var_name, term_params, ")")
          term_formula <- as.formula(formula_text)
          
          # Create model matrix
          mm <- model.matrix(term_formula, data = pred_data)
          
          # Get coefficients for this term
          coef_names <- grep(paste0("^", smoother_type, "\\(", var_name), names(coef(obj$model)), value = TRUE)
          term_coefs <- coef(obj$model)[coef_names]
          
          # Calculate fitted values for this term
          if (length(term_coefs) > 0) {
            pred <- mm[, -1, drop = FALSE] %*% term_coefs
          } else {
            pred <- rep(0, nrow(pred_data))
          }
        }
        
        # Determine which variable to use for aggregation (focus or year)
        agg_var <- if (obj$focus %in% obj$interactions[[term]]$vars) obj$focus else "year"
        
        # Aggregate by focus or year level
        infl <- aggregate(
          list(value = pred),
          by = list(level = pred_data[[obj$focus]]),
          mean
        )
        
        overall <- c(overall, with(infl, exp(mean(abs(value))) - 1))
        trend <- c(trend, with(infl, exp(cov(1:length(value), value) / var(1:length(value))) - 1))
        names(infl) <- c("level", term)
        obj$influences <- merge(obj$influences, infl, all.x = TRUE, by = "level")
      } else {
        # For regular terms, use standard approach
        term_col <- paste0("fit.", term)
        
        if (term_col %in% names(obj$preds)) {
          infl <- aggregate(
            list(value = obj$preds[[term_col]]),
            list(level = obj$preds[[obj$focus]]),
            mean
          )
          overall <- c(overall, with(infl, exp(mean(abs(value))) - 1))
          trend <- c(trend, with(infl, exp(cov(1:length(value), value) / var(1:length(value))) - 1))
          names(infl) <- c("level", term)
          obj$influences <- merge(obj$influences, infl, all.x = TRUE, by = "level")
        }
      }
    }
  }
  
  # Insert statistics into summary table
  obj$summary$overall <- overall
  obj$summary$trend <- trend
  
  return(obj)
}
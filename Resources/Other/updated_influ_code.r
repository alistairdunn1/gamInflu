# Updated initialize_influence function to handle splines package smoothers
initialize_influence <- function(obj) {
  # If no data was supplied
  if (is.null(obj$data)) {
    # See if the model has a data variable
    obj$data <- obj$model$data
    if (is.null(obj$data)) {
      # See if it is available as a model attribute
      obj$data <- attr(obj$model, "data")
    }
    if (is.null(obj$data)) {
      # For GAMs, use model$model
      if (inherits(obj$model, "gam")) {
        obj$data <- obj$model$model
      } else {
        stop("You need to provide the data that was fitted to for this type of model e.g. influence(model, data=mydata, ...)")
      }
    }
  }
  
  # Get term labels based on model type
  if (inherits(obj$model, "gam")) {
    # For GAMs, handle smooth terms differently
    parameterization <- obj$model$pterms
    obj$terms <- attr(parameterization, "term.labels")
    
    # Add smooth terms from the model
    smooth_terms <- names(obj$model$smooth)
    if (length(smooth_terms) > 0) {
      # Store interaction information
      obj$interactions <- list()
      # Store smoother information (including non-mgcv smoothers)
      obj$smoothers <- list()
      
      for (i in 1:length(obj$model$smooth)) {
        smooth_label <- obj$model$smooth[[i]]$label
        obj$terms <- c(obj$terms, smooth_label)
        
        # Check for interactions with year or focus
        if (inherits(obj$model$smooth[[i]], "tensor.smooth") || 
            inherits(obj$model$smooth[[i]], "tp") || 
            inherits(obj$model$smooth[[i]], "ti")) {
          
          # Get variables involved in this smooth
          vars <- obj$model$smooth[[i]]$term
          
          # Check if this is an interaction involving the focus term or year
          if (length(vars) > 1) {
            # Store interaction information
            obj$interactions[[smooth_label]] <- list(
              type = class(obj$model$smooth[[i]])[1],
              vars = vars,
              by = obj$model$smooth[[i]]$by,
              label = smooth_label
            )
          }
        }
        
        # Store smoother information for all smooth terms
        obj$smoothers[[smooth_label]] <- list(
          type = class(obj$model$smooth[[i]])[1],
          vars = obj$model$smooth[[i]]$term,
          label = smooth_label
        )
      }
    }
  } else {
    # For GLMs and other models
    obj$terms <- attr(obj$model$terms, "term.labels")
    
    # Initialize interaction and smoother tracking
    obj$interactions <- list()
    obj$smoothers <- list()
    
    # Check for non-GAM smoothers (ns, bs, etc.) and interactions in coefficient names
    coef_names <- names(coef(obj$model))
    
    # Function to extract smoother type and variable from coefficient names
    extract_smoother_info <- function(coef_name) {
      # Match smoothers like ns(x, df=4), bs(x, df=4), etc.
      # Both with and without numeric suffix
      smoother_match <- regexpr("^(ns|bs|poly|lo|pspline|rcs)\\(([^,]+)[^)]*\\)($|[0-9]+$)", coef_name)
      
      if (smoother_match > 0) {
        match_text <- regmatches(coef_name, smoother_match)
        
        # Extract smoother type (ns, bs, etc.)
        smoother_type <- gsub("\\(.*$", "", match_text)
        
        # Extract variable name
        var_match <- regexpr("\\(([^,]+)", match_text)
        var_name <- substr(regmatches(match_text, var_match), 2, 1000)
        
        return(list(
          type = smoother_type,
          var = var_name,
          full_name = gsub("([0-9]+$)", "", match_text)  # Remove any numeric suffix
        ))
      }
      
      return(NULL)
    }
    
    # Track detected smoothers
    found_smoothers <- list()
    
    # Check each coefficient name for smoothers
    for (coef_name in coef_names) {
      smoother_info <- extract_smoother_info(coef_name)
      
      if (!is.null(smoother_info)) {
        full_name <- smoother_info$full_name
        
        # Only add each unique smoother once
        if (!(full_name %in% names(found_smoothers))) {
          found_smoothers[[full_name]] <- smoother_info
        }
      }
    }
    
    # Add found smoothers to terms list and smoother info
    for (term_name in names(found_smoothers)) {
      smoother_info <- found_smoothers[[term_name]]
      
      # Add to terms list if not already there
      if (!(term_name %in% obj$terms)) {
        obj$terms <- c(obj$terms, term_name)
      }
      
      # Add smoother information
      obj$smoothers[[term_name]] <- list(
        type = smoother_info$type,
        vars = smoother_info$var,
        label = term_name
      )
    }
    
    # Check for interactions with year or focus in GLM terms
    for (term in obj$terms) {
      if (grepl(":", term, fixed = TRUE)) {
        # This is an interaction term
        parts <- strsplit(term, ":", fixed = TRUE)[[1]]
        obj$interactions[[term]] <- list(
          type = "parametric",
          vars = parts,
          label = term
        )
      }
    }
  }
  
  # Set response and focus if not provided
  if (is.null(obj$response)) {
    obj$response <- as.character(formula(obj$model)[[2]])
  }
  
  if (is.null(obj$focus)) {
    obj$focus <- obj$terms[1]
  }
  
  # Initialize labels and orders
  obj$labels <- setNames(as.list(obj$terms), obj$terms)
  obj$orders <- setNames(as.list(rep("asis", length(obj$terms))), obj$terms)
  
  return(obj)
}

# Updated get_coeffs function to handle non-mgcv smoothers
get_coeffs <- function(obj, model = obj$model, term = obj$focus) {
  UseMethod("get_coeffs", model)
}

# Updated get_coeffs.default to handle splines
get_coeffs.default <- function(obj, model = obj$model, term = obj$focus) {
  # Special handling for interaction terms
  if (term %in% names(obj$interactions)) {
    interaction_info <- obj$interactions[[term]]
    if ("year" %in% interaction_info$vars || obj$focus %in% interaction_info$vars) {
      # For interactions involving year or focus variable, use predicted effects
      pred <- get_term_effects(obj, model, term)
      return(pred)
    }
  }
  
  # Check if term is a smoother from splines package (ns, bs, etc.)
  if (term %in% names(obj$smoothers) && obj$smoothers[[term]]$type %in% c("ns", "bs", "poly", "lo", "pspline", "rcs")) {
    # For spline smoothers, use predicted effects
    pred <- get_term_effects(obj, model, term)
    return(pred)
  }
  
  # Standard handling for regular parametric terms
  coeffs <- coefficients(model)
  rows <- startsWith(names(coeffs), term)
  c(0, coeffs[rows])
}

# Updated get_coeffs.glm to handle splines
get_coeffs.glm <- function(obj, model = obj$model, term = obj$focus) {
  # Special handling for interaction terms
  if (term %in% names(obj$interactions)) {
    interaction_info <- obj$interactions[[term]]
    if ("year" %in% interaction_info$vars || obj$focus %in% interaction_info$vars) {
      # For interactions involving year or focus variable, use predicted effects
      pred <- get_term_effects(obj, model, term)
      return(pred)
    }
  }
  
  # Check if term is a smoother from splines package (ns, bs, etc.)
  if (term %in% names(obj$smoothers) && obj$smoothers[[term]]$type %in% c("ns", "bs", "poly", "lo", "pspline", "rcs")) {
    # For spline smoothers, use predicted effects
    pred <- get_term_effects(obj, model, term)
    return(pred)
  }
  
  # Standard handling for regular parametric terms
  coeffs <- coefficients(model)
  rows <- startsWith(names(coeffs), term)
  c(0, coeffs[rows])
}

# Updated get_coeffs.gam to handle splines
get_coeffs.gam <- function(obj, model = obj$model, term = obj$focus) {
  # Check if the term is a smooth term in GAM
  if (any(sapply(model$smooth, function(s) s$label == term))) {
    # For smooth terms, extract coefficients from predict
    pred <- get_term_effects(obj, model, term)
    return(pred)
  } else if (term %in% names(obj$smoothers) && obj$smoothers[[term]]$type %in% c("ns", "bs", "poly", "lo", "pspline", "rcs")) {
    # For spline smoothers from splines package, use predicted effects
    pred <- get_term_effects(obj, model, term)
    return(pred)
  } else {
    # For parametric terms
    coeffs <- coefficients(model)
    rows <- startsWith(names(coeffs), term)
    if (sum(rows) > 0) {
      c(0, coeffs[rows])
    } else {
      # Handle case where term might be part of an interaction but not a main effect
      pred <- get_term_effects(obj, model, term)
      return(pred)
    }
  }
}

# Updated get_term_effects to handle splines package smoothers
get_term_effects <- function(obj, model, term) {
  # Check if we're dealing with a focus term that's involved in interactions
  is_interaction <- term %in% names(obj$interactions)
  is_spline_smoother <- term %in% names(obj$smoothers) && obj$smoothers[[term]]$type %in% c("ns", "bs", "poly", "lo", "pspline", "rcs")
  has_interactions_with_focus <- any(sapply(obj$interactions, function(x) obj$focus %in% x$vars))
  has_interactions_with_year <- any(sapply(obj$interactions, function(x) "year" %in% x$vars))
  
  # Create a prediction grid
  if (is_interaction || is_spline_smoother || (term == obj$focus && has_interactions_with_focus) || 
      (term == "year" && has_interactions_with_year)) {
    # For interactions, splines, or focus terms with interactions
    pred_data <- create_prediction_grid(obj, term)
    
    # For splines package smoothers (ns, bs, etc.), we need special handling
    if (is_spline_smoother) {
      # Get the variable name and smoother type
      var_name <- obj$smoothers[[term]]$vars
      smoother_type <- obj$smoothers[[term]]$type
      
      # Create formula just for this term to get predictions
      # Extract any parameters from the term name (e.g., df from ns(x, df=4))
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
      
      # Create model matrix for just this term
      mm <- model.matrix(term_formula, data = pred_data)
      
      # Get coefficients for this term
      coef_names <- grep(paste0("^", smoother_type, "\\(", var_name), names(coef(model)), value = TRUE)
      term_coefs <- coef(model)[coef_names]
      
      # Calculate fitted values for this term
      term_effects <- mm[, -1, drop = FALSE] %*% term_coefs
      
      # Average by focus variable
      agg_effects <- aggregate(
        term_effects, 
        by = list(focus_level = pred_data[[obj$focus]]), 
        FUN = mean
      )
      
      effects <- agg_effects$x
    } else {
      # Get predictions for non-spline terms
      preds <- predict(model, newdata = pred_data, type = "terms", se.fit = TRUE)
      
      # Extract predictions for the term
      if (term %in% colnames(preds$fit)) {
        # Direct term effect
        effects <- preds$fit[, term]
      } else if (is_interaction) {
        # Interaction effect
        interaction_vars <- obj$interactions[[term]]$vars
        
        # For tensor smooths or ti() terms
        if (obj$interactions[[term]]$type %in% c("tensor.smooth", "tp", "ti")) {
          effects <- preds$fit[, term]
        } else {
          # For parametric interactions, we need to aggregate over the non-focus variables
          focus_var <- if(obj$focus %in% interaction_vars) obj$focus else 
                      if("year" %in% interaction_vars) "year" else interaction_vars[1]
          
          # Group by focus variable and average other effects
          agg_effects <- aggregate(
            preds$fit[, term], 
            by = list(focus_level = pred_data[[focus_var]]), 
            FUN = mean
          )
          effects <- agg_effects$x
        }
      } else {
        # No direct effect, might be marginal effect from interactions
        # Sum all effects related to the term
        related_terms <- colnames(preds$fit)[grepl(term, colnames(preds$fit))]
        if (length(related_terms) > 0) {
          effects_mat <- preds$fit[, related_terms, drop = FALSE]
          effects <- rowSums(effects_mat)
        } else {
          # Fallback
          effects <- rep(0, nrow(pred_data))
        }
      }
    }
    
    # Aggregate effects by focus variable
    if (term == obj$focus || (term == "year" && obj$focus != "year")) {
      # For the focus term or year, return the effects directly
      unique_effects <- unique(effects)
      return(c(0, unique_effects))
    } else {
      # For other terms, return the aggregated effects
      unique_effects <- unique(effects)
      return(c(0, unique_effects))
    }
  } else {
    # For non-interaction terms or terms not involving the focus or year
    if (inherits(model, "gam") && any(sapply(model$smooth, function(s) s$label == term))) {
      # For smooth terms
      newdata <- obj$data
      term_effect <- predict(model, newdata, type = "terms", terms = term)
      # Return unique levels with 0 prepended
      unique_levels <- unique(term_effect)
      return(c(0, unique_levels))
    } else {
      # For parametric terms
      coeffs <- coefficients(model)
      rows <- startsWith(names(coeffs), term)
      return(c(0, coeffs[rows]))
    }
  }
}

# Updated get_ses.default to handle splines
get_ses.default <- function(obj, model = obj$model, term = obj$focus) {
  # Check if term is a spline smoother
  is_spline_smoother <- term %in% names(obj$smoothers) && obj$smoothers[[term]]$type %in% c("ns", "bs", "poly", "lo", "pspline", "rcs")
  
  if (is_spline_smoother) {
    # For spline smoothers, we need to calculate SE based on the prediction
    pred_data <- create_prediction_grid(obj, term)
    
    # Get model matrix for this term
    var_name <- obj$smoothers[[term]]$vars
    smoother_type <- obj$smoothers[[term]]$type
    
    # Extract parameters from term name
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
    
    # Create model matrix for just this term
    mm <- model.matrix(term_formula, data = pred_data)[, -1, drop = FALSE]  # Remove intercept
    
    # Extract coefficient names for this term
    coef_names <- grep(paste0("^", smoother_type, "\\(", var_name), names(coef(model)), value = TRUE)
    
    # Extract relevant part of vcov matrix
    V <- vcov(model)[coef_names, coef_names, drop = FALSE]
    
    # Calculate SEs for predictions
    pred_ses <- sqrt(diag(mm %*% V %*% t(mm)))
    
    # Aggregate by focus level
    agg_ses <- aggregate(
      pred_ses,
      by = list(focus_level = pred_data[[obj$focus]]),
      FUN = mean
    )
    
    # Return appropriate number of SEs
    ses <- agg_ses$x
    return(ses)
  }
  
  # Standard handling for non-spline terms
  if (inherits(model, "survreg")) {
    V <- model$var
  } else {
    V <- vcov(model)
  }
  
  rows <- startsWith(row.names(V), term)
  if (sum(rows) > 0) {
    V <- V[rows, rows]
    n <- sum(rows) + 1
    Q <- matrix(-1 / n, nrow = n, ncol = n - 1)
    Q[-1, ] <- Q[-1, ] + diag(rep(1, n - 1))
    V0 <- (Q %*% V) %*% t(Q)
    se <- sqrt(diag(V0))
  } else {
    # Fallback for cases where we can't find the term directly
    se <- rep(mean(diag(V)), length(get_coeffs(obj, model, term)) - 1)
  }
  
  return(se)
}

# Updated get_ses.glm to handle splines
get_ses.glm <- function(obj, model = obj$model, term = obj$focus) {
  # Check if term is a spline smoother
  is_spline_smoother <- term %in% names(obj$smoothers) && obj$smoothers[[term]]$type %in% c("ns", "bs", "poly", "lo", "pspline", "rcs")
  
  # Special handling for interaction terms
  if (term %in% names(obj$interactions) || is_spline_smoother) {
    if (term %in% names(obj$interactions)) {
      interaction_info <- obj$interactions[[term]]
      if ("year" %in% interaction_info$vars || obj$focus %in% interaction_info$vars) {
        # For interactions involving year or focus variable, use predicted standard errors
        pred_data <- create_prediction_grid(obj, term)
        preds <- predict(model, newdata = pred_data, type = "terms", se.fit = TRUE)
        
        # Extract SEs for the term
        if (term %in% colnames(preds$se.fit)) {
          ses <- preds$se.fit[, term]
        } else {
          # For interactions, aggregate SEs
          ses <- rep(mean(diag(vcov(model))), length(get_coeffs(obj, model, term)) - 1)
        }
        
        # Return SEs
        return(ses)
      }
    }
    
    if (is_spline_smoother) {
      # Use spline-specific SE calculation
      return(get_ses.default(obj, model, term))
    }
  }
  
  # Standard handling for main effects
  V <- summary(model)$cov.scaled
  rows <- startsWith(row.names(V), term)
  V <- V[rows, rows]
  n <- sum(rows) + 1
  Q <- matrix(-1 / n, nrow = n, ncol = n - 1)
  Q[-1, ] <- Q[-1, ] + diag(rep(1, n - 1))
  V0 <- (Q %*% V) %*% t(Q)
  se <- sqrt(diag(V0))
  se
}

# Updated get_ses.gam to handle splines
get_ses.gam <- function(obj, model = obj$model, term = obj$focus) {
  # Check if term is a spline smoother from splines package
  is_spline_smoother <- term %in% names(obj$smoothers) && obj$smoothers[[term]]$type %in% c("ns", "bs", "poly", "lo", "pspline", "rcs")
  
  # Check if the term is a smooth term or involved in interactions
  if (any(sapply(model$smooth, function(s) s$label == term)) || term %in% names(obj$interactions) || 
      term == "year" || is_spline_smoother) {
    
    if (is_spline_smoother) {
      # Use spline-specific SE calculation
      return(get_ses.default(obj, model, term))
    }
    
    # For GAM smooth terms or interactions, extract standard errors from predict
    pred_data <- create_prediction_grid(obj, term)
    preds <- predict(model, newdata = pred_data, type = "terms", se.fit = TRUE)
    
    # Extract SEs for the term
    if (term %in% colnames(preds$se.fit)) {
      ses <- unique(preds$se.fit[, term])
    } else {
      # For terms not directly in the output (e.g., parts of interactions)
      # Use a conservative estimate based on the average SE across all terms
      all_ses <- as.vector(preds$se.fit)
      ses <- rep(mean(all_ses, na.rm = TRUE), length(get_coeffs(obj, model, term)) - 1)
    }
    
    return(ses)
  } else {
    # For parametric terms not in interactions
    V <- summary(model)$cov.scaled
    rows <- startsWith(row.names(V), term)
    
    if (sum(rows) > 0) {
      V <- V[rows, rows]
      n <- sum(rows) + 1
      Q <- matrix(-1 / n, nrow = n, ncol = n - 1)
      Q[-1, ] <- Q[-1, ] + diag(rep(1, n - 1))
      V0 <- (Q %*% V) %*% t(Q)
      se <- sqrt(diag(V0))
    } else {
      # Fallback for cases where we can't find the term directly
      se <- rep(mean(diag(V)), length(get_coeffs(obj, model, term)) - 1)
    }
    
    return(se)
  }
}

# Updated create_prediction_grid to handle splines
create_prediction_grid <- function(obj, term) {
  # Initialize an empty list to store the grid variables
  grid_vars <- list()
  
  # Always include the focus variable - ensuring factor handling
  focus_var <- obj$data[[obj$focus]]
  if (is.factor(focus_var)) {
    grid_vars[[obj$focus]] <- levels(focus_var)
  } else if (is.numeric(focus_var)) {
    grid_vars[[obj$focus]] <- seq(min(focus_var), max(focus_var), length.out = 20)
  } else {
    grid_vars[[obj$focus]] <- unique(focus_var)
  }
  
  # Check if term is a spline smoother from splines package
  is_spline_smoother <- term %in% names(obj$smoothers) && 
                        obj$smoothers[[term]]$type %in% c("ns", "bs", "poly", "lo", "pspline", "rcs")
  
  if (is_spline_smoother) {
    # For splines package smoothers, include the spline variable in the grid
    var_name <- obj$smoothers[[term]]$vars
    var_data <- obj$data[[var_name]]
    
    if (!(var_name %in% names(grid_vars))) {
      if (is.factor(var_data)) {
        grid_vars[[var_name]] <- levels(var_data)
      } else if (is.numeric(var_data)) {
        # For numeric variables used in splines, we need more points for a smooth curve
        grid_vars[[var_name]] <- seq(min(var_data), max(var_data), length.out = 50)
      } else {
        grid_vars[[var_name]] <- unique(var_data)
      }
    }
  }
  
  # If the term is an interaction, include all variables involved
  if (term %in% names(obj$interactions)) {
    interaction_vars <- obj$interactions[[term]]$vars
    for (var in interaction_vars) {
      if (var != obj$focus && !(var %in% names(grid_vars))) {
        var_data <- obj$data[[var]]
        if (is.factor(var_data)) {
          grid_vars[[var]] <- levels(var_data)
        } else if (is.numeric(var_data)) {
          grid_vars[[var]] <- seq(min(var_data), max(var_data), length.out = 10)
        } else {
          grid_vars[[var]] <- unique(var_data)
        }
      }
    }
  }
  
  # For all other terms involved in interactions with the focus or year
  if (term == obj$focus || term == "year") {
    for (interaction_name in names(obj$interactions)) {
      interaction_info <- obj$interactions[[interaction_name]]
      if (obj$focus %in% interaction_info$vars || "year" %in% interaction_info$vars) {
        for (var in interaction_info$vars) {
          if (var != obj$focus && var != "year" && !(var %in% names(grid_vars))) {
            var_data <- obj$data[[var]]
            if (is.factor(var_data)) {
              grid_vars[[var]] <- levels(var_data)
            } else if (is.numeric(var_data)) {
              grid_vars[[var]] <- seq(min(var_data), max(var_data), length.out = 10)
            } else {
              grid_vars[[var]] <- unique(var_data)
            }
          }
        }
      }
    }
  }
  
  # Handle year variable separately if it's not the focus
  if ("year" %in% names(obj$data) && !"year" %in% names(grid_vars) && "year" != obj$focus) {
    year_data <- obj$data[["year"]]
    if (is.factor(year_data)) {
      grid_vars[["year"]] <- levels(year_data)
    } else if (is.numeric(year_data)) {
      grid_vars[["year"]] <- seq(min(year_data), max(year_data), length.out = 10)
    } else {
      grid_vars[["year"]] <- unique(year_data)
    }
  }
  
  # Add all other variables needed for prediction
  all_vars <- names(obj$data)
  for (var in all_vars) {
    if (!(var %in% names(grid_vars)) && var != obj$response) {
      var_data <- obj$data[[var]]
      # Use median/mode values for variables not directly involved
      if (is.factor(var_data)) {
        # Use the most common level (mode)
        grid_vars[[var]] <- names(sort(table(var_data), decreasing = TRUE)[1])
      } else if (is.numeric(var_data)) {
        # Use median for numeric variables
        grid_vars[[var]] <- median(var_data)
      } else {
        # Use the first value for other types
        grid_vars[[var]] <- var_data[1]
      }
    }
  }
  
  # Create the grid using expand.grid
  pred_grid <- do.call(expand.grid, grid_vars)
  
  # Ensure all factors are preserved as factors with proper levels
  for (var in names(pred_grid)) {
    if (is.factor(obj$data[[var]])) {
      pred_grid[[var]] <- factor(pred_grid[[var]], levels = levels(obj$data[[var]]))
    }
  }
  
  return(pred_grid)
}

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
  
  # Add unstandardized index
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
  
  # Add standardized index
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
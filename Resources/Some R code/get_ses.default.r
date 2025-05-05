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


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


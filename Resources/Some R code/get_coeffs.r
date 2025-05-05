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


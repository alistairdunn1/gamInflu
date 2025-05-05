#' Extract standard errors for coefficients from a GLM
#'
#' @param obj The influence object
#' @param model The GLM model
#' @param term The term to extract SEs for
#' @return Standard errors for the coefficients
#' @export
get_ses.glm <- function(obj, model = obj$model, term = obj$focus) {
  # Special handling for interaction terms
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

#' Extract standard errors for coefficients from a GAM
#'
#' @param obj The influence object
#' @param model The GAM model
#' @param term The term to extract SEs for
#' @return Standard errors for the coefficients
#' @export
get_ses.gam <- function(obj, model = obj$model, term = obj$focus) {
  # Check if the term is a smooth term or involved in interactions
  if (any(sapply(model$smooth, function(s) s$label == term)) || term %in% names(obj$interactions)) {
    # For smooth terms or interactions, extract standard errors from predict
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

#' Default method for get_ses
#'
#' @param obj The influence object
#' @param model The model
#' @param term The term to extract SEs for
#' @return Standard errors for the coefficients
#' @export
get_ses.default <- function(obj, model = obj$model, term = obj$focus) {
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


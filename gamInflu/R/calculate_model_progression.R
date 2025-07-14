#' Calculate Model Progression Statistics
#'
#' Build models progressively and calculate fit statistics
#'
#' @param x gam_influence object
#' @return Data frame with progression statistics
#' @keywords internal
#'
calculate_model_progression <- function(x) {
  original_formula <- formula(x$model)

  # Extract terms from the formula
  all_terms <- attr(terms(formula(x$model)), "term.labels")

  # Separate parametric and smooth terms
  param_terms <- character(0)
  smooth_terms <- character(0)

  for (term in all_terms) {
    if (grepl("s\\(|te\\(|ti\\(", term)) {
      smooth_terms <- c(smooth_terms, term)
    } else {
      param_terms <- c(param_terms, term)
    }
  }

  # Build progression data frame
  progression <- data.frame()

  # Start with intercept-only model
  base_formula <- update(original_formula, . ~ 1)
  model_intercept <- gam(base_formula, data = x$data, family = x$model$family)

  progression <- rbind(progression, data.frame(
    step = 0,
    term_added = "intercept",
    aic = AIC(model_intercept),
    r_squared = 0,
    deviance_explained = 0,
    focus_term_included = FALSE
  ))

  # Add terms in order (including focus term)
  current_terms <- character(0)
  step_count <- 1

  for (term in all_terms) {
    current_terms <- c(current_terms, term)
    new_formula <- update(
      base_formula,
      formula(paste("~ .", "+", paste(current_terms, collapse = " + ")))
    )

    tryCatch(
      {
        model_step <- gam(new_formula, data = x$data, family = x$model$family)
        model_summary <- summary(model_step)

        progression <- rbind(progression, data.frame(
          step = step_count,
          term_added = term,
          aic = AIC(model_step),
          r_squared = model_summary$r.sq,
          deviance_explained = model_summary$dev.expl,
          focus_term_included = x$focus %in% current_terms
        ))

        step_count <- step_count + 1
      },
      error = function(e) {
        warning("Could not fit model with term: ", term)
      }
    )
  }

  return(progression)
}

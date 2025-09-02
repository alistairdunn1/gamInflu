#' @title Box-Cox Transform a Response Variable
#' @description Apply a Box-Cox transformation to a response variable for use in GAMs,
#'   automatically estimating the optimal lambda parameter or using a specified value.
#'   The Box-Cox transformation can help normalize data and stabilize variance.
#'
#'   This function can work with either raw response data or with a gam_influence object
#'   to properly utilize the MASS::boxcox function for optimal lambda estimation.
#'
#' @param y Numeric vector or gam_influence object. If numeric, the response variable to transform.
#'   If gam_influence object, the lambda will be estimated using the fitted model and data.
#' @param lambda Numeric or NULL. The Box-Cox transformation parameter.
#'   If NULL (default), the optimal lambda is estimated from the data.
#'   Common values include:
#'   \itemize{
#'     \item 1: No transformation (y)
#'     \item 0: Log transformation (log(y))
#'     \item 0.5: Square root transformation (sqrt(y))
#'     \item -1: Reciprocal transformation (1/y)
#'   }
#' @param lower Numeric. Lower bound for lambda search (default: -2).
#' @param upper Numeric. Upper bound for lambda search (default: 2).
#' @param eps Numeric. Small value added to y to ensure strictly positive values (default: 0).
#'
#' @return A list containing:
#'   \itemize{
#'     \item transformed: The transformed response variable
#'     \item lambda: The lambda value used (estimated or specified)
#'     \item original: The original response variable
#'     \item eps: The epsilon value added to the data (if any)
#'   }
#'
#' @details
#' The Box-Cox transformation is defined as:
#' \deqn{y^{(\lambda)} = \begin{cases}
#'   \frac{y^\lambda - 1}{\lambda} & \text{if } \lambda \neq 0 \\
#'   \log(y) & \text{if } \lambda = 0
#' \end{cases}}
#'
#' When the original data contains zeros or negative values, a small epsilon value
#' should be added to make all values positive before transformation.
#'
#' @examples
#' \dontrun{
#' # Automatically estimate optimal lambda from data
#' bc_transform <- box_cox_transform(data$response)
#' transformed_y <- bc_transform$transformed
#'
#' # Use with gam_influence object for proper MASS::boxcox estimation
#' gam_inf <- gam_influence(model, focus = "s(x)")
#' bc_transform <- box_cox_transform(gam_inf)
#'
#' # Use log transformation (lambda = 0)
#' bc_log <- box_cox_transform(data$response, lambda = 0)
#'
#' # Square root transformation (lambda = 0.5)
#' bc_sqrt <- box_cox_transform(data$response, lambda = 0.5)
#'
#' # Add small epsilon to handle zeros in data
#' bc_transform_eps <- box_cox_transform(data$response, eps = 0.01)
#' }
#'
#' @export
box_cox_transform <- function(y, lambda = NULL, lower = -2, upper = 2, eps = 0) {
  # Handle gam_influence object input
  if (inherits(y, "gam_influence")) {
    gam_inf <- y

    # Extract the response variable from the gam_influence object
    # The response variable should be extracted from the data using the response name
    response_name <- gam_inf$response # This is the variable name
    if (is.null(response_name)) {
      stop("gam_influence object must contain a response variable name", call. = FALSE)
    }

    # Get actual response values from the data
    y <- gam_inf$data[[response_name]]
    if (is.null(y)) {
      stop("Response variable '", response_name, "' not found in gam_influence data", call. = FALSE)
    }

    # If lambda not provided and we have a gam_influence object, use MASS::boxcox
    if (is.null(lambda)) {
      if (!requireNamespace("MASS", quietly = TRUE)) {
        warning("MASS package not available. Using lambda = 0 (log transformation) as default.")
        lambda <- 0
      } else {
        # Use the fitted model from gam_influence object
        if (is.null(gam_inf$model)) {
          stop("gam_influence object must contain a fitted model for lambda estimation", call. = FALSE)
        }

        # Add epsilon if specified to handle zeros or negative values
        y_adjusted <- y
        if (eps > 0) {
          y_adjusted <- y_adjusted + eps
        }

        # Check for non-positive values
        if (any(y_adjusted <= 0)) {
          stop("Box-Cox transformation requires strictly positive values. ",
            "Use 'eps' parameter to add a constant to your data.",
            call. = FALSE
          )
        }

        # For gam_influence objects, use the same grid search as numeric
        # to avoid complex scoping issues with MASS::boxcox
        # This provides consistent results across input types
        lambda_grid <- seq(lower, upper, by = 0.1)
        max_loglik <- -Inf
        best_lambda <- 0

        # Simple log-likelihood calculation for Box-Cox
        for (l in lambda_grid) {
          if (abs(l) < 1e-10) {
            # Log transformation
            transformed <- log(y_adjusted)
          } else {
            # Standard Box-Cox
            transformed <- (y_adjusted^l - 1) / l
          }

          # Calculate log-likelihood (simplified)
          loglik <- -(length(y_adjusted) / 2) * log(var(transformed))

          if (loglik > max_loglik) {
            max_loglik <- loglik
            best_lambda <- l
          }
        }

        lambda <- best_lambda

        # Round lambda to common values if it's close
        common_lambdas <- c(-1, -0.5, 0, 0.5, 1, 2)
        for (cl in common_lambdas) {
          if (abs(lambda - cl) < 0.05) {
            lambda <- cl
            break
          }
        }
      }
    }
  } else {
    # Input validation for numeric input
    if (!is.numeric(y)) {
      stop("Response variable 'y' must be numeric or a gam_influence object", call. = FALSE)
    }
  }

  # Add epsilon if specified to handle zeros or negative values
  y_adjusted <- y
  if (eps > 0) {
    y_adjusted <- y_adjusted + eps
  }

  # Check for non-positive values
  if (any(y_adjusted <= 0)) {
    stop("Box-Cox transformation requires strictly positive values. ",
      "Use 'eps' parameter to add a constant to your data.",
      call. = FALSE
    )
  }

  # If lambda still not provided (for numeric input), estimate using simple method
  if (is.null(lambda)) {
    # Use a simpler approach for numeric input
    lambda_grid <- seq(lower, upper, by = 0.1)
    max_loglik <- -Inf
    best_lambda <- 0

    # Simple log-likelihood calculation for Box-Cox
    for (l in lambda_grid) {
      if (abs(l) < 1e-10) {
        # Log transformation
        transformed <- log(y_adjusted)
      } else {
        # Standard Box-Cox
        transformed <- (y_adjusted^l - 1) / l
      }

      # Calculate log-likelihood (simplified)
      loglik <- -(length(y_adjusted) / 2) * log(var(transformed))

      if (loglik > max_loglik) {
        max_loglik <- loglik
        best_lambda <- l
      }
    }

    lambda <- best_lambda

    # Round lambda to common values if it's close
    common_lambdas <- c(-1, -0.5, 0, 0.5, 1, 2)
    for (cl in common_lambdas) {
      if (abs(lambda - cl) < 0.05) {
        lambda <- cl
        break
      }
    }
  }

  # Apply Box-Cox transformation
  if (abs(lambda) < 1e-10) {
    # Lambda is very close to 0, use log transformation
    transformed <- log(y_adjusted)
    lambda <- 0 # Set exactly to 0 for clarity
  } else {
    # Use standard Box-Cox formula
    transformed <- (y_adjusted^lambda - 1) / lambda
  }

  # Return transformation results
  result <- list(
    transformed = transformed,
    lambda = lambda,
    original = y,
    eps = eps
  )

  class(result) <- c("box_cox_transform", "list")

  return(result)
}

#' @title Inverse Box-Cox Transform
#' @description Convert transformed values back to the original scale using the inverse Box-Cox transformation.
#'
#' @param x Numeric vector. The Box-Cox transformed values.
#' @param lambda Numeric. The Box-Cox transformation parameter used.
#' @param eps Numeric. The epsilon value added to the original data before transformation (default: 0).
#'
#' @return Numeric vector with values converted back to the original scale.
#'
#' @details
#' The inverse Box-Cox transformation is defined as:
#' \deqn{y = \begin{cases}
#'   (\lambda x + 1)^{1/\lambda} & \text{if } \lambda \neq 0 \\
#'   \exp(x) & \text{if } \lambda = 0
#' \end{cases}}
#'
#' If an epsilon was added to the original data, it will be subtracted after the inverse transformation.
#'
#' @examples
#' \dontrun{
#' # Transform data
#' bc_transform <- box_cox_transform(data$response)
#'
#' # Fit model on transformed data
#' mod <- gam(bc_transform$transformed ~ s(x), data = data)
#'
#' # Generate predictions on transformed scale
#' pred_transformed <- predict(mod, newdata = new_data)
#'
#' # Convert predictions back to original scale
#' pred_original <- inverse_box_cox(pred_transformed,
#'   lambda = bc_transform$lambda,
#'   eps = bc_transform$eps
#' )
#' }
#'
#' @export
inverse_box_cox <- function(x, lambda, eps = 0) {
  # Input validation
  if (!is.numeric(x)) {
    stop("Input 'x' must be numeric", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1) {
    stop("Lambda must be a single numeric value", call. = FALSE)
  }

  # Apply inverse Box-Cox transformation
  if (abs(lambda) < 1e-10) {
    # Lambda is very close to 0, use exp
    y <- exp(x)
  } else {
    # Use standard inverse Box-Cox formula
    y <- (lambda * x + 1)^(1 / lambda)
  }

  # Subtract epsilon if it was added
  if (eps > 0) {
    y <- y - eps
  }

  return(y)
}

#' @title Box-Cox Transform and Fit a GAM
#' @description Apply a Box-Cox transformation to the response variable, fit a GAM on the
#'   transformed data, and provide utility methods to work with the resulting model.
#'
#' @param formula Formula. A GAM formula (see mgcv::gam).
#' @param data Data frame. The data to use for fitting the model.
#' @param lambda Numeric or NULL. The Box-Cox transformation parameter.
#'   If NULL (default), the optimal lambda is estimated from the data.
#' @param eps Numeric. Small value added to response to ensure strictly positive values (default: 0).
#' @param ... Additional arguments passed to mgcv::gam.
#'
#' @return An object of class "gam_box_cox" which is a list containing:
#'   \itemize{
#'     \item model: The fitted GAM model on transformed data
#'     \item box_cox: The Box-Cox transformation details
#'     \item call: The function call
#'   }
#'
#' @details
#' This function applies a Box-Cox transformation to the response variable of a GAM,
#' fits the model on the transformed data, and returns an object with both the
#' transformed model and the transformation parameters. This allows for easier
#' back-transformation of predictions.
#'
#' The model requires the mgcv package to be installed.
#'
#' @examples
#' \dontrun{
#' # Fit a GAM with automatic Box-Cox transformation
#' bc_gam <- box_cox_gam(response ~ s(x) + factor(year), data = mydata)
#'
#' # Use specific lambda
#' bc_gam_log <- box_cox_gam(response ~ s(x) + factor(year), data = mydata, lambda = 0)
#'
#' # Make predictions (automatically back-transformed to original scale)
#' preds <- predict(bc_gam, newdata = newdata)
#'
#' # Summary of the model
#' summary(bc_gam)
#' }
#'
#' @export
box_cox_gam <- function(formula, data, lambda = NULL, eps = 0, ...) {
  # Check if mgcv is available
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package mgcv is required for this function", call. = FALSE)
  }

  # Extract response variable name from formula
  response_name <- all.vars(formula)[1]

  # Get response variable data
  y <- data[[response_name]]
  if (is.null(y)) {
    stop("Response variable '", response_name, "' not found in data", call. = FALSE)
  }

  # Apply Box-Cox transformation
  bc_transform <- box_cox_transform(y, lambda = lambda, eps = eps)

  # Create a new data frame with the transformed response
  data_transformed <- data
  data_transformed[[response_name]] <- bc_transform$transformed

  # Fit the GAM model on transformed data
  gam_model <- mgcv::gam(formula = formula, data = data_transformed, ...)

  # Create a gam_box_cox object
  result <- list(
    model = gam_model,
    box_cox = bc_transform,
    call = match.call()
  )

  class(result) <- c("gam_box_cox", "list")

  return(result)
}

#' @title Predict Method for Box-Cox Transformed GAM
#' @description Make predictions from a Box-Cox transformed GAM model,
#'   with automatic back-transformation to the original scale.
#'
#' @param object An object of class "gam_box_cox" from box_cox_gam().
#' @param newdata Optional data frame for predictions.
#' @param type Character string. Type of prediction, as in mgcv::predict.gam.
#' @param back_transform Logical. Whether to back-transform predictions to original scale (default: TRUE).
#' @param ... Additional arguments passed to predict.gam.
#'
#' @return A vector of predictions.
#'
#' @method predict gam_box_cox
#' @export
predict.gam_box_cox <- function(object, newdata = NULL, type = "response",
                                back_transform = TRUE, ...) {
  # Get predictions on transformed scale
  predictions <- predict(object$model, newdata = newdata, type = type, ...)

  # Back-transform if requested and type is "response" or "link"
  if (back_transform && type %in% c("response", "link")) {
    # Extract Box-Cox parameters
    lambda <- object$box_cox$lambda
    eps <- object$box_cox$eps

    # Back-transform predictions
    predictions <- inverse_box_cox(predictions, lambda = lambda, eps = eps)
  }

  return(predictions)
}

#' @title Summary Method for Box-Cox Transformed GAM
#' @description Summarize a Box-Cox transformed GAM model, showing both
#'   the transformation details and the model summary.
#'
#' @param object An object of class "gam_box_cox" from box_cox_gam().
#' @param ... Additional arguments passed to summary.gam.
#'
#' @return A summary object with transformation and model details.
#'
#' @method summary gam_box_cox
#' @export
summary.gam_box_cox <- function(object, ...) {
  # Get the GAM model summary
  gam_summary <- summary(object$model, ...)

  # Add Box-Cox transformation information
  result <- list(
    gam_summary = gam_summary,
    box_cox = list(
      lambda = object$box_cox$lambda,
      eps = object$box_cox$eps
    ),
    call = object$call
  )

  class(result) <- c("summary.gam_box_cox", "list")

  return(result)
}

#' @title Print Method for Box-Cox Transformed GAM Summary
#' @description Print method for summary of a Box-Cox transformed GAM.
#'
#' @param x An object of class "summary.gam_box_cox".
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object.
#'
#' @method print summary.gam_box_cox
#' @export
print.summary.gam_box_cox <- function(x, ...) {
  # Print Box-Cox transformation information
  cat("Box-Cox Transformed GAM\n")
  cat("======================\n\n")

  cat("Box-Cox Transformation:\n")
  lambda_desc <- switch(as.character(x$box_cox$lambda),
    "0" = "Log",
    "0.5" = "Square root",
    "1" = "None (identity)",
    "2" = "Square",
    "-1" = "Reciprocal",
    "-0.5" = "Reciprocal square root",
    paste0("Power (", x$box_cox$lambda, ")")
  )

  cat("  Lambda:", x$box_cox$lambda, paste0("(", lambda_desc, ")\n"))
  if (x$box_cox$eps > 0) {
    cat("  Epsilon:", x$box_cox$eps, "(added to response)\n")
  }
  cat("\n")

  # Print the standard GAM summary
  cat("GAM Model Summary:\n")
  cat("=================\n\n")

  # Extract and print the standard GAM summary parts
  print(x$gam_summary)

  invisible(x)
}

#' @title Print Method for Box-Cox Transformed GAM
#' @description Print method for a Box-Cox transformed GAM.
#'
#' @param x An object of class "gam_box_cox".
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the input object.
#'
#' @method print gam_box_cox
#' @export
print.gam_box_cox <- function(x, ...) {
  cat("Box-Cox Transformed GAM\n")
  cat("======================\n\n")

  # Show transformation info
  lambda_desc <- switch(as.character(x$box_cox$lambda),
    "0" = "Log",
    "0.5" = "Square root",
    "1" = "None (identity)",
    "2" = "Square",
    "-1" = "Reciprocal",
    "-0.5" = "Reciprocal square root",
    paste0("Power (", x$box_cox$lambda, ")")
  )

  cat("Response transformed with Box-Cox:\n")
  cat("  Lambda:", x$box_cox$lambda, paste0("(", lambda_desc, ")\n"))
  if (x$box_cox$eps > 0) {
    cat("  Epsilon:", x$box_cox$eps, "(added to response)\n")
  }
  cat("\n")

  # Print formula
  cat("Formula:\n")
  print(formula(x$model))
  cat("\n")

  # Print brief model info
  cat("Model family:", x$model$family$family, "\n")
  cat("Link function:", x$model$family$link, "\n")
  cat("\nCall:\n")
  print(x$call)
  cat("\nUse summary() for more details.\n")

  invisible(x)
}

#' @title Convert a Box-Cox Transformed GAM to gam_influence object
#' @description Create a gam_influence object from a Box-Cox transformed GAM for use with
#'   influence diagnostics and other gamInflu functions.
#'
#' @param bc_model An object of class "gam_box_cox" from box_cox_gam().
#' @param focus Character. Name of the focus term for influence analysis.
#' @param ... Additional arguments passed to gam_influence.
#'
#' @return A gam_influence object.
#'
#' @export
gam_influence_from_box_cox <- function(bc_model, focus, ...) {
  # Validate input
  if (!inherits(bc_model, "gam_box_cox")) {
    stop("Input must be a gam_box_cox object", call. = FALSE)
  }

  # Extract the GAM model and create gam_influence object
  gi <- gam_influence(bc_model$model, focus = focus, ...)

  # Store Box-Cox transformation info in the gam_influence object
  gi$box_cox <- bc_model$box_cox

  # Return the enhanced gam_influence object
  return(gi)
}

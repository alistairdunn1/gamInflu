# A simplified Box-Cox transformation function
box_cox_simple <- function(y, lambda = NULL) {
  # Basic input validation
  if (!is.numeric(y)) {
    stop("Response variable 'y' must be numeric")
  }

  if (any(y <= 0)) {
    stop("Box-Cox transformation requires strictly positive values")
  }

  # If lambda not provided, estimate it with a simple approach
  if (is.null(lambda)) {
    # We'll use a simple grid search for the best lambda
    lambda_grid <- seq(-2, 2, by = 0.1)
    max_loglik <- -Inf
    best_lambda <- 0

    # Simple log-likelihood function for Box-Cox
    for (l in lambda_grid) {
      if (abs(l) < 1e-10) {
        # Log transformation
        transformed <- log(y)
      } else {
        # Standard Box-Cox
        transformed <- (y^l - 1) / l
      }

      # Calculate log-likelihood (simplified)
      loglik <- -(length(y) / 2) * log(var(transformed))

      if (loglik > max_loglik) {
        max_loglik <- loglik
        best_lambda <- l
      }
    }

    lambda <- best_lambda
  }

  # Apply the transformation
  if (abs(lambda) < 1e-10) {
    # Lambda is very close to 0, use log transformation
    transformed <- log(y)
    lambda <- 0 # Set exactly to 0 for clarity
  } else {
    # Use standard Box-Cox formula
    transformed <- (y^lambda - 1) / lambda
  }

  # Return results
  result <- list(
    transformed = transformed,
    lambda = lambda,
    original = y,
    eps = 0
  )

  return(result)
}

# Simple inverse Box-Cox function
inverse_box_cox_simple <- function(x, lambda) {
  if (abs(lambda) < 1e-10) {
    # Lambda is very close to 0, use exp
    y <- exp(x)
  } else {
    # Use standard inverse Box-Cox formula
    y <- (lambda * x + 1)^(1 / lambda)
  }

  return(y)
}

# Test with some basic data
set.seed(123)
y <- exp(rnorm(100)) # Generate some positive right-skewed data

# Test 1: Simple transformation with auto lambda
cat("Test 1: Box-Cox transformation with automatic lambda\n")
result <- box_cox_simple(y)
cat("Estimated lambda:", result$lambda, "\n")

# Test 2: Specific lambda values
cat("\nTest 2: Specific lambda values\n")
lambdas <- c(0, 0.5, 1, 2)
for (lambda in lambdas) {
  cat("Using lambda =", lambda, "\n")
  result <- box_cox_simple(y, lambda = lambda)
  cat("  - Mean of transformed data:", mean(result$transformed), "\n")
  cat("  - SD of transformed data:", sd(result$transformed), "\n")
}

# Test 3: Inverse transformation
cat("\nTest 3: Inverse transformation\n")
y_transform <- box_cox_simple(y, lambda = 0.5)$transformed
y_inverse <- inverse_box_cox_simple(y_transform, lambda = 0.5)
cat("Original mean:", mean(y), "\n")
cat("Inverse transform mean:", mean(y_inverse), "\n")
cat("Mean absolute difference:", mean(abs(y - y_inverse)), "\n")

cat("\nTest completed successfully\n")

# Comprehensive verification of Box-Cox transformation calculations
# This script tests the Box-Cox implementation against known theoretical values

# Source the Box-Cox functions
source("c:/Users/alist/OneDrive/Projects/Software/gamInflu/gamInflu/R/box_cox_transform.R")

# Function to manually calculate Box-Cox transformation for verification
manual_box_cox <- function(y, lambda) {
  if (abs(lambda) < 1e-10) {
    # Lambda = 0: log transformation
    return(log(y))
  } else {
    # Standard Box-Cox formula: (y^lambda - 1) / lambda
    return((y^lambda - 1) / lambda)
  }
}

# Function to manually calculate inverse Box-Cox transformation
manual_inverse_box_cox <- function(x, lambda) {
  if (abs(lambda) < 1e-10) {
    # Lambda = 0: exponential
    return(exp(x))
  } else {
    # Inverse Box-Cox formula: (lambda * x + 1)^(1/lambda)
    return((lambda * x + 1)^(1 / lambda))
  }
}

cat("=== Box-Cox Transformation Verification ===\n\n")

# Test data
test_values <- c(1, 2, 4, 8, 16)
test_lambdas <- c(-1, -0.5, 0, 0.5, 1, 2)

cat("1. Testing Forward Transformation\n")
cat("=================================\n")

for (lambda in test_lambdas) {
  cat(sprintf("Lambda = %g:\n", lambda))

  # Test with our function
  our_result <- box_cox_transform(test_values, lambda = lambda)
  our_transformed <- our_result$transformed

  # Test with manual calculation
  manual_transformed <- manual_box_cox(test_values, lambda)

  # Compare results
  max_diff <- max(abs(our_transformed - manual_transformed))

  cat("  Input values: ", paste(test_values, collapse = ", "), "\n")
  cat("  Our function: ", paste(round(our_transformed, 6), collapse = ", "), "\n")
  cat("  Manual calc:  ", paste(round(manual_transformed, 6), collapse = ", "), "\n")
  cat("  Max difference: ", sprintf("%.2e", max_diff), "\n")

  if (max_diff < 1e-12) {
    cat("  ✓ PASSED\n\n")
  } else {
    cat("  ✗ FAILED - Difference too large\n\n")
  }
}

cat("2. Testing Inverse Transformation\n")
cat("=================================\n")

for (lambda in test_lambdas) {
  cat(sprintf("Lambda = %g:\n", lambda))

  # First, transform the values
  transformed_values <- manual_box_cox(test_values, lambda)

  # Test our inverse function
  our_inverse <- inverse_box_cox(transformed_values, lambda = lambda)

  # Test manual inverse
  manual_inverse <- manual_inverse_box_cox(transformed_values, lambda)

  # Compare with original values
  max_diff_our <- max(abs(our_inverse - test_values))
  max_diff_manual <- max(abs(manual_inverse - test_values))

  cat("  Original values:     ", paste(test_values, collapse = ", "), "\n")
  cat("  Our inverse result:  ", paste(round(our_inverse, 6), collapse = ", "), "\n")
  cat("  Manual inverse:      ", paste(round(manual_inverse, 6), collapse = ", "), "\n")
  cat("  Our diff from orig:  ", sprintf("%.2e", max_diff_our), "\n")
  cat("  Manual diff from orig:", sprintf("%.2e", max_diff_manual), "\n")

  if (max_diff_our < 1e-12) {
    cat("  ✓ PASSED\n\n")
  } else {
    cat("  ✗ FAILED - Inverse transformation error too large\n\n")
  }
}

cat("3. Testing Special Cases\n")
cat("========================\n")

# Test lambda = 1 (should be identity minus 1)
cat("Lambda = 1 (should be y - 1):\n")
y_test <- c(1, 2, 3, 4, 5)
result_1 <- box_cox_transform(y_test, lambda = 1)
expected_1 <- y_test - 1
diff_1 <- max(abs(result_1$transformed - expected_1))
cat("  Input:     ", paste(y_test, collapse = ", "), "\n")
cat("  Result:    ", paste(result_1$transformed, collapse = ", "), "\n")
cat("  Expected:  ", paste(expected_1, collapse = ", "), "\n")
cat("  Max diff:  ", sprintf("%.2e", diff_1), "\n")
if (diff_1 < 1e-12) cat("  ✓ PASSED\n") else cat("  ✗ FAILED\n")

# Test lambda = 0 (should be log)
cat("\nLambda = 0 (should be log(y)):\n")
result_0 <- box_cox_transform(y_test, lambda = 0)
expected_0 <- log(y_test)
diff_0 <- max(abs(result_0$transformed - expected_0))
cat("  Input:     ", paste(y_test, collapse = ", "), "\n")
cat("  Result:    ", paste(round(result_0$transformed, 6), collapse = ", "), "\n")
cat("  Expected:  ", paste(round(expected_0, 6), collapse = ", "), "\n")
cat("  Max diff:  ", sprintf("%.2e", diff_0), "\n")
if (diff_0 < 1e-12) cat("  ✓ PASSED\n") else cat("  ✗ FAILED\n")

# Test lambda = 0.5 (should be 2*(sqrt(y) - 1))
cat("\nLambda = 0.5 (should be 2*(sqrt(y) - 1)):\n")
result_05 <- box_cox_transform(y_test, lambda = 0.5)
expected_05 <- 2 * (sqrt(y_test) - 1)
diff_05 <- max(abs(result_05$transformed - expected_05))
cat("  Input:     ", paste(y_test, collapse = ", "), "\n")
cat("  Result:    ", paste(round(result_05$transformed, 6), collapse = ", "), "\n")
cat("  Expected:  ", paste(round(expected_05, 6), collapse = ", "), "\n")
cat("  Max diff:  ", sprintf("%.2e", diff_05), "\n")
if (diff_05 < 1e-12) cat("  ✓ PASSED\n") else cat("  ✗ FAILED\n")

# Test lambda = -1 (should be -1 + 1/y = (1-y)/y)
cat("\nLambda = -1 (should be 1 - 1/y):\n")
result_neg1 <- box_cox_transform(y_test, lambda = -1)
expected_neg1 <- 1 - 1 / y_test
diff_neg1 <- max(abs(result_neg1$transformed - expected_neg1))
cat("  Input:     ", paste(y_test, collapse = ", "), "\n")
cat("  Result:    ", paste(round(result_neg1$transformed, 6), collapse = ", "), "\n")
cat("  Expected:  ", paste(round(expected_neg1, 6), collapse = ", "), "\n")
cat("  Max diff:  ", sprintf("%.2e", diff_neg1), "\n")
if (diff_neg1 < 1e-12) cat("  ✓ PASSED\n") else cat("  ✗ FAILED\n")

cat("\n4. Testing Round-Trip Transformation\n")
cat("====================================\n")

# Test that transform -> inverse gives back original values
set.seed(42)
test_data <- runif(10, min = 0.1, max = 10) # Random positive values

for (lambda in test_lambdas) {
  # Forward transformation
  transformed <- box_cox_transform(test_data, lambda = lambda)
  # Inverse transformation
  recovered <- inverse_box_cox(transformed$transformed, lambda = lambda)

  # Check if we get back the original
  max_error <- max(abs(recovered - test_data))

  cat(sprintf("Lambda = %g: Max round-trip error = %s", lambda, sprintf("%.2e", max_error)))
  if (max_error < 1e-12) {
    cat(" ✓ PASSED\n")
  } else {
    cat(" ✗ FAILED\n")
  }
}

cat("\n5. Testing Edge Cases\n")
cat("=====================\n")

# Test with epsilon
cat("Testing with epsilon parameter:\n")
y_with_zero <- c(0, 1, 2, 3)
eps_val <- 0.1
result_eps <- box_cox_transform(y_with_zero, lambda = 0.5, eps = eps_val)
# Manual calculation: transform (y + eps)
y_adjusted <- y_with_zero + eps_val
expected_eps <- manual_box_cox(y_adjusted, 0.5)
diff_eps <- max(abs(result_eps$transformed - expected_eps))

cat("  Input (with zeros):   ", paste(y_with_zero, collapse = ", "), "\n")
cat("  Epsilon:              ", eps_val, "\n")
cat("  Adjusted input:       ", paste(y_adjusted, collapse = ", "), "\n")
cat("  Result:               ", paste(round(result_eps$transformed, 6), collapse = ", "), "\n")
cat("  Expected:             ", paste(round(expected_eps, 6), collapse = ", "), "\n")
cat("  Max diff:             ", sprintf("%.2e", diff_eps), "\n")
if (diff_eps < 1e-12) cat("  ✓ PASSED\n") else cat("  ✗ FAILED\n")

# Test inverse with epsilon
cat("\nTesting inverse with epsilon:\n")
inverse_eps <- inverse_box_cox(result_eps$transformed, lambda = 0.5, eps = eps_val)
expected_inverse_eps <- y_with_zero
diff_inverse_eps <- max(abs(inverse_eps - expected_inverse_eps))

cat("  Inverse result:       ", paste(round(inverse_eps, 6), collapse = ", "), "\n")
cat("  Expected (original):  ", paste(expected_inverse_eps, collapse = ", "), "\n")
cat("  Max diff:             ", sprintf("%.2e", diff_inverse_eps), "\n")
if (diff_inverse_eps < 1e-12) cat("  ✓ PASSED\n") else cat("  ✗ FAILED\n")

cat("\n=== Verification Complete ===\n")
